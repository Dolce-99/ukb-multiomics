""" Training and testing of the model
"""
import numpy as np
import torch
import torch.nn.functional as F
from models import init_model_dict, init_optim
from utils import cal_sample_weight, gen_adj_mat_tensor, gen_test_adj_mat_tensor, cal_adj_mat_parameter

cuda = True if torch.cuda.is_available() else False


def prepare_trte_data(labels_tr,labels_te, data_tr_list, data_te_list, view_list):
    num_view = len(view_list)
    num_tr = data_tr_list[0].shape[0]
    num_te = data_te_list[0].shape[0]
    data_mat_list = []
    for i in range(num_view):
        data_mat_list.append(np.concatenate((data_tr_list[i], data_te_list[i]), axis=0))
    data_tensor_list = []
    for i in range(len(data_mat_list)):
        data_tensor_list.append(torch.FloatTensor(data_mat_list[i]))
        if cuda:
            data_tensor_list[i] = data_tensor_list[i].cuda()
    idx_dict = {}
    idx_dict["tr"] = list(range(num_tr))
    idx_dict["te"] = list(range(num_tr, (num_tr+num_te)))
    data_train_list = []
    data_all_list = []
    for i in range(len(data_tensor_list)):
        data_train_list.append(data_tensor_list[i][idx_dict["tr"]].clone())
        data_all_list.append(torch.cat((data_tensor_list[i][idx_dict["tr"]].clone(),
                                       data_tensor_list[i][idx_dict["te"]].clone()),0))
    labels = np.concatenate((labels_tr, labels_te))
    
    return data_train_list, data_all_list, idx_dict, labels


def gen_trte_adj_mat(data_tr_list, data_trte_list, trte_idx, adj_parameter):
    adj_metric = "cosine" # cosine distance
    adj_train_list = []
    adj_test_list = []
    for i in range(len(data_tr_list)):
        adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, data_tr_list[i], adj_metric)
        adj_train_list.append(gen_adj_mat_tensor(data_tr_list[i], adj_parameter_adaptive, adj_metric))
        adj_test_list.append(gen_test_adj_mat_tensor(data_trte_list[i], trte_idx, adj_parameter_adaptive, adj_metric))
    
    return adj_train_list, adj_test_list


def train_epoch(data_list, adj_list, label, sample_weight, model_dict, optim_dict, train_VCDN=True):
    loss_dict = {}
    total_loss = 0
    criterion = torch.nn.MSELoss(reduction='none')
    for m in model_dict:
        model_dict[m].train()    
    num_view = len(data_list)
    for i in range(num_view):
        optim_dict["C{:}".format(i+1)].zero_grad()
        ci_loss = 0
        ci = model_dict["C{:}".format(i+1)](model_dict["E{:}".format(i+1)](data_list[i],adj_list[i]))
        ci_loss = torch.mean(torch.mul(criterion(ci, label),sample_weight))
        ci_loss.backward()
        optim_dict["C{:}".format(i+1)].step()
        loss_dict["C{:}".format(i+1)] = ci_loss.detach().cpu().numpy().item()
        total_loss += ci_loss.detach().cpu().numpy().item()
    if train_VCDN and num_view >= 2:
        optim_dict["C"].zero_grad()
        c_loss = 0
        ci_list = []
        for i in range(num_view):
            ci_list.append(model_dict["C{:}".format(i+1)](model_dict["E{:}".format(i+1)](data_list[i],adj_list[i])))
        c = model_dict["C"](ci_list)    
        c_loss = torch.mean(torch.mul(criterion(c, label),sample_weight))
        c_loss.backward()
        optim_dict["C"].step()
        loss_dict["C"] = c_loss.detach().cpu().numpy().item()
        total_loss += c_loss.detach().cpu().numpy().item()

    return loss_dict, total_loss
    

def test_epoch(data_list, adj_list, te_idx, model_dict):
    for m in model_dict:
        model_dict[m].eval()
    num_view = len(data_list)
    ci_list = []
    for i in range(num_view):
        ci_list.append(model_dict["C{:}".format(i+1)](model_dict["E{:}".format(i+1)](data_list[i],adj_list[i])))
    if num_view >= 2:
        c = model_dict["C"](ci_list)    
    else:
        c = ci_list[0]
    c = c[te_idx,:]
    prob = c.data.cpu().numpy()
    return prob


def train_test(labels_tr, labels_te, data_list_tr, data_list_te, 
               view_list, num_class,
               lr_e_pretrain, lr_e, lr_c, 
               num_epoch_pretrain, num_epoch,adj_parameter, dim_he_list):
    test_inverval = 10
    num_view = len(view_list)
    dim_hvcdn = 16
    adj_parameter = adj_parameter
    dim_he_list = dim_he_list # number of hidden units in each layer
    data_tr_list, data_trte_list, trte_idx, labels_trte = prepare_trte_data(labels_tr, labels_te, data_list_tr, data_list_te, view_list)
    labels_tr_tensor = torch.FloatTensor(labels_trte[trte_idx["tr"]])
    if len(labels_tr_tensor.shape) == 1:
        labels_tr_tensor = labels_tr_tensor.unsqueeze(1)
    sample_weight_tr = cal_sample_weight(labels_trte[trte_idx["tr"]], num_class)
    sample_weight_tr = torch.FloatTensor(sample_weight_tr)
    labels_tr_tensor = labels_tr_tensor.view(-1, 1)
    sample_weight_tr = sample_weight_tr.view(-1, 1)
    if cuda:
        labels_tr_tensor = labels_tr_tensor.cuda()
        sample_weight_tr = sample_weight_tr.cuda()
    adj_tr_list, adj_te_list = gen_trte_adj_mat(data_tr_list, data_trte_list, trte_idx, adj_parameter)
    dim_list = [x.shape[1] for x in data_tr_list]
    model_dict = init_model_dict(num_view, num_class, dim_list, dim_he_list, dim_hvcdn)
    for m in model_dict:
        if cuda:
            model_dict[m].cuda()
    
    print("\nPretrain GCNs...")
    optim_dict = init_optim(num_view, model_dict, lr_e_pretrain, lr_c)
    for epoch in range(num_epoch_pretrain):
        loss_dict, total_loss = train_epoch(data_tr_list, adj_tr_list, labels_tr_tensor, sample_weight_tr,
                    model_dict, optim_dict, train_VCDN=False)
        print(f"Pretrain Epoch: {epoch}, Total Loss: {total_loss}")

    print("\nTraining...")
    optim_dict = init_optim(num_view, model_dict, lr_e, lr_c)
    previous_loss = None
    for epoch in range(num_epoch+1):
        loss_dict, total_loss = train_epoch(data_tr_list, adj_tr_list, labels_tr_tensor, sample_weight_tr, model_dict, optim_dict)
        if previous_loss is not None and abs(total_loss - previous_loss) < 1e-6:
            print(f"Early stopping at epoch {epoch}")
            break
        previous_loss = total_loss
        if epoch % test_inverval == 0:
            print(f"\nTest: Epoch {epoch}, Total Loss: {total_loss}")
    te_prob = test_epoch(data_trte_list, adj_te_list, trte_idx["te"], model_dict)
    return te_prob