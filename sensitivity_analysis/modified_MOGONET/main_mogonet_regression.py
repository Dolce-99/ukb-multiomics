from train_test_regression import train_test
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--trait",type=str)
parser.add_argument("--res_suffix",type=str)
parser.add_argument("--omic",type=str,required=False)
parser.add_argument('--large_cohort', action='store_true', required=False)
parser.add_argument('--num_epoch_pretrain', type=int, required=True)
parser.add_argument('--num_epoch', type=int, required=True)
parser.add_argument('--adj_parameter', type=int, required=False)
parser.add_argument('--random_state',type=int,required=False)
args = parser.parse_args()
trait = args.trait
res_suffix = args.res_suffix
omic = args.omic
large_cohort = args.large_cohort
rs = args.random_state

def load_data(res_suffix,trait,omic,large_cohort, rs):
    # read metabolies
    meta = pd.read_csv("ukb_metabofilter_complete_covars.txt",delimiter="\t")
    meta_field = pd.read_csv("metabolites.data.field.txt", delimiter = "\t")
    map_col_dict = pd.Series(meta_field.Abbreviation.values, index = meta_field['Field ID']).to_dict()

    def update_col_names(col_name):
        match = re.search(r'\d+', col_name)
        if match:
            number = int(match.group(0))
            name = map_col_dict.get(number, None)
            if name:
                return map_col_dict[number]
        return col_name

    meta.columns = [update_col_names(col) for col in meta.columns]
    meta_keep = meta.loc[:,['n_eid'] + [name for name in map_col_dict.values() if name in meta.columns]]
    meta_keep = meta_keep.rename(columns={'n_eid':'f.eid'})
    print("\nFinish reading metabolites")

    # read proteomics
    prote = pd.read_csv("proteomics_full.tab",sep="\t").rename(columns={'eid':'f.eid'})
    print("\nFinish reading proteomics")

    trait_res = trait + '_surv_res'
    trait_split = trait + '_split'
    if large_cohort:
        res = pd.read_csv(f"surv_res_{res_suffix}_all.tab", delimiter="\t",usecols=['f.eid',trait_res])
        idsplit = pd.read_csv(f"idsplit_full_rs{rs}.tab", delimiter="\t", usecols = ['f.eid',trait]).rename(columns={trait:trait_split})
    else:
        res = pd.read_csv(f"surv_res_{res_suffix}_overlap.tab", delimiter="\t",usecols=['f.eid',trait_res])
        idsplit = pd.read_csv(f"idsplit_overlap_rs{rs}.tab", delimiter="\t", usecols = ['f.eid',trait]).rename(columns={trait:trait_split})
    print(f"\nFinish reading residuals and idsplit")

    data_return = []
    if large_cohort and omic =="meta":
        # using full samples
        comb = pd.merge(meta_keep,res,on='f.eid',how='inner').dropna()
        labels = comb[['f.eid',trait_res]]
        comb_feid = comb[['f.eid']]
        meta_data = comb.drop(columns=[trait_res])
        idsplit = pd.merge(comb_feid,idsplit,on='f.eid',how='left')
        if labels.shape[0] == meta_data.shape[0]:
            data_return.append(meta_data)
            print(f"\n Using full metabolomics data, sample size: {meta_data.shape[0]}")
        else:
            raise ValueError(f"Error: inconsistent shape, labels:{labels.shape[0]}, data:{meta_data.shape[0]}")
    elif large_cohort and omic == "prote":
        comb = pd.merge(prote,res,on='f.eid',how='inner').dropna()
        labels = comb[['f.eid',trait_res]]
        comb_feid = comb[['f.eid']]
        prote_data = comb.drop(columns=[trait_res])
        idsplit = pd.merge(comb_feid,idsplit,on='f.eid',how='left')
        if labels.shape[0] == prote_data.shape[0]:
            data_return.append(prote_data)
            print("\n Using full proteomics data")
        else:
            raise ValueError(f"Error: inconsistent shape, labels:{labels.shape[0]}, data:{prote_data.shape[0]}")
    elif not large_cohort:
        # only using individuals that with both proteomics and metabolites
        meta_prote = pd.merge(meta_keep, prote,on='f.eid',how='inner')
        meta_prote_feid = meta_prote[['f.eid']]
        labels = pd.merge(meta_prote_feid,res,on='f.eid',how='inner').dropna() 
        comb_feid = labels[['f.eid']]
        idsplit = pd.merge(comb_feid,idsplit,on='f.eid',how='left')

        meta_data = pd.merge(comb_feid,meta_keep,on='f.eid',how='left')
        prote_data = pd.merge(comb_feid,prote,on='f.eid',how='left')
        if omic is None:
            data_return.append(meta_data)
            data_return.append(prote_data)
            print("\n Using overlaped individuals with meta+prote omics data")
        elif omic == "meta":
            data_return.append(meta_data)
            print("\n Using overlaped individuals with metabolomics only data")
        elif omic == "prote":
            data_return.append(prote_data)
            print("\n Using overlaped individuals with proteomics only data")


    return data_return, labels, comb_feid, idsplit


if __name__ == "__main__":
    if omic is None:
        view_list = [1,2]
    else:
        view_list = [1]
    num_epoch_pretrain = args.num_epoch_pretrain
    num_epoch = args.num_epoch
    lr_e_pretrain = 1e-3
    lr_e = 5e-4
    lr_c = 1e-3
    num_class=1
    
    data_return, labels, data_eid, idsplit = load_data(res_suffix,trait,omic,large_cohort, rs)
    print("\n Finish loading data")

    final_output = pd.DataFrame()
    
    trait_split = trait + '_split'
    trait_res = trait + '_surv_res'

    for fold in range(1,6):
        tr_eid = idsplit[idsplit[trait_split] != fold]['f.eid']
        te_eid = idsplit[idsplit[trait_split] == fold]['f.eid']
        labels_tr = pd.merge(tr_eid,labels,on='f.eid',how='left')[trait_res].to_numpy()
        labels_te = pd.merge(te_eid,labels,on='f.eid',how='left')[trait_res].to_numpy()
        data_list_tr, data_list_te = [],[]
        for data in data_return:
            tr_data = pd.merge(tr_eid,data,on='f.eid',how='left').drop(columns=['f.eid']).to_numpy()
            data_list_tr.append(tr_data)
            te_data = pd.merge(te_eid,data,on='f.eid',how='left').drop(columns=['f.eid']).to_numpy()
            data_list_te.append(te_data)

        pred = train_test(labels_tr,labels_te,data_list_tr,data_list_te, 
                view_list, num_class,
                lr_e_pretrain, lr_e, lr_c, 
                num_epoch_pretrain, num_epoch, args.adj_parameter, [400,400,200])
        print(pred)
        pred_df = pd.DataFrame(pred)
        pred_df.rename(columns={0:'pred'},inplace=True)
        pred_df = pred_df.reset_index(drop=True)
        te_eid = te_eid.reset_index(drop=True)
        output_df = pd.concat([te_eid,pred_df],axis=1).reset_index(drop=True)
        final_output = pd.concat([final_output,output_df],axis=0)
    
    if large_cohort and omic == "meta":
        omic_display = "full_meta"
    elif large_cohort and omic == "prote":
        omic_display = "full_prote"
    elif not large_cohort:
        if omic is None:
            omic_display = "mp"
        else:
            omic_display = "mp_" + omic

    final_output.to_csv(f"{omic_display}_{trait}_{res_suffix}_rs{rs}.txt", index=False,sep="\t")

         