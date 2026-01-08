import pandas as pd

Cardiometabolic = ['MACE','Heart_Failure','CHD','PAD', 'Atrial_Fibrillation', 'T2_Diabetes']
Cancers = ['Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer']
eye = ['Cataracts', 'Glaucoma']
lung = ['COPD', 'Asthma']
dementia = ['Dementia']
fractures = ['Fractures']
metabolic = ['Renal_Disease', 'Liver_Disease']

def make_tables_norm(trait,covarset):
    block_weights = pd.DataFrame()
    feat = pd.read_csv(f'/proj/yunligrp/users/djw/UKB_multiomics/revision1/02.new_features/mp_{trait}_{covarset}_weighted_abs.txt',sep="\t")
    feat = feat[feat['Feature'] != 'Batch'].reset_index().drop(columns=['index'])
    feat['norm_score'] = feat['Final_Score'] / feat['Final_Score'].std()
    meta_block_weight = feat[feat['Block']=='meta']['Block_Weight'].values[0]
    prot_block_weight = feat[feat['Block']=='prote']['Block_Weight'].values[0]
    weight_fold = prot_block_weight / meta_block_weight
    block_weight = {'metab': meta_block_weight, 'prot': prot_block_weight, 'weight_fold': weight_fold}
    block_weights = pd.DataFrame([block_weight])
    block_weights['Trait'] = trait
    block_weights['Predictor_set'] = covarset
    selected = feat[feat['norm_score'].abs() > 1.96]
    return_selected = selected[['Feature', 'Block','norm_score']].copy()
    return_selected['Feature'] = return_selected['Feature'].str.replace('_INVN','')
    return_selected['Trait'] = trait
    return_selected['Predictor_set'] = covarset
    return_selected['Omic'] = return_selected['Block'].replace({'meta': 'Metabolomics', 'prote': 'Proteomics'})
    return_selected.sort_values('norm_score', key=abs, ascending=False, inplace=True)
    return_selected = return_selected[['Trait', 'Predictor_set','Feature', 'Omic', 'norm_score']]
    return return_selected, block_weights

res_list = []
weight_list = []
for t in Cardiometabolic + Cancers + eye + lung + dementia + fractures + metabolic:
    tmp_df1, weight_df1 = make_tables_norm(trait=t, covarset='age_sex')
    res_list.append(tmp_df1)
    weight_list.append(weight_df1)
    tmp_df2, weight_df2 = make_tables_norm(trait=t, covarset='ASCVD')
    res_list.append(tmp_df2)
    weight_list.append(weight_df2)
    tmp_df3, weight_df3 = make_tables_norm(trait=t, covarset='PANEL')
    res_list.append(tmp_df3)
    weight_list.append(weight_df3)
block_weights_df = pd.concat(weight_list, axis=0, ignore_index=True)
res_df = pd.concat(res_list, axis=0, ignore_index=True)
res_df.to_csv('all_selected_features_norm.txt', sep="\t", index=False)


# save tables for single omics
def make_tables_norm_single_omics(trait,covarset,omic):
    feat = pd.read_csv(f'/proj/yunligrp/users/djw/UKB_multiomics/revision1/02.new_features/mp_{trait}_{covarset}_weighted_abs.txt',sep="\t")
    feat = feat[(feat['Feature'] != 'Batch') & (feat['Block'] == omic)].reset_index().drop(columns=['index'])
    feat['norm_score'] = feat['Final_Score'] / feat['Final_Score'].std()
    selected = feat[feat['norm_score'].abs() > 1.96]
    return_selected = selected[['Feature', 'norm_score']].copy()
    return_selected['Feature'] = return_selected['Feature'].str.replace('_INVN','')
    return_selected['Trait'] = trait
    return_selected['Predictor_set'] = covarset
    return_selected.sort_values('norm_score', key=abs, ascending=False, inplace=True)
    return_selected = return_selected[['Trait', 'Predictor_set','Feature', 'norm_score']]
    return return_selected

protein_list = []
metab_list = []
for t in Cardiometabolic + Cancers + eye + lung + dementia + fractures + metabolic:
    tmp_df1 = make_tables_norm_single_omics(trait=t, covarset='age_sex',omic='prote')
    protein_list.append(tmp_df1)
    tmp_df2 = make_tables_norm_single_omics(trait=t, covarset='ASCVD',omic='prote')
    protein_list.append(tmp_df2)
    tmp_df3 = make_tables_norm_single_omics(trait=t, covarset='PANEL',omic='prote')
    protein_list.append(tmp_df3)
    tmp_df4 = make_tables_norm_single_omics(trait=t, covarset='age_sex', omic='meta')
    metab_list.append(tmp_df4)
    tmp_df5 = make_tables_norm_single_omics(trait=t, covarset='ASCVD', omic='meta')
    metab_list.append(tmp_df5)
    tmp_df6 = make_tables_norm_single_omics(trait=t, covarset='PANEL', omic='meta')
    metab_list.append(tmp_df6)
protein_df = pd.concat(protein_list, axis=0, ignore_index=True)
protein_df.to_csv('all_selected_protein_features_norm.txt', sep="\t", index=False)
metab_df = pd.concat(metab_list, axis=0, ignore_index=True)
metab_df.to_csv('all_selected_metab_features_norm.txt', sep="\t", index=False)
