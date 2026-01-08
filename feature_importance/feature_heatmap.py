import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--height', type=float, default=40) # 5.3
parser.add_argument('--width', type=float, default=20) # 17
parser.add_argument('--top_n', type=int, default=5)
parser.add_argument('--type', type=str)
parser.add_argument('--covarset', type=str)
args = parser.parse_args()

Cardiometabolic = ['MACE','Heart_Failure','CHD','PAD', 'Atrial_Fibrillation', 'T2_Diabetes']
Cancers = ['Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer']
eye = ['Cataracts', 'Glaucoma']
lung = ['COPD', 'Asthma']
dementia = ['Dementia']
fractures = ['Fractures']
metabolic = ['Renal_Disease', 'Liver_Disease']

feature_block_map = {}

def norm_top(df, trait, top_n):
    std = df['Final_Score'].std()
    df['norm'] = df['Final_Score']/ std
    df['abs_norm'] = df['norm'].abs()
    df.sort_values('abs_norm', ascending=False, inplace=True)
    if top_n > df[df['abs_norm'] > 1.96].shape[0]:
        df_top = df[df['abs_norm'] > 1.96]
    else:
        df_top = df.head(top_n).get(['Feature','norm'])
    df.rename(columns={'norm': f'norm_{trait}'}, inplace=True)
    df.drop('abs_norm', axis=1, inplace=True)
    return df_top['Feature'], df[['Feature', f'norm_{trait}']]

all_dfs = []
top_feats = None
for trait in globals().get(args.type):
    feat = pd.read_csv(f'mp_{trait}_{args.covarset}_weighted_abs.txt',sep="\t")
    feat = feat[feat['Feature'] != 'Batch'].reset_index().drop(columns=['index'])
    tmp_top_feat, all = norm_top(feat, trait, args.top_n)
    current_map = dict(zip(feat['Feature'], feat['Block']))
    feature_block_map.update(current_map)
    top_feats = tmp_top_feat if top_feats is None else pd.merge(top_feats,tmp_top_feat,on='Feature',how='outer')
    all_dfs.append(all)
for df in all_dfs:
    top_feats = pd.merge(top_feats,df,on='Feature',how='left')
print(top_feats)

def draw_single_heatmap(df,feat_map):
    all_values = df.drop('Feature',axis=1)
    vmin = all_values.min().min()
    vmax = all_values.max().max()

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    df['avg'] = df.iloc[:,1:].abs().mean(axis=1)
    df = df.sort_values('avg',ascending=False)
    df.drop('avg',axis=1,inplace=True)
    df.columns = df.columns.str.replace('norm_', '')
    df['Feature'] = df['Feature'].str.replace('_INVN', '')
    df_hmp = df.set_index('Feature').copy().T

    hm = sns.heatmap(df_hmp, ax = ax,vmin = vmin, vmax = vmax,
            cmap="bwr",  
            center = 0,
            annot=True,
            annot_kws={"fontsize":12},      
            fmt=".2f",       
            cbar=False,      
            cbar_kws={"shrink": 0.7}
        )
    ax.minorticks_off()
    ax.tick_params(axis='x', rotation=30, labelsize=24)
    color_dict = {'prote': "#060606", 'meta': "#40AC22"}
    clean_feat_map = {k.replace('_INVN', ''): v for k, v in feat_map.items()}
    for label in ax.get_xticklabels():
        feature_name = label.get_text()
        block_type = clean_feat_map.get(feature_name)
        if block_type in color_dict:
            label.set_color(color_dict[block_type])
        label.set_ha('right')         
        label.set_rotation_mode('anchor')
    ax.tick_params(axis='y', labelsize=24,rotation=0)
    ax.set_xlabel("")
    cbar_ax = fig.add_axes([0.93,0.4, 0.02, 0.3]) 
    cbar = fig.colorbar(hm.collections[0], cax=cbar_ax)
    cbar.ax.tick_params(labelsize=24,pad=10)
    cbar.ax.set_title('Normalized\nImportance\nScore', fontsize=18, pad=10)
    plt.gcf().subplots_adjust(left=0.15, right=0.90, bottom=0.2, top=0.9, hspace=0.3, wspace=0.3)
    plt.savefig(f'heatmap_phase2_{args.covarset}_{args.type}_sqrtP_abs_norm.pdf')
    plt.close()

draw_single_heatmap(top_feats,feature_block_map)


