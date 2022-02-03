import os
import sys
import pandas as pd

import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns

def OutputMatrix(df, output_fname):
    subj_dict = dict()
    positions = sorted(set(df['Pos']))
    for i in range(len(df)):
        ind_id = df['Individual'][i]
        if ind_id not in subj_dict:
            subj_dict[ind_id] = dict()
        subj_dict[ind_id][df['Pos'][i]] = df['Fraction'][i]
    fh = open(output_fname, 'w')
    fh.write('Subject\t' + '\t'.join(positions) + '\n')
    for ind in sorted(subj_dict):
        fractions = []
        for pos in positions:
            f = 0
            if pos in subj_dict[ind]:
                f = subj_dict[ind][pos]
            fractions.append(str(f))
        fh.write(ind + '\t' + '\t'.join(fractions) + '\n')
    fh.close()    

def OutputPositionMatrix(pos_df, sorted_nucls, output_fname):
    matrix = []
    for i in range(len(pos_df)):
        row = []
        for n in sorted_nucls:
            row.append(pos_df[n][i])
        matrix.append(row)
    g = sns.clustermap(matrix, col_cluster = False, vmin = 0, vmax = 1, cmap = 'coolwarm', xticklabels = [str(n_ind + 1) + ':' + n for n_ind, n in enumerate(sorted_nucls)], yticklabels = [])
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 16)
    plt.savefig(output_fname)
    plt.close()
    plt.clf()

input_dir = sys.argv[1]
output_dir = sys.argv[2]

ind_dirs = [f for f in os.listdir(input_dir) if f.find('alleles') != -1]
dfs = []

nucl_cols = ['NumAs', 'NumCs', 'NumGs', 'NumTs']
for ind_dir in ind_dirs:
    ind_id = ind_dir.split('_')[0]
    full_path = os.path.join(input_dir, ind_dir)
    files = [f for f in os.listdir(full_path) if f.find('IGHV') != -1]
    for f in files:
        vgene = f.split('_')[0]
        file_path = os.path.join(full_path, f)
        df = pd.read_csv(file_path, sep = '\t')
        df['Gene'] = vgene
        df['Individual'] = ind_id
        df['NumSeqs'] = [(df['NumAs'][i] + df['NumCs'][i] + df['NumGs'][i] + df['NumTs'][i]) for i in range(len(df))]
        for col in nucl_cols:
            df[col[3]] = [float(df[col][i]) / df['NumSeqs'][i] for i in range(len(df))]
        dfs.append(df)

df = pd.concat(dfs)
genes = set(df['Gene'])
nucls = ['A', 'C', 'G', 'T']
output_dfs = []
for gene in genes:
    print(gene)
    if gene not in ['IGHV1-7', 'IGHV1-10', 'IGHV1-14', 'IGHV1-17', 'IGHV1-20', 'IGHV1-21', 'IGHV1-27']:
        continue
    gene_df = df.loc[df['Gene'] == gene]
    positions = sorted(set(gene_df['Position']))
    selected_df = {'Pos' : [], 'Fraction' : [], 'Individual' : []}
    for p in positions:
        pos_df = gene_df.loc[gene_df['Position'] == p].reset_index()
        nucl_dict = dict()
        for col in nucls:
            nucl_dict[col] = sum(pos_df[col])
        sorted_nucls = [n for n in sorted(nucl_dict, key = lambda x : nucl_dict[x], reverse = True)]
        nucl1 = sorted_nucls[0]
        nucl2 = sorted_nucls[1]
        pos_id = str(p) + '_' + nucl1 + nucl2
        if min(pos_df[nucl1]) > 0.55:
            continue
# report matrix
        OutputPositionMatrix(pos_df, sorted_nucls, os.path.join(output_dir, 'pos_' + gene + '_' + pos_id + '.pdf'))
# report ratios R 
        fractions = []
        for i in range(len(pos_df)):
            fractions.append(pos_df[nucl1][i] / (pos_df[nucl2][i] + pos_df[nucl1][i]))
        plt.figure()
        plt.hist(fractions)
        plt.xlim(0, 1.05)
        plt.xlabel('ratio $R$', fontsize = 14)
        plt.ylabel('# subjects', fontsize = 14)
        plt.savefig(os.path.join(output_dir, 'hist_' + gene + '_' + pos_id + '.pdf'))
        plt.clf()
# add to DF
        for i in range(len(pos_df)):
            selected_df['Pos'].append(gene + ':' + pos_id)
            selected_df['Fraction'].append(pos_df[nucl1][i] / (pos_df[nucl2][i] + pos_df[nucl1][i]))
            selected_df['Individual'].append(pos_df['Individual'][i])
    selected_df = pd.DataFrame(selected_df)
    print(set(selected_df['Pos']))
    output_fname = os.path.join(output_dir, 'df_' + gene + '.txt')
    OutputMatrix(selected_df, output_fname)
    output_dfs.append(pd.read_csv(output_fname, sep = '\t'))

output_df = output_dfs[0]
for df in output_dfs[1 :]:
    output_df = output_df.merge(df, on = 'Subject')
output_df.to_csv(os.path.join(output_dir, 'merged.txt'), sep = '\t', index = False)
