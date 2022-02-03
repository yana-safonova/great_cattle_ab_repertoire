import os
import sys
import pandas as pd

import numpy as np
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from kneed import KneeLocator

import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns

def ComputePValues(df, target_column, sep_column):
    states = list(set(df[target_column]))
    usage_dict = dict()
    for s in states:
        usage_dict[s] = []
    for i in range(len(df)):
        if not pd.isnull(df[sep_column][i]):
            usage_dict[df[target_column][i]].append(df[sep_column][i])
    pi_stats = ''
    if len(states) <= 1:
        return 'NA'
    if len(states) == 2:
        pi_stats = stats.f_oneway(usage_dict[states[0]], usage_dict[states[1]])
    elif len(states) == 3:
        pi_stats = stats.f_oneway(usage_dict[states[0]], usage_dict[states[1]], usage_dict[states[2]])
    else:
        pi_stats = ('', 'NA')
    return pi_stats[1]

def SelectLoadings(loading_df, max_num = 15):
    load_dict = dict()
    for i in range(len(loading_df)):
        load_dict[i] = abs(loading_df['PC1'][i]) * abs(loading_df['PC2'][i])
    selected_pos = []
    for pos_ind in sorted(load_dict, key = lambda x : load_dict[x], reverse = True):
        selected_pos.append(pos_ind)
        if len(selected_pos) == max_num:
            break
    return selected_pos

def ComputeAssociations(clusters, df, cols):
    selected_cols = []
    max_pval = 0.05 / len(cols)
    for c_ind, c in enumerate(cols):
        c_df = {'Cluster' : [], 'Frac' : []}
        for i in range(len(df)):
            c_df['Cluster'].append(clusters[i])
            c_df['Frac'].append(df[c][i])
        c_df = pd.DataFrame(c_df)
        pval = ComputePValues(c_df, 'Cluster', 'Frac')
        if pval == 'NA':
            continue
        if pval < max_pval:
            selected_cols.append(c_ind)
    return selected_cols

def OutputUsageMatrix(usage_df, labels, color_list, output_fname):
    print(labels)
    matrix = []
    subj_colors = []
    pos_columns = sorted([c for c in usage_df.columns if c not in ['Subject', 'Individual']])
    ind_labels = dict()
    for l_ind, l in enumerate(labels):
        ind_labels[l_ind] = l
    for ind in sorted(ind_labels, key = lambda x : ind_labels[x]):
        row = [usage_df[c][ind] for c in pos_columns]
        matrix.append(row)
        subj_colors.append(color_list[ind_labels[ind]])
    sns.clustermap(matrix, cmap = 'coolwarm', row_colors = subj_colors, yticklabels = [], xticklabels = [], row_cluster = False) #, vmin = 0, vmax = 1)
    plt.savefig(output_fname)
    plt.clf()
    plt.close()

usage_txt = sys.argv[1]
output_dir = sys.argv[2]
num_loadings = 100
if len(sys.argv) == 4:
    num_loadings = int(sys.argv[3])

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

usage_df = pd.read_csv(usage_txt, sep = '\t')
usage_lines = open(usage_txt).readlines()
matrix = []
subjects = []
vgenes = usage_lines[0].strip().split()[1 :]
for l in usage_lines[1 :]:
    splits = l.strip().split()
    matrix.append([float(s) for s in splits[1: ]])
    subjects.append(splits[0])

matrix = np.array(matrix)
matrix = StandardScaler().fit_transform(matrix)
pca = PCA(n_components = 2)
pca_result = pca.fit_transform(matrix)

pca_df = pd.DataFrame(pca_result)
loading_df = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'])
loading_df['Pos'] = vgenes
#print(loading_df)

loading_positions = SelectLoadings(loading_df, num_loadings)

colors = ['red', 'green', 'blue', 'yellow', 'orange', 'violet', 'black', 'pink', 'grey']
num_k = range(1, 10)
inertias = []
for k in num_k:
    print(k)
    model = KMeans(n_clusters = k)
    kmeans = model.fit(pca_df.iloc[:,:3])
    labels = kmeans.labels_
    subj_colors = [colors[l] for l in labels]
    plt.scatter(pca_df[0], pca_df[1], alpha=.5, color=subj_colors)
    OutputUsageMatrix(usage_df, labels, colors, os.path.join(output_dir, 'matrix_' + str(k) + '.pdf'))
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    print(labels)
    selected_pos_ind = ComputeAssociations(labels, usage_df, vgenes)
    print('Indices: ', [vgenes[ind] for ind in selected_pos_ind])
    for pos_ind in sorted(loading_positions, key = lambda x : loading_df['Pos'][x]):
        x = loading_df['PC1'][pos_ind]
        y = loading_df['PC2'][pos_ind]
        label = loading_df['Pos'][pos_ind].split('_')[0][4 :]
        if k == 1:
            print(label, x, y)
        plt.arrow(0, 0, x * 10, y * 10, width = 0.005, ec = '#8b8b8b', fc = '#8b8b8b')
#        plt.text(x * 2, y * 2, label)
    plt.savefig(os.path.join(output_dir, '_PCA_' + str(k) + '_loadings.pdf'))
    plt.clf()
    plt.scatter(pca_df[0], pca_df[1], alpha=.5, color=[colors[l] for l in labels])
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig(os.path.join(output_dir, '_PCA_' + str(k) + '.pdf'))
    plt.clf()    
    inertias.append(model.inertia_)
    fh = open(os.path.join(output_dir, 'PCA_' + str(k) + '.txt'), 'w')
    fh.write('Subject\tCluster\n')
    for ind, c_id in zip(subjects, labels):
        fh.write(ind + '\t' + str(c_id) + '\n')
    fh.close()

kneedle = KneeLocator(num_k, inertias, curve='convex', direction='decreasing')
print('elbow: ' + str(kneedle.elbow) + ', knee: ' + str(kneedle.knee))

plt.plot([kneedle.knee] * 2, [min(inertias), max(inertias)], color = 'b')
plt.plot(num_k, inertias, '-o', color='black')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(num_k)
plt.savefig(os.path.join(output_dir, 'PCA_interia.pdf'))
