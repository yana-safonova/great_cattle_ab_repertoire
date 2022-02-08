import os
import sys
import shutil
import pandas as pd

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

def ReadAndCombineIndividualDataframes(ind_dirs):
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
    return df

def ComputeGenePositionDataframe(gene, gene_df, num_f1_frac = 0.55):
    nucls = ['A', 'C', 'G', 'T']
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
        if min(pos_df[nucl1]) > num_f1_frac:
            continue
        # add to DF
        for i in range(len(pos_df)):
            selected_df['Pos'].append(gene + ':' + pos_id)
            selected_df['Fraction'].append(pos_df[nucl1][i] / (pos_df[nucl2][i] + pos_df[nucl1][i]))
            selected_df['Individual'].append(pos_df['Individual'][i])
    selected_df = pd.DataFrame(selected_df)
    return selected_df

def main(input_dir, output_dir):
    #### reading nucleotide frequencies
    ind_dirs = [f for f in os.listdir(input_dir) if f.find('alleles') != -1]
    df = ReadAndCombineIndividualDataframes(ind_dirs)
    #### preparing output dir
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)
    #### iterating over all V genes
    genes = set(df['Gene'])
    output_dfs = []
    for gene in genes:
        print('Processing ' + gene + '...')
        if gene not in ['IGHV1-7', 'IGHV1-10', 'IGHV1-14', 'IGHV1-17', 'IGHV1-20', 'IGHV1-21', 'IGHV1-27']:
            continue
        gene_df = df.loc[df['Gene'] == gene]
        selected_df = ComputeGenePositionDataframe(gene, gene_df) 
        output_fname = os.path.join(output_dir, 'df_' + gene + '.txt')
        OutputMatrix(selected_df, output_fname)
        output_dfs.append(pd.read_csv(output_fname, sep = '\t'))
    #### merge gene dataframes
    output_df = output_dfs[0]
    for df in output_dfs[1 :]:
        output_df = output_df.merge(df, on = 'Subject')
    output_df.to_csv(os.path.join(output_dir, 'merged_gsv.txt'), sep = '\t', index = False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('ERROR: invalid number of arguments')
        print('usage: python variation_dir output_dir')
        sys.exit(1)
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_dir, output_dir)
