import glob
import pandas as pd
import numpy as np
import os
from scipy.io import mmwrite
from scipy.sparse import coo_matrix, vstack, diags
import scipy as sp
import progressbar
import argparse

def ismember(a, b):
    """
    Function acts like matlab's ismember function. Returns indicies to map from a to b.
    """
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, None) for itm in a]


def main(score_dir):

    region_name = score_dir.split('/')[-1]
    print(f'''Processing {region_name}\n\t Reading unique snps''')
    # Read in SNPs
    snps = pd.read_csv(score_dir + '/unique.snps', header=None, sep='\t')
    snps.dropna(inplace=True)
    # For duplicate SNPs - for some reason uniq doesn't seem to collect
    if len(snps.drop_duplicates()) != len(snps):
        snps = snps.drop_duplicates()
        snps.to_csv(score_dir + '/unique.snps', index=False, header=False, sep='\t')
    score_files = glob.glob(score_dir + '/*.score')
    genes = [file.split('/')[-1].split('.')[0] for file in score_files]
    genes_copy = [file.split('/')[-1].split('.')[0] for file in score_files]
    # Create copy to allow items to be removed from gene list while iterating

    # Create string version of SNPs to allow for correct matching of multialleic variants
    snps_string = snps.apply(lambda x: ''.join(x), axis=1)

    # # Construct Sparse Weights Matrix
    col = []
    row = []
    vals = []
    row_ind = 0
    print('\t Iterating over gene score files')
    # For progress bar
    with progressbar.ProgressBar(max_value=len(genes)) as bar:
        for i, (gene, score_file) in enumerate(zip(genes_copy, score_files)):
            bar.update(i)
            # If file contains no lines
            try:
                score = pd.read_csv(score_file, sep='\t', header=None)
            except:
                genes.remove(gene)
                continue
            # Remove any NAs in last column that make weights get read in as strings
            if (score.loc[:,3].dtype != 'float64') & (score.loc[:,3].dtype != 'int64'):
                keep_ind = np.logical_not(score.loc[:, 3].str.contains('NA'))
                score = score.loc[keep_ind, :]
                score.loc[:, 3] = score.loc[:, 3].astype(float)

            # Remove entries that have TRUE in them - don't know why they are there
            keep_ind = ~((score==True).sum(axis=1).astype(bool))
            score = score.loc[keep_ind, :]
            score = score.dropna()

            # In case there is no lines in score file
            if len(score)==0:
                genes.remove(gene)
                continue

            # Append value to vals list
            vals = vals + score.loc[:, 3].values.tolist()
            # Convert score to string for index matching
            score_string = score.iloc[:, :3].apply(lambda x: ''.join(x), axis=1)
            # ind = ismember(score_string, snps_string)
            ind = np.where(score_string.isin(snps_string))[0].tolist()
            # Append row and col indicies
            col = col + ind
            row = row + (np.ones(len(ind))*row_ind).tolist()
            row_ind += 1

    print('\t Converting weights to sparse matrix')
    # Convert rows and columns to arrays
    row = np.array(row)
    col = np.array(col)
    vals = np.array(vals)
    W_TWAS = coo_matrix((vals, (row, col)), shape=(row_ind, len(snps)))

    # Remove all zero SNPs (columns)
    C = sp.sparse.find(W_TWAS)[1]
    keep_snps = (np.in1d(np.arange(W_TWAS.shape[1]), C))
    W_TWAS = W_TWAS.tocsr()[:, keep_snps].tocoo()
    snps = snps.iloc[keep_snps,:]

    out_dir = score_dir 
    print(f'''Saving as .mm file in {out_dir}''')
    # if not (os.path.isdir(out_dir)):
    #     os.makedirs(out_dir)
    # mmwrite(out_dir+'/W_TWAS.mm', W_TWAS)
    sp.sparse.save_npz(out_dir + '/W_TWAS.npz', W_TWAS, 'csr')
    pd.Series(genes).to_csv(out_dir + '/gene_names.txt', index=False, header=False)
    snps.to_csv(out_dir + '/snps.txt', index=False, header=False, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Function converts a series of TWAS score files to a scipy sparse matrix')
    parser.add_argument('-score_dir', help='Directory which contains score files', type=str, default=None)

    opt = parser.parse_args()

    main(opt.score_dir)

    