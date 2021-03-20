import os
import argparse
import pandas as pd
import numpy as np
import scipy.sparse
from pandas_plink import read_plink1_bin, read_plink

def __gene_reduce(gene_list, gene_map):
    """" 
    Private fuction for use in gene_map, reduces gene_list and gene_map intersection of genes. 

    Parameters
    ----------
    gene_list: pd.Series
    gene_map: pd.DataFrame  in biomart file
    
    """

    # Find intersection
    intersect_genes = list(set(gene_list) & set(gene_map.index))
    gene_list_ind = gene_list.isin(intersect_genes)
    # Raise Error if there seems to be a poor mapping to gene list
    if (sum(gene_list_ind)/len(gene_list_ind))<0.99:
        raise Exception('Some genes are not found in mapping - make sure that you are using the correct g_build')
    return(gene_list_ind, gene_map.loc[intersect_genes, :])


def gene_map(gene_list, g_build, multimatch_filt=None):
    """
    Function takes gene list and remaps, if input list is ensble it will map to gene names otherwise it will do the inverse. 
    
    Parameters
    ----------
    gene_list : list like
        A list of genes to be remapped.
    g_build : int
        Value of 37 or 38 depending on the ensmble genome build   
    multimatch_filt : dict:
        Dictionary for how to handle many to one matching
    Returns
    ----------
    tuple
        element_1 - a boolian list indicating elements of gene_list present in genome build.
        element_2 - remapped values
    """
    gene_list = pd.Series(gene_list)
    # Decide if mapping is ensmble to gene names or the other way around
    if gene_list.str.startswith('ENSG').sum() == len(gene_list):
        input_gene_type = 'ensbmle'
        index_col = 'Gene stable ID'
        new_col = 'Gene name'
        split_dot = lambda x: x.split('.')[0] # Function to get rid of version number
        gene_list = gene_list.apply(split_dot)
    else:
        input_gene_type = 'names'
        index_col = 'Gene name'
        new_col = 'Gene stable ID'
        
    gene_map_file = f'''{os.path.dirname(os.path.realpath(__file__))}/Gene_Ensmbl_Map/GRCh{g_build}.p13.biomart.tsv'''
    gene_map = pd.read_csv(gene_map_file, sep='\t', index_col=index_col)
    gene_map = gene_map.loc[~gene_map.index.duplicated(keep='first')]  # drop rows with duplicated indices
    
    # Match genes
    gene_list_ind, gene_map = __gene_reduce(gene_list, gene_map)

    # If there are duplicates and multimatch_filt is defined
    if (multimatch_filt != None) and (gene_map[new_col].duplicated().sum()!=0):
        if 'gene_version' in multimatch_filt.keys():
            if multimatch_filt['gene_version'] == 'newest':
                ancending=False
            elif multimatch_filt['gene_version'] == 'oldest':
                ancending=True
            dup_ind = (gene_map[new_col].duplicated(keep=False)).values
            dup_genes = gene_map.iloc[dup_ind, :]
            dup_genes = dup_genes.sort_values('Version (gene)', ascending=ancending)
            remove_genes = dup_genes.drop_duplicates(subset=new_col, keep='last').index
            gene_map = gene_map.drop(index=remove_genes)
            gene_list_ind, gene_map = __gene_reduce(gene_list, gene_map)
            
        # Transcript filter e.g. protein_coding
        if 'transcript_filter' in multimatch_filt.keys() and (gene_map[new_col].duplicated().sum()!=0):
            # Make transcript type list
            trascript_filter = multimatch_filt['transcript_filter']
            if isinstance(trascript_filter, str):
                trascript_filter = [trascript_filter]
                
            dup_ind = (gene_map[new_col].duplicated(keep=False)).values
            dup_ind_good = gene_map.loc[dup_ind, 'Transcript type'].isin(trascript_filter)
            remove_genes = dup_ind_good.index[~dup_ind_good]
            gene_map = gene_map.drop(index=remove_genes)
            
            # Repeat lines above
            gene_list_ind, gene_map = __gene_reduce(gene_list, gene_map)
            
    reorder_i = np.where(gene_list[gene_list_ind].isin(gene_map.index))[0]
    return (gene_list_ind, gene_map.iloc[reorder_i, :][new_col])


def TWAS_project(genetics, twas_sparse_dir, output_file, g_build):
    """
    This function takes individual level genetics plink files (.bed) and projects it using the 
    TWAS weights defined by using sparse W_TWAS.npz, gene_names.txt and snps.txt files in twas_sparse_dir. 
    These weight files are created by score_to_sparse.py in S1_create_sparse_weights.sh.
    
    Parameters
    ----------
    genetics : str
        Defines path to plink genetics file (e.g. /path/to/all_chroms.bed).
    twas_sparse_dir : str
        Path to directory that contains W_TWAS.npz, gene_names.txt and snps.txt files, created by score_to_sparse.py
    output_file : str
        Path to out file for saving TWAS projection, should end in '.tsv' or '.hdf'
    g_build : int
        Value of 37 or 38 depending on the ensmble genome build .
    Returns
    ----------
    None
    """
    precision=np.float32
    log_file = output_file.split('.')[0] + '.log'
    if os.path.isfile(log_file): os.remove(log_file)

    with open(log_file, 'a') as log_file_id:
        log_file_id.write(f'''Reading TWAS_weights: {twas_sparse_dir}\n''')
    twas_snps = pd.read_csv(f'''{twas_sparse_dir}/snps.txt''', header=None, sep='\t')
    twas_genes = pd.read_csv(f'''{twas_sparse_dir}/gene_names.txt''', header=None).iloc[:, 0]
    W_twas = scipy.sparse.load_npz(f'''{twas_sparse_dir}/W_TWAS.npz''')
    W_twas = W_twas.tocsc() # Much faster computation for matrix multiplication
    W_twas = W_twas.astype(precision)
    twas_snps = twas_snps.iloc[:, [0,2,1]] # Ensures that we count copies of second allele (as is shown in TWAS documentation)
    twas_snps.loc[:, 'string'] = twas_snps.apply(lambda x: ''.join(x), axis=1)

    # Read in Genetics
    with open(log_file, 'a') as log_file_id:
        log_file_id.write(f'''Reading genetics {genetics}\n''')
    G  = read_plink1_bin(genetics)
    
    # Concatonate strings with allele codes to deal with multi allele
    G_snp_str = pd.Series(np.char.add(np.char.add(G.snp.values, G.a0.values), G.a1.values))
    G_snp_str_flipped = pd.Series(np.char.add(np.char.add(G.snp.values, G.a1.values), G.a0.values))
    # Rename G variant coordinates to be 'rsidA0A1' format
    G = G.assign_coords(variant=('variant',G_snp_str))
    G = G.assign_coords(variant_flipped=('variant',G_snp_str_flipped))

    # Preallocate TWAS results dataframe
    n_indiv = len(G.sample.values)
    n_genes = len(twas_genes)
    twas_results = pd.DataFrame(np.ones((n_indiv, n_genes), dtype=precision)*np.nan,
                            index=G.sample.values,
                            columns=twas_genes)
  
    # Process in 1000 individual chunks
    windows = np.arange(0, len(G.sample), 1000).tolist() + [len(G.sample)]
    for i in range(len(windows)-1):
        with open(log_file, 'a') as log_file_id:
            log_file_id.write('Processing window of individuals ' + str(windows[i]) + ':' + str(windows[i+1]) + '\n')
        print('Processing window of individuals ' + str(windows[i]) + ':' + str(windows[i+1]))
        # Match snps
        match_snps = pd.Series(G.variant.values).isin(twas_snps.string)
        match_snps_flipped = pd.Series(G.variant_flipped.values).isin(twas_snps.string)
        all_match_snps = (match_snps | match_snps_flipped)

        sub_G = G[windows[i]:windows[i+1], all_match_snps]
        sub_G = sub_G.compute()

        # Flip SNPs that are wrong way round
        wrong_way_snps = ~pd.Series(sub_G.variant.values).isin(twas_snps.string)
        # For SNPs already been read into memory
        sub_G[:, wrong_way_snps] = (sub_G[:, wrong_way_snps].values*-1)+2 # Flip values
        # Flip a0 and a1 labels
        new_a0 = sub_G.a0.values.copy() # Deep copy
        new_a0[wrong_way_snps] = sub_G.a1.values[wrong_way_snps]
        new_a1 = sub_G.a1.values.copy() # Deep copy
        new_a1[wrong_way_snps] = sub_G.a0.values[wrong_way_snps]
        sub_G = sub_G.assign_coords(a0=('variant', new_a0), a1=('variant', new_a1))
        sub_G_snp_str = pd.Series(np.char.add(np.char.add(sub_G.snp.values, sub_G.a0.values), sub_G.a1.values))
        sub_G = sub_G.assign_coords(variant=('variant',sub_G_snp_str))

        # Remove TWAS SNPs that are not Genetics
        twas_ind = twas_snps.string.isin(sub_G.variant.values)
        twas_snps = twas_snps.iloc[twas_ind.values, :]
        W_twas = W_twas[:, twas_ind.values]
        
        # Reorder W_twas
        # reorder_i = ismember(sub_G.variant.values, twas_snps.string.values)
        reorder_i = np.where(pd.Series(sub_G.variant.values).isin(twas_snps.string.values))[0].tolist()
        twas_snps = twas_snps.iloc[reorder_i, :]
        W_twas = W_twas[:, reorder_i]
        
        ## Trim on missingness 
        print('\t Triming on missingness')   
        # Remove SNPs with more than 10% missingness
        keep_snps = ((100*sub_G.isnull().sum(axis=0).values/sub_G.shape[0])<10)
        # Remove individuals with more than 5% missingness
        keep_indiv = ((100*sub_G.isnull().sum(axis=1).values/sub_G.shape[1])<5)
        # Restrict Weight Matrix
        W_twas = W_twas[:, keep_snps]
        twas_snps = twas_snps.iloc[keep_snps, :]
        # Restrict genetics
        sub_G = sub_G[keep_indiv, keep_snps]
        # Replace reamining null in sub_G with zero 
        sub_G = sub_G.fillna(0)
        
        print('\t Imputing gene expression')
        res = sub_G * W_twas.T
        twas_results.loc[sub_G.sample.values, :] = np.array(res.values).astype(precision)
        
    # Rename genes from ensemble to gene symbols
    twas_ind, new_names = gene_map(twas_results.columns, 
                            g_build=g_build, 
                            multimatch_filt={'transcript_filter':'protein_coding','gene_version': 'newest'})
    twas_results = twas_results.iloc[:, twas_ind.values]
    twas_results.columns = new_names
    # Need to change precission - also make sure that each dimension is in similar range
    print('Saving final results')
    if '.tsv' in twas_results:
        twas_results.to_csv(output_file, sep='\t')
    else:
        twas_results.to_hdf(output_file, key='twas')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TWAS Projection')

    parser.add_argument('-genetics', help='Path to genetics with # replacing chromosome', type=str)
    parser.add_argument('-twas_sparse_dir', help='Directory of TWAS sparse form', type=str, default=None)
    parser.add_argument('-g_build', help='Genome build (37 or 38), for converting Enesemble IDs to gene names', type=int, default=37)
    parser.add_argument('-output_file', help='Output file for TWAS projection, should end in .tsv or .hdf', type=str, default='none')

    opt = parser.parse_args()

    # Print arguments to console
    for key, term in opt.__dict__.items():
        print(str(key) + ':\t' + str(term) + "\n")

    twas_sparse_dir=opt.twas_sparse_dir

    genetics = opt.genetics.replace('#', '*')
    TWAS_project(genetics, twas_sparse_dir, opt.output_file, opt.g_build)