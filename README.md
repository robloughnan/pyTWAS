# Introduction

This small python/linux library attempts to streamline TWAS imputation of transcripts from [FUSION](http://gusevlab.org/projects/fusion/) (see section  "Individual-level predictors"). The library first creates a sparse scipy matrix in `S1_create_sparse_weights.sh` and then uses this to matrix multiply in chunks across a sample of inividuals using `S2_TWAS_projection.py`.

# Usage

## Creating a Sparse Matrix

To create a sparse matrix run:

``S1_create_sparse_weights.sh /path/to/score_dir``

where `/path/to/score_dir` is the path to the path the to the weights directory where FUSION formatted weigths (.wgt.RDat) have been downloaded. This script converts the `.wgt.RDat` files to `.score` files, it then creates a `unique.snps` file containing the unique snps across all genes. Finally it runs `score_to_sparse.py` to convert the individual scores file to a weights scipy file (`W_TWAS.npz` a sparse matrix of the weights, `gene_names.txt` for the rows and `snps.txt` for the columns). 

`~/Programs/fusion_twas-master/utils/make_score.R` in line 12 needs to be changed to the path to that script from FUSION - can't seem to get `MAKESCORE` variable to work (probably should just repeat that function in this repository). 

N.B. This function is currently very slow (about 7 hours) but can be pre-computed for different refferences. I believe the main slowdown is the SNP matching performed on `line 78` of `score_to_sparse.py` which has an implicit for loop - although I need to line profile to verify that this is the slowdown. 

## Performing TWAS Projection at the Individual Level

To project individual genetics onto TWAS genes:

``python S2_TWAS_projection.py -genetics /path/to/genetics.bed \
-twas_sparse_dir /path/to/score_dir \
-output_file /path/to/output.h5``

`/path/to/score_dir` is the same directory used in previous step (this is where the sparse weight files are). The ouput specificed can either have extention `.h5` (for faster reading) or `.tsv`. 

This function takes about 10 minutes to process 10k individuals. Could be worth trying to speed up by correctly deploying the lazy loading and dask features of `pandas_plink` or use `plinkio` instead. 

