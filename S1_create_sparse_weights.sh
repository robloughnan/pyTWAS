#!/bin/bash


SCORE_DIR=$1
MAKESCORE='~/Programs/fusion_twas-master/utils/make_score.R'
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


# Make score files
echo "Making score files"
cd $SCORE_DIR
find . -name "*.wgt.RDat" -type f -print0 | xargs -0 -L 1 bash -c 'Rscript ~/Programs/fusion_twas-master/utils/make_score.R $0 > "${0%.wgt.RDat}.score"'

# Make unique snps file
echo "Making unique snps file"
awk '{print $1"\t"$2"\t"$3}' *.score |  uniq > unique.snps 

# Create scipy-sparse matrix
echo "Making scipy stasfile"
cd $SCRIPTDIR
python score_to_sparse.py -score_dir $SCORE_DIR