#!/bin/bash
python $TRUST_ROOT/bin/KSH/preprocessor.py test_sigma.prm.P test_sigma.prm

# Generation des cas tests paralleles
sed 's/nproc_i 1/nproc_i 2/;s/nproc_j 1/nproc_j 2/;s/nproc_k 1/nproc_k 2/;' test_sigma_seq.data > test_sigma_par8.data 
