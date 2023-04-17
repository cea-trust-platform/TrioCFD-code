#!/bin/bash
python $TRUST_ROOT/bin/KSH/preprocessor.py ijkft_stat_diph_AI.prm.P ijkft_stat_diph_AI.prm

# Generation des cas tests paralleles
# sed 's/nproc_i 1/nproc_i 2/;s/nproc_j 1/nproc_j 2/;s/nproc_k 1/nproc_k 2/' ijkft_stat_diph_AI.data > ijkft_stat_diph_AI_par8.data
