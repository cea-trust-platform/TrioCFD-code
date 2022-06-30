#!/bin/bash
python $TRUST_ROOT/bin/KSH/preprocessor.py transport_perio.prm.P transport_perio.prm

# Generation des cas tests paralleles
sed 's/nproc_i 1/nproc_i 2/;s/nproc_j 1/nproc_j 2/;s/nproc_k 1/nproc_k 2/;s/ConvMultiSpherePerio/ConvMultiSpherePerio_par8/' ijkft_ConvectionMultiSphere_seq.data >ijkft_ConvectionMultiSphere_par8.data
sed 's/nproc_i 1/nproc_i 2/;s/nproc_j 1/nproc_j 2/;s/nproc_k 1/nproc_k 2/;s/ConvMultiSpherePerio/ConvMultiSpherePerio_par8/' ijkft_ConvectionMultiSphere_reprise.data >ijkft_ConvectionMultiSphere_par8_reprise.data
