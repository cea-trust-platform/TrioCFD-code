#!/bin/sh
python $TRUST_ROOT/bin/KSH/preprocessor.py perf_mg.data.P  perf_mg.data
python $TRUST_ROOT/bin/KSH/preprocessor.py test_AMG.data.P test_AMG.data
python $TRUST_ROOT/bin/KSH/preprocessor.py perf_defo180_g0.1.data.P perf_defo180_g0.1.data


