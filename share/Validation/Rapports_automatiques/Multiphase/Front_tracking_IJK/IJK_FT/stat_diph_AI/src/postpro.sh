#!/bin/bash
source $project_directory/env_for_PyTools.sh
python analytical.py &> analytical.log
python compare.py &> compare.log
