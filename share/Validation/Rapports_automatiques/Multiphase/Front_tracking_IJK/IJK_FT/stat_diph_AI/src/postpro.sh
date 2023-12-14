#!/bin/bash

if [ -d "${project_directory}/share/PyTools" ]
then
  # Using 'external' PyTools
  source $project_directory/env_for_PyTools.sh
else
# Using 'internal' PyTools
  export PYTHONPATH=${PYTHONPATH}:${project_directory}/share/bin/PyTools
fi

python analytical.py &> analytical.log
python compare.py &> compare.log
