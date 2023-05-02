#!/bin/bash
# Extracts timings from the triou result files, writes result to stdout
# ./extract_result.sh TRIOU_STDOUT_FILE
triou_stdout=$1
grep "clock: Total ex" $triou_stdout | awk '{print "Total time", $4}'
awk 'BEGIN{count=0} $0~"clock Ax=B" { count=count+1; if (count>2) {print "Solver time",$(NF-1)} }' $triou_stdout