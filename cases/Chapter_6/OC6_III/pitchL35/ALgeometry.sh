#!/bin/bash
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

rm -rf turbineGeometry
#rm log.foamListTimes
runApplication foamListTimes

region='turbine'
components=('blade1' 'blade2' 'blade3' 'tower')

python3.8 scripts/ALgeometry.py $region "${components[@]}"