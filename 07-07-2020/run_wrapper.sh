#!/bin/bash

## These variables need reset every time
## REMEMBER to insert jobname in last line
jobpath='02_model_development_kg/PCR_case_model/cluster_scripts/2020_03_18_accumVar_highBeta_sameWk/'


## Pass jobname, output and error paths from wrapper script
sbatch --job-name=sameWk --export=jj=$jobname,oo=$outpath,ee=$epath,sp=$jobpath ./$jobpath/array_script.sbatch
