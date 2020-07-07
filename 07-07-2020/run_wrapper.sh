#!/bin/bash

## Pass jobname, output and error paths from wrapper script
sbatch --export=rs=RSCRIPTPATH ./test_script.sbatch
