#!/bin/bash

## Pass jobname, output and error paths from wrapper script
sbatch --job-name=perfectly_specified --export=rs=02-test-perfectly_specified.R ./test_script.sbatch
sbatch --job-name=del-unif --export=rs=02-test_misspec-delay-uniform.R ./test_script.sbatch
sbatch --job-name=del-pars --export=rs=02-test-misspec-delay-pars.R ./test_script.sbatch
sbatch --job-name=gi --export=rs=02-test-misspec-gi.R ./test_script.sbatch
sbatch --job-name=inc --export=rs=02-test-misspec-inc.R ./test_script.sbatch
sbatch --job-name=tr-57 --export=rs=02-test-truncate-57.R ./test_script.sbatch
sbatch --job-name=tr-97 --export=rs=02-test-truncate-97.R ./test_script.sbatch

