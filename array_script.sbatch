#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --output=messages/%x_%a.out
#SBATCH --error=messages/%x_%a.err
#SBATCH --nodes=1
#SBATCH --partition=cobey
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=500
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=kgostic
#SBATCH --job-name=perfectly_specified

module load R/4.0.0

## Print some useful variables created by SLURM
echo "SLURM_JOB_NAME is" $SLURM_JOB_NAME
echo "SLURM_JOB_ID is" $SLURM_JOB_ID
echo "SLURM_ARRAY_TASK_ID is" $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID is" $SLURM_ARRAY_JOB_ID

sleep 3


Rscript test-perfectly_specified.R  

sacct -X --format="JobName,JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem,Start,End" -j $SLURM_ARRAY_JOB_ID > messages/$SLURM_JOB_NAME_runtimes.txt
