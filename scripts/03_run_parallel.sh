#!/bin/bash
#SBATCH --job-name=ms2_cluster_$SLURM_ARRAY_TASK_ID
#SBATCH --cpus-per-task=2
#SBATCH --partition test
#SBATCH --array=1-5
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
  
julia -t $SLURM_CPUS_PER_TASK 03_run_cluster.jl $SLURM_ARRAY_TASK_ID