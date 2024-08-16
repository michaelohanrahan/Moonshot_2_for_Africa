#!/bin/bash
#SBATCH --job-name=ms2_cluster
#SBATCH --cpus-per-task=1
#SBATCH --partition test
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
  
julia -t $SLURM_CPUS_PER_TASK 03_run_cluster.jl $SLURM_ARRAY_TASK_ID