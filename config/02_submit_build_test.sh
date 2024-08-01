#!/bin/bash -l
#SBATCH --job-name=moonshot2_hydromt_wflow_build_all_models                                    # Job name
#SBATCH --output=/p/moonshot2-casestudy/Wflow/h7/hydromt_wflow_build_%j.log      # Standard output and error log
#SBATCH --time=48:00:00                                                                 # Job duration (hh:mm:ss)
#SBATCH --partition 16vcpu
#SBATCH --ntasks=16                                                                 # Number of tasks (analyses) to run
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

conda init
conda activate hydromt-wflow
python /p/moonshot2-casestudy/Wflow/africa/scripts/02_hydromt_wflow_build.py
