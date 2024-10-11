#!/bin/bash -l
#SBATCH --job-name=moonshot2_hydromt_wflow_build_all_models
#SBATCH --output=/p/moonshot2-casestudy/Wflow/africa/data/0-log/h7_model_building/hydromt_wflow_build_%j.log
#SBATCH --time=6:00:00
#SBATCH --partition 24vcpu # 192GB RAM
#SBATCH --ntasks=24
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

conda init
conda activate hydromt-wflow
python /p/moonshot2-casestudy/Wflow/africa/src/2-build/2_hydromt_wflow_build.py