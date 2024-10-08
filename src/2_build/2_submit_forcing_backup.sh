#!/bin/bash -l
#SBATCH --job-name=moonshot2_backup_forcing
#SBATCH --output=/p/moonshot2-casestudy/Wflow/africa/data/0-log/h7_model_building/hydromt_wflow_forcing_%j.log
#SBATCH --time=6:00:00
#SBATCH --partition 24vcpu # 192GB RAM
#SBATCH --ntasks=24
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

conda init
conda activate hydromt-wflow
# hydromt update wflow /p/moonshot2-casestudy/Wflow/africa/data/3-input/wflow_build/1844 -i /p/moonshot2-casestudy/Wflow/workshop-prep/2_hydromt-backup-forcing-1844.yml -d /p/moonshot2-casestudy/Wflow/africa/config/deltares_data_custom.yml -vvv
# conda deactivate
# conda activate hydromt-wflow
hydromt update wflow /p/moonshot2-casestudy/Wflow/africa/data/3-input/wflow_build/1814 -i /p/moonshot2-casestudy/Wflow/workshop-prep/2_hydromt-backup-forcing-1814.yml -d /p/moonshot2-casestudy/Wflow/africa/config/deltares_data_custom.yml -vvv
conda deactivate
conda activate hydromt-wflow
hydromt update wflow /p/moonshot2-casestudy/Wflow/africa/data/3-input/wflow_build/1816 -i /p/moonshot2-casestudy/Wflow/workshop-prep/2_hydromt-backup-forcing-1816.yml -d /p/moonshot2-casestudy/Wflow/africa/config/deltares_data_custom.yml -vvv

