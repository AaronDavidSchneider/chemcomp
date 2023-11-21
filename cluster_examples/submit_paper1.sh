#!/bin/bash -l
#sbatch --cpus-per-task=6 submit_old_short.sh -j jobs/planet_vs_noplanet.yaml -c config/config_no_mig.yaml -d 0 -o 1
#sbatch --cpus-per-task=18 submit_old_short.sh -j jobs/growthtracks_carbon.yaml -c config/config_C.yaml -d 0 -o 1
sbatch --cpus-per-task=18 submit_old_short.sh -j jobs/growthtracks.yaml -c config/config_paper1.yaml -d 0 -o 1
sbatch --cpus-per-task=4 submit_old.sh -j jobs/default.yaml -c config/config_paper1.yaml -d 0 -o 1
sbatch --cpus-per-task=96 submit_old.sh -j jobs/exo3_thorngren.yaml -c config/config_paper1.yaml -d 0 -o 1
#sbatch --cpus-per-task=96 submit_old.sh -j jobs/exo3_thorngren_C.yaml -c config/config_C.yaml -d 0 -o 1
