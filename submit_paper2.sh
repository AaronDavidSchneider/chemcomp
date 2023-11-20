#!/bin/bash -l
sbatch --cpus-per-task=3 submit_old_short.sh -j jobs/paper2_disk_C.yaml -c config/config_paper2_C.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/paper2_growth_C.yaml -c config/config_paper2_C.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/Sat.yaml -c config/config_Sat.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/Jup.yaml -c config/config_Jup.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/Sat_more_C.yaml -c config/config_Sat_more_C.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/Jup_more_C.yaml -c config/config_Jup_more_C.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/paper2_growth_more_C.yaml -c config/config_paper2_more_C.yaml -d 0 -o 1
sbatch --cpus-per-task=3 submit_old_short.sh -j jobs/paper2_disk_more_C.yaml -c config/config_paper2_more_C.yaml -d 0 -o 1
