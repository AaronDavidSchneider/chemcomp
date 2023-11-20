#!/bin/bash -l
sbatch --cpus-per-task=3 submit_old_short.sh -j jobs/paper1_ref_disk.yaml -c config/config_paper1.yaml -d 0 -o 1
#sbatch --cpus-per-task=63 submit_old_short.sh -j jobs/paper1_ref_planets_nomig.yaml -c config/config_no_mig.yaml -d 0 -o 1
sbatch --cpus-per-task=9 submit_old_short.sh -j jobs/paper1_ref_planets_mig.yaml -c config/config_paper1.yaml -d 0 -o 1
#sbatch --cpus-per-task=63 submit_old_short.sh -j jobs/paper1_ref_planets_nomig_nogap.yaml -c config/config_no_mig.yaml -d 0 -o 1
