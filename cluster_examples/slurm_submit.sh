#!/bin/bash -l
# NOTE the -l flag!
#
# Name of job
#SBATCH -J par_study_accretion-$1

# Path for stdout
#SBATCH -o accretion-%j.output

# Path for stderr
#SBATCH -e accretion-%j.err

# Get status mails about your jobs
#SBATCH --mail-user schneider@mpia.de
#SBATCH --mail-type=ALL

# Request 4 weeks run time
#SBATCH -t 671:0:0

# Allocate 4 GB of Ram per node
# --mem-per-cpu can be used to allocate RAM per CPU
#SBATCH --mem=4000

#SBATCH --partition=four-wks

#SBATCH --ntasks-per-node=1   # only start 1 task via srun because Python multiprocessing starts more tasks internally!
#SBATCH --cpus-per-task=56    # assign 56 cores to that task to this task

# avoid overbooking of the cores which might occur via NumPy/MKL threading
export OMP_NUM_THREADS=1

# Run your code in partition debug

DELETE=""
OVERWRITE=""
CONFIG_FILE=""
JOB_FILE=""

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
  -c | --config_file)
    CONFIG_FILE="-c $2"
    shift # past argument
    shift # past value
    ;;
  -j | --job_file)
    JOB_FILE="-j $2"
    shift # past argument
    shift # past value
    ;;
  -d | --delete)
    DELETE="-d $2"
    shift # past argument
    shift # past value
    ;;
  -o | --overwrite)
    OVERWRITE="-o $2"
    shift # past argument
    shift # past value
    ;;
  esac
done

conda activate chemcomp
chemcomp_pipeline $CONFIG_FILE $JOB_FILE $DELETE $OVERWRITE
