#!/bin/bash

PARTITION=$1
SETTING=$2
NSEEDS=$3

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin

# note set cpus-per-task at 5 sample sizes * 3 inflation * 2 effect protect = 30

# * iterate over 1000 seeds in array

# use all of david's nodes for now, leave one core free per node

sbatch --array=79,81,140,142,149,153,197,230,392,399,441,470,496,566,615,622,623,714,743,776,802,840,904,963,977 \
  --partition=$PARTITION \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task=5 \
  --mem-per-cpu=6G \
  --job-name=vegrowth_sim \
  --output=/projects/dbenkes/allison/vegrowth_analysis/scratch/${SETTING}_%A_%a.out \
  --export=SETTING=$SETTING,NSEEDS=$NSEEDS \
  --wrap "/apps/R/4.4.0/bin/Rscript run_analysis_generic.R"
