#!/bin/sh -f
#$ -l h_data=4G,h_rt=23:59:00
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID.$TASK_ID
#$ -j y
#$ -pe shared 1
#$ -t 1-100

ncores=1

# Get the worker name
#PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`
PROJ=$JOB_NAME
WORKERID=`printf %06d $SGE_TASK_ID`
WORKPATH=${PROJ}_${WORKERID}
CWDIR=`pwd`

echo $WORKPATH starts at $CWDIR on `date`

# Set up xtb environment
source /u/home/z/zisheng/.bashrc
source /u/project/ana/zisheng/program/xtb-6.3.3/build/assets/config_env.bash
export XTBPATH=/u/home/z/zisheng/project-ana/program/xtb-6.3.3
export OMP_NUM_THREADS=$ncores,1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_STACKSIZE=4G
export MKL_NUM_THREADS=$ncores

# Run the job by worker script
cd $WORKPATH
time python -u $CWDIR/../worker-xtb.py $WORKPATH.vasp
cd $CWDIR
rm -rf $WORKPATH

