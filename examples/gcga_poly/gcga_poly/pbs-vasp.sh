#!/bin/bash

Nprocs=44
PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`

# must have the format of s00...01.vasp
VASPPOS=$1
WORKDIR=s$1
WORKERID=$1
JOBSUB=sub-${WORKDIR}.sh
CWD=`pwd`/$WORKDIR

mkdir $WORKDIR
cd $WORKDIR

cat<<'EOF' > $JOBSUB
#!/bin/bash
#PBS -A AFOSR35083MAV
#PBS -l select=1:ncpus=44:mpiprocs=44
#PBS -l walltime=24:00:00
#PBS -q standard
##PBS -N test_debug
#PBS -j oe
#PBS -M n
#PBS -m n
EOF

cat<<EOF >> $JOBSUB

export Nprocs=$Nprocs
export VASPPOS=$VASPPOS
export WORKDIR=$WORKDIR
export WORKERID=$WORKERID
export JOBSUB=$JOBSUB
export CWD=$CWD
EOF

cat<<'EOF' >> $JOBSUB

cd $CWD # enter the directory of the worker job

# RUNNING ENV FOR GOCIA
source ~/.bashrc
conda activate dp #activate the conda env for gocia if it is not the base

# SET UP VASP ENV
. $MODULESHOME/init/bash 2> /dev/null
module unload PrgEnv-cray
module load intel
module load PrgEnv-intel

cp ../input.py .
python3 -u ../ga-worker.py >> ../ga.log
rm input.py

EOF

qsub -N ${PROJ}-${WORKERID} $JOBSUB
cd ..
