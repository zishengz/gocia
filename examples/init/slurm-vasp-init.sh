#!/bin/bash

Nprocs=16
PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`

# must have the format of s0001.vasp
VASPPOS=$1
WORKDIR=`echo $VASPPOS | sed 's/.vasp//g'`
WORKERID=`echo $WORKDIR | sed 's/s//g'`
JOBSUB=sub-${WORKDIR}.sh

mkdir $WORKDIR
cd $WORKDIR
mv ../$VASPPOS POSCAR

cat<<'EOF' > $JOBSUB
#!/bin/sh
#SBATCH -A m3200
#SBATCH -N 2
#SBATCH -C knl
#SBATCH -q regular
#SBATCH -t 24:00:00

EOF

cat<<EOF >> $JOBSUB

export Nprocs=$Nprocs
export VASPPOS=$VASPPOS
export WORKDIR=$WORKDIR
export WORKERID=$WORKERID
export JOBSUB=$JOBSUB
EOF

cat<<'EOF' >> $JOBSUB
echo RUNNING ON $SLURM_SUBMIT_HOST
echo $SLURM_JOB_NODELIST

# ACTIVATE CONDA ENV
source /global/homes/z/zisheng/.bashrc

# SET UP VASP ENV
module load vasp-tpc/20170629-knl
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=4

python -u ../worker.py >> ../ga.log

EOF

sbatch -J ${PROJ}-${WORKERID} < $JOBSUB
cd ..
