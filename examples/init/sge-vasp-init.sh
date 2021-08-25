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
#!/bin/sh -f
#$ -l h_data=3G,h_rt=23:59:00
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID
#$ -j y
EOF

cat<<EOF >> $JOBSUB
#$ -pe dc* $Nprocs

export Nprocs=$Nprocs
export VASPPOS=$VASPPOS
export WORKDIR=$WORKDIR
export WORKERID=$WORKERID
export JOBSUB=$JOBSUB
EOF

cat<<'EOF' >> $JOBSUB
echo RUNNING ON $HOSTNAME
cat $PE_HOSTFILE

# ACTIVATE CONDA ENV
source /u/home/z/zisheng/.bashrc

# SET UP VASP ENV
source /u/local/Modules/default/init/modules.sh
module load intel/17.0.1
module load intelmpi/5.0.0
export VASPHOME=/u/project/ana/hczhai/cnsi/program/VASP-VTST-SOL/vasp.5.4.1/bin

cp ../input.py .
python -u ../init-worker.py
rm input.py

EOF

qsub -N ${PROJ}-${WORKERID} $JOBSUB
cd ..
