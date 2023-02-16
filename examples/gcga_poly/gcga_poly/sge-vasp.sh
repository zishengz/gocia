#!/bin/bash

Nprocs=8
PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`

# must have the format of s00...01.vasp
VASPPOS=$1
WORKDIR=s$1
WORKERID=$1
JOBSUB=sub-${WORKDIR}.sh

mkdir $WORKDIR
cd $WORKDIR

cat<<'EOF' > $JOBSUB
#!/bin/sh -f
#$ -l h_data=3G,h_rt=23:59:00#,h=!n7638  #,h=!m*,h=!n6638,h=!7638
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID
#$ -j y
EOF

cat<<EOF >> $JOBSUB
##$ -pe dc* $Nprocs
#$ -pe shared $Nprocs
##$ -pe node 1

export Nprocs=$Nprocs
export VASPPOS=$VASPPOS
export WORKDIR=$WORKDIR
export WORKERID=$WORKERID
export JOBSUB=$JOBSUB
EOF

cat<<'EOF' >> $JOBSUB

# only for exclusive mode
export Nprocs=$NSLOTS
#export Nprocs=`cat /proc/cpuinfo|grep processor|wc -l`

echo RUNNING ON $HOSTNAME
cat $PE_HOSTFILE
echo number of cores: $Nprocs

# ACTIVATE CONDA ENV
source /u/home/z/zisheng/.bashrc

# SET UP VASP ENV
source /u/local/Modules/default/init/modules.sh
module purge
module load IDRE intel/17.0.7
export VASPHOME=/u/project/ana/shared/vasp/vasp.5.4.4-tpc/bin

export PYTHONPATH=$PYTHONPATH:/u/project/ana/zisheng/program/pypkg/gocia/
python -u ../ga-worker.py
echo `date` '|' `pwd |rev|awk -F "/" '{print $1}'|rev` '<<' `cat label` >> ../ga.log

EOF

qsub -N ${PROJ}-${WORKERID} $JOBSUB
cd ..
