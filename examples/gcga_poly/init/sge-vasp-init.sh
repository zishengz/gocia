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
mv ../${WORKDIR}.fragments fragments
mv ../$VASPPOS POSCAR

cat<<'EOF' > $JOBSUB
#!/bin/sh -f
#$ -q !pod_apollo4.q
#$ -l h_data=2G,h_rt=23:59:59
##$ -l h_data=4G,h_vmem=16G,h_rt=23:59:59
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID
#$ -j y
EOF

cat<<EOF >> $JOBSUB
##$ -pe dc* $Nprocs
#$ -pe shared $Nprocs

export Nprocs=$Nprocs
export VASPPOS=$VASPPOS
export WORKDIR=$WORKDIR
export WORKERID=$WORKERID
export JOBSUB=$JOBSUB
EOF

cat<<'EOF' >> $JOBSUB
echo RUNNING ON $HOSTNAME
cat $PE_HOSTFILE
export Nprocs=$NSLOTS

# ACTIVATE CONDA ENV
source /u/home/z/zisheng/.bashrc

# SET UP VASP ENV
source /u/local/Modules/default/init/modules.sh
module purge
module load IDRE intel/17.0.7
export VASPHOME=/u/project/ana/shared/vasp/vasp.5.4.4-tpc/bin

cp ../input.py .
python -u ../init-worker.py
rm input.py

EOF

qsub -N ${PROJ}-${WORKERID} $JOBSUB
cd ..
