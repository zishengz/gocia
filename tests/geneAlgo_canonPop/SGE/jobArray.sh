#!/bin/sh -f
#$ -l h_data=3G,h_rt=23:59:00
#$ -cwd
#$ -o LOG.$JOB_NAME.$JOB_ID
#$ -j y
#$ -pe dc* 8
#$ -t 1-200

nproc=8

# Get the worker name
#PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`
PROJ=$JOB_NAME
WORKERID=`printf %06d $SGE_TASK_ID`
WORKPATH=${PROJ}_${WORKERID}
CWDIR=`pwd`

echo $WORKPATH starts at $CWDIR on `date`

# SET UP ENV
VASPHOME=/u/project/ana/zisheng/program/vasp54/vasp.5.4.1/bin
source /u/home/z/zisheng/.bashrc
source /u/local/Modules/default/init/modules.sh
module load intel/16.0.2
module load intelmpi/5.0.0


mkdir $WORKPATH
cd $WORKPATH

python -u ../getKid.py
mkpot pos
cp ../KPOINTS .

for i in 1 2 3
do
    echo $PROJ opt-$i starts `date`
    cp ../INCAR-$i INCAR
    mpirun -n $nproc $VASPHOME/vasp_gam > out
    grep Elapse OUTCAR >> log
    grep F OSZICAR|tail -n $i >> log
    cp CONTCAR $CURDIR/out-$i.vasp
done

python -u ../addKid.py
rm CHG* WAVECAR vasprun* POTCAR OUTCAR DOSCAR EIGENVAL XDATCAR

cd $CWD
# rm -rf $WORKPATH

