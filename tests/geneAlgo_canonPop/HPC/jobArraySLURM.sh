#!/bin/bash -f
#SBATCH -p RM-shared
#SBATCH -t 11:59:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem-per-cpu=4gb
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o LOG.%j
#SBATCH --array=1-200

# Get the worker name
#PROJ=`pwd |rev|awk -F "/" '{print $1}'|rev`
PROJ=$SLURM_JOB_NAME
WORKERID=`printf %06d $SLURM_ARRAY_TASK_ID`
WORKPATH=${PROJ}_${WORKERID}
CWDIR=$SLURM_SUBMIT_DIR

echo $WORKPATH starts at $CWDIR on `date`

# SET UP ENV
source /home/zisheng/.bashrc
VASPHOME=/pylon5/ch4s8jp/zisheng/program/VASP-VTST-SOL/vasp.5.4.1/bin/

mkdir $WORKPATH
cd $WORKPATH

python -u ../getKid.py
cp POSCAR inp.vasp
mkpot pos
cp ../KPOINTS .

for i in 1 2 3
do
    echo $PROJ opt-$i starts `date`
    cp ../INCAR-$i INCAR
    mpirun -n $SLURM_NTASKS $VASPHOME/vasp_gam > out
    grep Elapse OUTCAR >> log
    grep F OSZICAR|tail -n $i >> log
    cp CONTCAR out-$i.vasp
    cp CONTCAR POSCAR
done

python -u ../addKid.py
rm CHG* WAVECAR vasprun* POTCAR OUTCAR DOSCAR EIGENVAL XDATCAR

cd $CWDIR
# rm -rf $WORKPATH

