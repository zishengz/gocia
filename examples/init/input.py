
pp_path = ''

# SGE: HOFFMAN2
#vasp_cmd = 'mpirun -n $Nprocs $VASPHOME/vasp_gam > out'

# SLURM: NERSC-CORI
vasp_cmd = 'srun -n32 -c16 --cpu_bind=cores vasp_gam > out'