
pp_path = '/u/project/ana/shared/vasp/potentials/potpaw_PBE/'

#zLim = [12, 20]


# SGE: HOFFMAN2
vasp_cmd = 'time mpirun -n $Nprocs $VASPHOME/vasp_gam > out'



# SLURM: NERSC-CORI
#vasp_cmd = 'srun -n64 -c16 --cpu_bind=cores vasp_gam > out'

