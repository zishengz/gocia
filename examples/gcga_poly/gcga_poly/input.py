import math

u_she = -1.0
partial_pressure = 0.05

popSize = 30

subsPot =  -190.091

chemPotDict = {
    # E(H2,gas)/2 - ln(10)*k_B*T*pH - \delta(ZPE&TS) - eU                   
    'H': -6.7596393/2 - 2.302585093*8.61733E-05 * 298.15 * 7 - 0.24 - u_she,
    'CO': -14.4503  + 8.617333262e-5 * 298 * math.log(partial_pressure)
}

zLim = [12, 20]

pp_path = '/u/project/ana/shared/vasp/potentials/potpaw_PBE/'

vasp_cmd = 'time mpirun -n $Nprocs $VASPHOME/vasp_gam > out'

