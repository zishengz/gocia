
popSize = 30

subsPot = -410.2023+12*(-5.8667)

chemPotDict = {
        'O': -6.268694104,  # From B2O3 - uB
        'H': -4.21,  # From C3H8 - C3H7
}

zLim = [8, 13]

pp_path = ''

vasp_cmd = 'mpirun -n $Nprocs $VASPHOME/vasp_gam > out'
