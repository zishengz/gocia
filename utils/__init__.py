from gacia.data import elemSymbol, covalRadii

__all__ = [

]

def get_bond_length(elemNum1, elemNum2, multiplier = 1):
    '''
    From BLDA based on covalent radii of two atoms.
    '''
    meanBl = (covalent_radii[elemNum1] + covalent_radii[elemNum2])
    return multiplier * np.random.normal(meanBl, 0.1)

def get_rand_direction():
    randVec = np.random.rand(3) - 0.5
    return randVec / np.sqrt(np.sum(randVec**2))

def get_rand_point(pbccell):
    randScale = np.random.rand(2)
    return np.dot(randScale, pbccell[:2])

    