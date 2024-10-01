
from gocia.data import atomNum, covalRadii
from ase.neighborlist import NeighborList
import numpy as np

def BLDA(elem1, elem2, sigma = 0.1, scale = 1):
    '''
    From BLDA based on covalent radii of two atoms.
    '''
    if type(elem1) is str and type(elem2) is str:
        elem1, elem2 = atomNum[elem1], atomNum[elem2]
    else:
        print(f'ERROR: elem1 and elem2 need to be str! now: {type(elem2)} and {type(elem2)}')
    covalBl = (covalRadii[elem1] + covalRadii[elem2])
    return np.random.normal(covalBl, sigma) * scale

def RMSD(atoms1, atoms2):
    r1, r2 = atoms1.get_positions(), atoms2.get_positions()
    return np.sqrt(((r1 - r2)**2).sum() / len(r1))

def rand_directionOld():
    randVec = np.array([1, 1, 1])
    while np.sqrt(np.sum(randVec**2)) >= 1:
        randVec = 2 * (np.random.rand(3) - 0.5)
    return randVec / np.sqrt(np.sum(randVec**2))

def rand_direction():
    '''
    Marsaglia method (https://mathworld.wolfram.com/SpherePointPicking.html)
    '''
    x1, x2 = 1, 1
    while x1**2 + x2**2 >= 1:
        x1, x2 = 2 * (np.random.rand(2) - 0.5)
    return np.array([
        2 * x1 * np.sqrt(1 - x1**2 - x2**2),
        2 * x2 * np.sqrt(1 - x1**2 - x2**2),
        1 - 2 * (x1**2 + x2**2)
    ])

def align_vectors(R, V1, V2):
    # Ensure V1 and V2 are unit vectors
    V1 = V1 / np.linalg.norm(V1)
    V2 = V2 / np.linalg.norm(V2)
    
    # Calculate the cross product and dot product of V1 and V2
    cross_prod = np.cross(V1, V2)
    dot_prod = np.dot(V1, V2)
    
    # Calculate the angle of rotation
    angle = np.arccos(dot_prod)
    
    # If the vectors are already aligned, return the original coordinates
    if np.isclose(angle, 0):
        return R
    
    # Calculate the rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -cross_prod[2], cross_prod[1]],
                  [cross_prod[2], 0, -cross_prod[0]],
                  [-cross_prod[1], cross_prod[0], 0]])
    
    I = np.eye(3)
    R_mat = I + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    
    # Rotate the set of coordinates
    R_rotated = np.dot(R, R_mat.T)
    
    return R_rotated

def rotate_around_vector(R, V2, A):
    # Ensure V2 is a unit vector
    V2 = V2 / np.linalg.norm(V2)
    
    # Compute the components of Rodrigues' rotation formula
    K = np.array([[0, -V2[2], V2[1]],
                  [V2[2], 0, -V2[0]],
                  [-V2[1], V2[0], 0]])
    
    I = np.eye(3)
    R_mat = I + np.sin(A) * K + (1 - np.cos(A)) * np.dot(K, K)
    
    # Rotate the set of coordinates
    R_rotated = np.dot(R, R_mat.T)
    
    return R_rotated

def rotate_around_point(R, P, V1, V2):
    # Ensure V1 and V2 are unit vectors
    V1 = V1 / np.linalg.norm(V1)
    V2 = V2 / np.linalg.norm(V2)
    
    # Calculate the cross product and dot product of V1 and V2
    cross_prod = np.cross(V1, V2)
    dot_prod = np.dot(V1, V2)
    
    # Calculate the angle of rotation
    angle = np.arccos(dot_prod)
    
    # If the vectors are already aligned, return the original coordinates
    if np.isclose(angle, 0):
        return R
    
    # Calculate the rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -cross_prod[2], cross_prod[1]],
                  [cross_prod[2], 0, -cross_prod[0]],
                  [-cross_prod[1], cross_prod[0], 0]])
    
    I = np.eye(3)
    R_mat = I + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    
    # Translate the points so that P is the origin
    R_translated = R - P
    
    # Rotate the translated points
    R_rotated = np.dot(R_translated, R_mat.T)
    
    # Translate the points back
    R_rotated += P
    
    return R_rotated

def is_withinPosLim(vec, xLim=None, yLim=None, zLim=None):
    if xLim is None: xLim = [-999, 999]
    if yLim is None: yLim = [-999, 999]
    if zLim is None: zLim = [-999, 999]
    condition = False
    if min(xLim) < vec[0] < max(xLim) and\
       min(yLim) < vec[1] < max(yLim) and\
       min(zLim) < vec[2] < max(zLim):
        condition = True
#    print(vec, condition)
    return condition
       
def rand_point_3D(cell, xLim=None, yLim=None, zLim=None):
    # for occasions where cell counts
    flag = False
    while not flag:
        myRand = np.dot(np.random.rand(3), cell)
        flag = is_withinPosLim(myRand, xLim, yLim, zLim)
    return myRand

def rand_point_box(xyzLims):
    # mathematically sampling in a box
    tmpRand = np.random.rand(3)
    return tmpRand * (xyzLims[:,1]-xyzLims[:,0]) + xyzLims[:,0]

def get_bondpairs(atoms, scale=1.0):
    """Get all pairs of bonding atoms

    Return all pairs of atoms which are closer than radius times the
    sum of their respective covalent radii.  The pairs are returned as
    tuples::

      (a, b, (i1, i2, i3))

    so that atoms a bonds to atom b displaced by the vector::

        _     _     _
      i c + i c + i c ,
       1 1   2 2   3 3

    where c1, c2 and c3 are the unit cell vectors and i1, i2, i3 are
    integers."""
    cutoffs = scale * covalRadii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
    nl.update(atoms)
    bondpairs = []
    for a in range(len(atoms)):
        indices, offsets = nl.get_neighbors(a)
        bondpairs.extend([(a, a2, offset)
                          for a2, offset in zip(indices, offsets)])
    return bondpairs

def get_bondVec(atoms, bondpair):
    bv = atoms.get_positions()[bondpair[1]] - \
        atoms.get_positions()[bondpair[0]]
    offset = np.dot(bondpair[2], atoms.get_cell(complete=True))
    return bv + offset


def get_coordList(index, bondpairs):
    tmp = []
    for bp in bondpairs:
        if index in [bp[0], bp[1]]:
            if index == bp[0]:
                tmp.append(bp[1])
            else:
                tmp.append(bp[0])
    return tmp

def chk_bondlength(atoms, radTol = 0.2):
    bondpairs = get_bondpairs(atoms, 0.85)
    flag = True
    for bp in bondpairs:
        bl = atoms.get_distance(bp[0],bp[1], mic=True)
        stdbl = covalRadii[atoms.get_atomic_numbers()[bp[0]]]+\
                covalRadii[atoms.get_atomic_numbers()[bp[1]]]
        if bl < stdbl * (1 - radTol):
            flag = False
            break
        else: continue
    return flag

def get_neighbors(atoms, ind, chkList=None, scale = 0.85):
    cutoffs = scale * covalRadii[atoms.numbers]
    nl = NeighborList(
        cutoffs=cutoffs,
        self_interaction=False,
        bothways=True,  # or only half the list
        )
    nl.update(atoms)
    nl = nl.get_neighbors(ind)
    myNl = [[],[]]
    for i in range(len(nl[0])):
        if chkList is not None:
            if nl[0][i] in chkList:
                myNl[0].append(nl[0][i])
                myNl[1].append(nl[1][i])
        else:
            myNl[0].append(nl[0][i])
            myNl[1].append(nl[1][i])
    return myNl

def get_coordStatus(atoms, scale=1.25):
    """Returns an array of coordination numbers and an array of existing bonds determined by
    distance and covalent radii.  By default a bond is defined as 120% of the combined radii
    or less. This can be changed by setting 'covalent_percent' to a float representing a 
    factor to multiple by (default = 1.2).
    If 'exclude' is set to an array,  these atomic numbers with be unable to form bonds.
    This only excludes them from being counted from other atoms,  the coordination
    numbers for these atoms will still be calculated,  but will be unable to form
    bonds to other excluded atomic numbers.
    """

    # Get all the distances
    distances = np.divide(atoms.get_all_distances(mic=True), scale)
    
    # Atomic Numbers
    numbers = atoms.numbers
    # Coordination Numbers for each atom
    cnList, blList, nbList = [], [], []
    cr = np.take(covalRadii, numbers)
    # Array of indices of bonded atoms.  len(bonded[x]) == cn[x]
    indices = list(range(len(atoms)))
    for i in indices:
        neibInd = []
        bondLen = []
        for ii in indices:
            # Skip if measuring the same atom
            if i == ii:
                continue
            if (cr[i] + cr[ii]) >= distances[i,ii]:
                neibInd.append(ii)
                bondLen.append(distances[i,ii])
        # Add this atoms bonds to the bonded list
        nbList.append(neibInd)
        blList.append(bondLen)
    for i in nbList:
        cnList.append(len(i))
    return np.array(cnList), nbList, blList

def cart2frac(cartPos, cell):
    cell_inv = np.linalg.inv(cell)
    return np.dot(cartPos, cell_inv)

def frac2cart(fracPos, cell):
    return np.dot(fracPos, cell)

def get_fragments_old(atoms, scale = 1.0):
    if len(atoms) == 1:
        return [[0]]    
    bps = [[bp[0], bp[1]] for bp in get_bondpairs(atoms, scale=scale)]
    frags = bps.copy()
    while sum([len(f) for f in frags]) > len(atoms):
        flag = False
        for i in range(len(frags)):
            for j in range(i+1,len(frags)):
                if len(set(frags[i]) & set(frags[j]))>0:
                    frags[i] = list(set(frags[i]) | set(frags[j]))
                    frags.remove(frags[j])
                    flag = True
                    break
            if flag == True:
                break

    if sum([len(f) for f in frags]) < len(atoms):
        #print('adding back single atoms...')
        frag_single = []
        for i in range(len(atoms)):
            if not any(i in l for l in frags):
                frag_single.append([i])
        frags += frag_single
    return frags

def get_unique_indices(frags):
    tmp = set([])
    for i in range(len(frags)):
        tmp = tmp | set(frags[i])
    return tmp

def get_fragments(atoms, scale = 1.0):
    if len(atoms) == 1:
        return [[0]]    
    bps = [[bp[0], bp[1]] for bp in get_bondpairs(atoms, scale=scale)]
    frags = bps.copy()

    # add back the single atoms that are not bonded to anything
    for i in range(len(atoms)):
        for j in frags:
            if i in j:
                continue
        frags.append([i])

    # iterate to merge all frags that share any atoms
    while sum([len(f) for f in frags]) > len(atoms):
        flag = False
        for i in range(len(frags)):
            for j in range(i+1,len(frags)):
                if len(set(frags[i]) & set(frags[j]))>0:
                    frags[i] = list(set(frags[i]) | set(frags[j]))
                    frags.remove(frags[j])
                    flag = True
                    break
            if flag == True:
                break

    return frags

def del_freeMol(atoms, list_keep=[0], scale = 1.0):
    tmpAtoms = atoms.copy()
    myfrags = get_fragments(tmpAtoms, scale=scale)
    deadList = []
    if len(myfrags) > 1:
        for l in myfrags:
            if not set(l) & set(list_keep):
                print(f'Remove fragment: {tmpAtoms[l].get_chemical_symbols()}, {l}')
                deadList+=l
        del tmpAtoms[deadList]
    return tmpAtoms, deadList


def detect_bond_between_adsFrag(atoms, fragList):
    # returns: [3, 8, [{1, 2}, {1, 3}]]
    fragAtoms = [atoms[f] for f in fragList]
    bondList = []
    for i in range(len(fragAtoms)):
        for j in range(len(fragAtoms)):
            if i<j:
                bonds_i = [set(l[:2]) for l in get_bondpairs(fragAtoms[i])]
                bonds_j = [set(l[:2]) for l in get_bondpairs(fragAtoms[j])]
                bonds_j = [set([i + len(fragAtoms[i]) for i in l]) for l in bonds_j]
                frags_combined = fragAtoms[i] + fragAtoms[j]
                bonds_merge = bonds_i + bonds_j
                bonds_c = [set(l[:2]) for l in get_bondpairs(frags_combined)]
                bonds_inter = [l for l in bonds_c if l not in bonds_merge]
                if len(bonds_inter) > 0:
                    bondList.append([i, j, bonds_inter])
    return bondList

def get_covalRadii(atoms):
        return [covalRadii[i] for i in atoms.numbers]

def get_contactMat(atoms, scale=1.0):
    '''
    returns the standard covalent bondlength between all atoms
    '''
    cvRad = get_covalRadii(atoms)
    contact = np.tile(cvRad, [len(cvRad), 1])
    contact += contact.T
    contact *= scale
    np.fill_diagonal(contact, 0)
    return contact

def has_badContact(atoms, tolerance=0):
    diff = atoms.get_all_distances(mic=True) - get_contactMat(atoms, scale=1-tolerance)
    return diff.min() < 0

def is_bonded(atoms, list1, list2, scale=1.1):
    flag = False
    cvRad = get_covalRadii(atoms)
    for i in list1:
        for j in list2:
            if atoms.get_distance(i, j, mic=True) <= scale * (cvRad[i] + cvRad[j]):
                flag = True
                break
        if flag == True:
            break
    return flag





