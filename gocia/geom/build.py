

from gocia import geom
from ase import Atoms
from ase.io.pov import get_bondpairs
from ase.build.tools import sort
import numpy as np
import random
import copy
from scipy.spatial.transform import Rotation as R

def grow_adatom(
    interfc, addElemList,
    xLim=None, yLim=None, zLim=None,
    sampZEnhance=None,
    bldaSigma=0.1, bldaScale=1, toler=0.5,
    doShuffle=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    sameElemPenalty=0,
    cnCount=False, cnToler=0.5,
    ljopt=False, ljstepsize=0.01, ljnsteps=400
):
    numAds = len(addElemList)
    badStructure = True
    n_place = 0
    n_attempts = 0
    while badStructure:
        if doShuffle:
            random.shuffle(addElemList)
        ind_curr = 0
        tmpInterfc = interfc.copy()
        if rattle:
            tmpInterfc.rattle(rattleStdev, zEnhance=rattleZEnhance)
        while len(tmpInterfc) < len(interfc) + numAds and ind_curr < len(addElemList):
            optList = tmpInterfc.get_optList()
            weights = np.ones(len(optList))
            # Same-element penalty
            # value of 0 has no effect on weighting
            # set to >0 to avoid same elem sampling:
            #       1 -> same 33%, diff 66%
            #       3 -> same 25%, diff 75%
            # set to -1 for pure cluster growth (naive implementation)
            optElem = [tmpInterfc.get_chemical_symbols()[i] for i in optList]
            mult = np.ones(len(optList)) +\
                np.array([sameElemPenalty if s !=
                          addElemList[ind_curr] else 0 for s in optElem])
            if np.count_nonzero(mult) == 0:
                mult = np.ones(len(optList))
            weights *= mult
            # Z-weighted sampling enhancement
            if sampZEnhance is not None:
                # sampAEnhance = 1 makes p(zmax) = p(zmin)*2
                optZ = tmpInterfc.get_pos()[:, 2][[i for i in optList]]
                if optZ.max() - optZ.min() != 0:
                    mult = 1 + (optZ - optZ.min()) / \
                        (optZ.max() - optZ.min()) * sampZEnhance
                    weights *= mult
            weights /= weights.sum()
            i = np.random.choice(optList, p=weights)
            # Coordination-adapted sampling
            if cnCount:
                cn = geom.get_coordStatus(tmpInterfc.get_allAtoms())[0]
                if cn.max() - cn.min() != 0:
                    cn = (cn[i] - cn.min()) / (cn.max() - cn.min())
                    if np.random.rand() < cn * cnToler:
                        continue
#             if addElemList[ind_curr] == tmpInterfc.get_chemical_symbols()[i]:
#                 if np.random.rand() >  1/2 + sameElemPenalty: continue
#             else:
#                 if np.random.rand() <= 1/2 + sameElemPenalty: continue
# #                if np.random.rand() < sameElemPenalty: continue
            coord = [0, 0, -1000]
            while not geom.is_withinPosLim(coord, xLim, yLim, zLim):
                growVec = geom.rand_direction()
                blda = geom.BLDA(
                    addElemList[ind_curr],
                    tmpInterfc.get_chemical_symbols()[i],
                    sigma=bldaSigma, scale=bldaScale
                )
                growVec *= blda
                coord = tmpInterfc.get_pos()[i] + growVec
#                print(blda, growVec, coord)
#                print(i, tmpInterfc.get_pos()[i])
                n_place += 1
#            print('pass')
            tmpInterfc.merge_adsorbate(Atoms(addElemList[ind_curr], [coord]))
            ind_curr += 1
        if tmpInterfc.has_badContact(tolerance=toler):
            badStructure = True
        else:
            badStructure = False
        n_attempts += 1
    tmpInterfc.sort()
    print('%i\tplacements| %i\ttabula rasa' % (n_place, n_attempts - 1))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc

def grow_frag(
    interfc,  growList,
    xLim=None, yLim=None, zLim=None,
    sampZEnhance=None,
    bldaSigma=0.1, bldaScale=1, toler=0.2,
    doShuffle=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    sameFragPenalty=0,
    cnCount=False, cnToler=0.5,
    ljopt=False, ljstepsize=0.01, ljnsteps=400
):
    # Translate fragments from strings into atoms objects
    # Atoms object MUST list atom connecting fragment to substrate at origin first (called the 'bridle' atom as steers whole fragment)
    frags_to_add = []
    for fragName in growList:
        if fragName == 'CO':
            frags_to_add.append(Atoms('CO',[(0, 0, 0),(0, 0, 1.15034)]))
        elif fragName == 'H':
            frags_to_add.append(Atoms('H',[(0,0,0)]))
        elif fragName =='H2O':
            frags_to_add.append(Atoms('OH2',[(0,0,0),(0.758602,0,0.504284),(-0.758602,0,0.504284)]))
        else:
            print('Unknown fragment: must add option to grow {} in grow_adFrag() in build.py'.format(fragName))

    # Initialize variables to track adding fragment until successful
    tmpInterfc = interfc.copy()
    numFrags = len(growList)
    numAds = len([atm for frag in frags_to_add for atm in frag])
    badStructure = True
    n_place = 0
    n_attempts = 0
    # Until we get a good structure . . .
    while badStructure:
        if doShuffle:
            random.shuffle(frags_to_add)
        ind_curr = 0
        tmpInterfc_test = interfc.copy()
        if rattle:
            tmpInterfc_test.rattle(rattleStdev, zEnhance=rattleZEnhance)
        # Keep going until added all necessary fragments 
        while len(tmpInterfc_test) < len(interfc) + numAds and ind_curr < numFrags:
            fragAtms = copy.deepcopy(frags_to_add[ind_curr])
            # Obtain candidate atoms to anchor fragment and give them weights before selecting one
            anchorChoices = tmpInterfc_test.get_bufList()# + tmpInterfc_test.get_bridList() # Might consider using get_topLayerList() 
            weights = np.ones(len(anchorChoices))
            # Same-fragment penalty
            # value of 0 has no effect on weighting
            # set to >0 to avoid same elem sampling:
            #       1 -> same 33%, diff 66%
            #       3 -> same 25%, diff 75%
            # set to -1 for pure cluster growth (naive implementation)
            elemChoices = [tmpInterfc_test.get_chemical_symbols()[i] for i in anchorChoices]
            mult = np.ones(len(anchorChoices)) + np.array([sameFragPenalty if s != growList[ind_curr] else 0 for s in elemChoices])
            if np.count_nonzero(mult) == 0:
                mult = np.ones(len(anchorChoices))
            weights *= mult
            # Z-weighted sampling enhancement
            if sampZEnhance is not None:
                # sampAEnhance = 1 makes p(zmax) = p(zmin)*2
                optZ = tmpInterfc_test.get_pos()[:, 2][[i for i in anchorChoices]]
                if optZ.max() - optZ.min() != 0:
                    mult = 1 + (optZ - optZ.min()) / (optZ.max() - optZ.min()) * sampZEnhance
                    weights *= mult
            weights /= weights.sum()
            # Pick atom off of which to grow fragment 
            i = np.random.choice(anchorChoices, p=weights)
            # Coordination-adapted sampling
            if cnCount:
                cn = geom.get_coordStatus(tmpInterfc_test.get_allAtoms())[0]
                if cn.max() - cn.min() != 0:
                    cn = (cn[i] - cn.min()) / (cn.max() - cn.min())
                    if np.random.rand() < cn * cnToler:
                        continue
                    #WON'T THIS CONTINUE EITHER WAY?
            coord = [0, 0, -1000]
            # Grow fragment a proper bond length away from chosen atom
            while not geom.is_withinPosLim(coord, xLim, yLim, zLim):
                growVec = geom.rand_direction()
                growVec[2] = np.abs(growVec[2])
                blda = geom.BLDA(
                    fragAtms.get_chemical_symbols()[0], #NEED TO IMPROVE : assumes first atom is bridle 
                    tmpInterfc_test.get_chemical_symbols()[i],
                    sigma=bldaSigma, scale=bldaScale
                )
                growVec *= blda
                ## TODO: add rotation of the frag according to growVec
                coord = tmpInterfc_test.get_pos()[i] + growVec
                #print(blda, growVec, coord)
                #print(i, tmpInterfc.get_pos()[i])
                n_place += 1

            # Might want to use other way to randomly rotate positiosn
            fragAtms.rotate(random.random()*np.pi/2,geom.rand_direction(),center=(0,0,0))
            fragAtms.set_positions(fragAtms.get_positions() + coord)
            # Add fragment then sort
            tmpInterfc_test.add_adsFrag(fragAtms) 
            tmpInterfc_test.sortAds_frag()
            ind_curr += 1
        # Reject if bad contacts between atoms
        if not tmpInterfc_test.has_badContact(tolerance=toler):
            badStructure = False
            tmpInterfc = tmpInterfc_test.copy()
        n_attempts += 1
        # prevent dead loop
        if n_attempts >= 10000:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    tmpInterfc.sortAds_frag()
    print('%i\tplacements| %i\ttabula rasa' % (n_place, n_attempts - 1)) # what does this line supposed to communicate?
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc 

def boxSample_adatom(
    interfc, addElemList,
    xyzLims,
    toler_BLmax=0,
    toler_BLmin=-0.2,
    toler_CNmax=4,
    toler_CNmin=1,
    bondRejList=None,
    bondMustList=None,
    doShuffle=True,
    constrainTop=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    ljopt=False, ljstepsize=0.01, ljnsteps=400,
    bondCutoff = 0.85
):
    numAds = len(addElemList)
    n_attempts = 0
    if doShuffle:
        random.shuffle(addElemList)
    ind_curr = 0
    tmpInterfc = interfc.copy()
    if rattle:
        tmpInterfc.rattle(rattleStdev, zEnhance=rattleZEnhance)
    while len(tmpInterfc) < len(interfc) + numAds and ind_curr < len(addElemList):
        n_attempts += 1
        optList = tmpInterfc.get_optList()
        testInterfc = tmpInterfc.copy()
        # get the current bond pairs
#        mySymb = testInterfc.get_chemical_symbols()
        if bondMustList is not None:
            myBPs_old = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
            myBPs_old = [[i[0], i[1]] for i in myBPs_old]
#        myBPs_old = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
        # add adsorbate
        newAdsCoord = geom.rand_point_box(xyzLims)
        testInterfc.merge_adsorbate(
            Atoms(addElemList[ind_curr], [newAdsCoord]))
        myContact = testInterfc.get_contactMat()[-1][:-1]
        myDist = testInterfc.get_allDistances()[-1][:-1]
        bondDiff = myDist - myContact  # actual distance - ideal covalent BL
        myBond = bondDiff[bondDiff > toler_BLmin]
        myBond = myBond[myBond < toler_BLmax]
        if len(bondDiff[bondDiff < toler_BLmin]) == 0 and toler_CNmin <= len(myBond) <= toler_CNmax:
            goodStruc = True
            if bondRejList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                for rj in bondRejList:
                    if rj in myBPs or [rj[1], rj[0]] in myBPs:
                        goodStruc = False
                        break
            if bondMustList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[i[0], i[1]] for i in myBPs]
                myBPs = [i for i in myBPs if i not in myBPs_old]
                if len(myBPs) == 0:
                    goodStruc = False
                    #break
                else:
                    myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                    print(myBPs)
                    print(any(x in myBPs for x in bondMustList))
                    if not any(x in myBPs for x in bondMustList):
                        goodStruc = False
                        #break
            if constrainTop:
                pos = testInterfc.get_pos()
                for i in testInterfc.get_bufList():
                    if np.linalg.norm(newAdsCoord - pos[i]) < 2:
                        if pos[i][2] > newAdsCoord[2]:
                            goodStruc = False
            if goodStruc:
                tmpInterfc = testInterfc
                ind_curr += 1
        # prevent dead loop
        if n_attempts >= 10000:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    tmpInterfc.sort()
    print('%i\tplacements' % (n_attempts))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc

    
def boxSample_frag(
    interfc, growList,
    xyzLims,
    toler_BLmax=0,
    toler_BLmin=-0.2,
    toler_CNmax=4,
    toler_CNmin=1,
    bondRejList=None,
    bondMustList=None,
    doShuffle=True,
    constrainTop=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    ljopt=False, ljstepsize=0.01, ljnsteps=400,
    bondCutoff = 0.85
):
    # Translate fragments from strings into atoms objects 
    # Atoms object MUST list atom connecting fragment to substrate at origin first (called the 'bridle' atom as steers whole fragment)
    frags_to_add = []
    for fragName in growList:
        if fragName == 'CO':
            frags_to_add.append(Atoms('CO',[(0, 0, 0),(0, 0, 1.15034)]))
        elif fragName == 'H':
            frags_to_add.append(Atoms('H',[(0,0,0)]))
        elif fragName =='H2O':
            frags_to_add.append(Atoms('OH2',[(0,0,0),(0.758602,0,0.504284),(-0.758602,0,0.504284)]))
        else:
            print('Unknown fragment: must add option to grow {} in grow_adFrag() in build.py'.format(fragName))

    # Initialize variables to track adding fragment until successful
    tmpInterfc = interfc.copy()
    numFrags = len(growList)
    numAds = len([atm for frag in frags_to_add for atm in frag])
    ind_curr = 0
    n_attempts = 0
    if doShuffle:
        random.shuffle(frags_to_add)
    if rattle:
        tmpInterfc.rattle(rattleStdev, zEnhance=rattleZEnhance)
    # Keep trying until interface has right number of fragments added
    while len(tmpInterfc) < len(interfc) + numAds and ind_curr < numFrags:
        n_attempts += 1
        testInterfc = tmpInterfc.copy()
        # Get the current bond pairs
        if bondMustList is not None:
            myBPs_old = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
            myBPs_old = [[i[0], i[1]] for i in myBPs_old]
        # Generate point where try to add fragment
        newFragCoord = geom.rand_point_box(xyzLims)
        tmpFrag = frags_to_add[ind_curr].copy()
        pos_ads = tmpFrag.get_positions()
        # Rotate fragment about z-axis
        rot = R.from_rotvec(- random.random()*np.pi * np.array([0, 0, 1]))
        pos_ads = rot.apply(pos_ads)
        # Rotate fragment by random angle (less than pi/3) about random axis 
        rot = R.from_rotvec(- random.random()*np.pi/3 * geom.rand_direction())
        pos_ads = rot.apply(pos_ads)
        # Add fragment to inferface
        tmpFrag.set_positions(pos_ads + newFragCoord)
        testInterfc.add_adsFrag(tmpFrag)

        # Get contact and distance matrices for bridle atom (atom at origin listed first so use length of fragment)
        myContact = testInterfc.get_contactMat()[-len(tmpFrag)][:-len(tmpFrag)] 
        myDist = testInterfc.get_allDistances()[-len(tmpFrag)][:-len(tmpFrag)]
        bondDiff = myDist - myContact  # actual distance - ideal covalent BL
        # Considered bonded if within bond limits
        myBond = bondDiff[bondDiff > toler_BLmin]
        myBond = myBond[myBond < toler_BLmax]
        # Considered good structure if no atoms are too close and if proper number of bonds
        if len(bondDiff[bondDiff < toler_BLmin]) == 0 and toler_CNmin <= len(myBond) <= toler_CNmax:
            goodStruc = True
            # Unless required to reject bonding with certain atoms
            if bondRejList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                for rj in bondRejList:
                    if rj in myBPs or [rj[1], rj[0]] in myBPs:
                        goodStruc = False
                        break
            # Or required to bond with certain atoms
            if bondMustList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[i[0], i[1]] for i in myBPs]
                myBPs = [i for i in myBPs if i not in myBPs_old]
                if len(myBPs) == 0:
                    goodStruc = False
                else:
                    myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                    print(myBPs)
                    print(any(x in myBPs for x in bondMustList))
                    if not any(x in myBPs for x in bondMustList):
                        goodStruc = False
            # Or constrained to adding on top
            if constrainTop:
                pos = testInterfc.get_pos()
                for i in testInterfc.get_bufList():
                    if np.linalg.norm(newFragCoord - pos[i]) < 2:
                        if pos[i][2] > newFragCoord[2]:
                            goodStruc = False
            if goodStruc:
                tmpInterfc = testInterfc
                ind_curr += 1
        # Prevent dead loop
        if n_attempts >= 10000:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    tmpInterfc.sortAds_frag()
    print('%i\tplacements' % (n_attempts))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc

# Below are for symmetric cell construction
def split_sym_mirror(atoms, z_mirror, z_range=0.5):
    slab_upper = atoms[[
        a.index for a in atoms if a.position[2] >= z_mirror+z_range]]
    slab_middle = atoms[[
        a.index for a in atoms if z_mirror-z_range < a.position[2] < z_mirror+z_range]]
    bottom_list = [a.index for a in atoms if a.position[2] <= z_mirror-z_range]
    if len(bottom_list)>0:
        slab_bottom = atoms[[
            a.index for a in atoms if a.position[2] <= z_mirror-z_range]]
    else:
        slab_bottom = []
    return slab_upper, slab_middle, slab_bottom


def get_sym_mirror(atoms, z_mirror, z_cell, z_range=0.5):
    slab_asym = atoms.copy()
    del slab_asym.constraints # otherwise they wont move
    s_up, s_mid, s_btm = split_sym_mirror(slab_asym, z_mirror)
    s_btm_sym = s_up.copy()
    s_btm_sym.set_positions(
        2*np.array([0, 0, z_mirror])-np.array([-1, -1, 1])*s_btm_sym.positions)
    # Then put them together
    slab_sym = s_btm_sym.copy()
    slab_sym.extend(s_mid)
    slab_sym.extend(s_up)
    slab_sym = sort(slab_sym, tags=slab_sym.positions[:, 2])
    mycell = slab_sym.get_cell(complete=True)
    mycell[2][2] = z_cell
    slab_sym.set_cell(mycell)
    return slab_sym

def split_sym_center(atoms, z_split):
    return atoms[[a.index for a in atoms if a.position[2] >= z_split]]

def get_sym_center(atoms, pos_center, z_cell):
    slab_asym = atoms.copy()
    s_up = split_sym_center(slab_asym, pos_center[2])
    s_btm_sym = s_up.copy()
    # Otherwise the position will not be changed
    del s_btm_sym.constraints
    s_btm_sym.set_positions(
        2*np.array(pos_center)-np.array([1, 1, 1])*s_btm_sym.positions)
    # Then put them together
    slab_sym = s_btm_sym.copy()
    slab_sym.extend(s_up)
    slab_sym = sort(slab_sym, tags=slab_sym.positions[:, 2])
    mycell = slab_sym.get_cell(complete=True)
    mycell[2][2] = z_cell
    slab_sym.set_cell(mycell)
    slab_sym.wrap()
    return slab_sym




# TODO Adsorption site finder
