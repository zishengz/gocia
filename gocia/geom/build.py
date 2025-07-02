

from gocia import geom
from ase import Atoms
from ase.io.pov import get_bondpairs
from ase.build.tools import sort
import numpy as np
import random
import copy
from scipy.spatial.transform import Rotation as R
from ase.io import read, write

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
    #print('addElemList', type(addElemList))
    if type(addElemList) is str:
        addElemList_tmp = [addElemList]
    elif type(addElemList) is list:
        addElemList_tmp = addElemList.copy()
    #print('addElemList_tmp', addElemList_tmp)
    numAds = len(addElemList_tmp)
    badStructure = True
    n_place = 0
    n_attempts = 0
    while badStructure:
        if doShuffle:
            random.shuffle(addElemList_tmp)
        ind_curr = 0
        tmpInterfc = interfc.copy()
        if rattle:
            tmpInterfc.rattle(rattleStdev, zEnhance=rattleZEnhance)
        while len(tmpInterfc) < len(interfc) + numAds and ind_curr < len(addElemList_tmp):
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
                          addElemList_tmp[ind_curr] else 0 for s in optElem])
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
            # i = np.random.choice(optList, p=weights)
            # # Coordination-adapted sampling
            # if cnCount:
            #     cn = geom.get_coordStatus(tmpInterfc.get_allAtoms())[0]
            #     if cn.max() - cn.min() != 0:
            #         cn = (cn[i] - cn.min()) / (cn.max() - cn.min())
            #         if np.random.rand() < cn * cnToler:
            #             continue
#             if addElemList_tmp[ind_curr] == tmpInterfc.get_chemical_symbols()[i]:
#                 if np.random.rand() >  1/2 + sameElemPenalty: continue
#             else:
#                 if np.random.rand() <= 1/2 + sameElemPenalty: continue
# #                if np.random.rand() < sameElemPenalty: continue
            coord = [0, 0, -1000]
            while not geom.is_withinPosLim(coord, xLim, yLim, zLim):
                i = np.random.choice(optList, p=weights)
                # Coordination-adapted sampling
                if cnCount:
                    cn = geom.get_coordStatus(tmpInterfc.get_allAtoms())[0]
                    if cn.max() - cn.min() != 0:
                        cn = (cn[i] - cn.min()) / (cn.max() - cn.min())
                        if np.random.rand() < cn * cnToler:
                            continue
                growVec = geom.rand_direction()
                blda = geom.BLDA(
                    addElemList_tmp[ind_curr],
                    tmpInterfc.get_chemical_symbols()[i],
                    sigma=bldaSigma, scale=bldaScale
                )
                growVec *= blda
                coord = tmpInterfc.get_pos()[i] + growVec
#                print(blda, growVec, coord)
#                print(i, tmpInterfc.get_pos()[i])
                n_place += 1
#            print('pass')
            tmpInterfc.merge_adsorbate(Atoms(addElemList_tmp[ind_curr], [coord]))
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
    bldaSigma=0.1, bldaScale=1, toler=0.5,
    doShuffle=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    sameFragPenalty=0,
    cnCount=False, cnToler=0.5,
    bondRejList=None,
    ljopt=False, ljstepsize=0.01, ljnsteps=400,
):
    # Translate fragments from strings into atoms objects
    # Atoms object MUST list atom connecting fragment to substrate at origin first (called the 'bridle' atom as steers whole fragment)
    frags_to_add = []
    for fragName in growList:
        if fragName == 'CO':
            frags_to_add.append(Atoms('CO',[(0, 0, 0),(0, 0, 1.15034)]))
         elif fragName == 'O':
            frags_to_add.append(Atoms('O',[(0,0,0)]))
        elif fragName == 'Rh':
            frags_to_add.append(Atoms('Rh',[(0,0,0)]))
        elif fragName == 'Ti':
            frags_to_add.append(Atoms('Ti',[(0,0,0)]))
        elif fragName == 'H':
            frags_to_add.append(Atoms('H',[(0,0,0)]))
        elif fragName =='H2O':
            frags_to_add.append(Atoms('OH2',[(0,0,0),(0.758602,0,0.504284),(-0.758602,0,0.504284)]))
        elif fragName == 'HO':
            frags_to_add.append(Atoms('OH',[(0,0,0), (0,0,0.9)]))
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
            #anchorChoices = tmpInterfc_test.get_topBufLayerList()

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
            
            if tmpInterfc_test.get_pos()[i][2] < min(zLim) - 2:
                continue

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
                ## TODO: add rotation of the frag according to growVec
                coord = tmpInterfc_test.get_pos()[i] + growVec * blda
                #print(blda, growVec, coord)
                #print(i, tmpInterfc.get_pos()[i])
                n_place += 1
                #print(coord, xLim, yLim, zLim, geom.is_withinPosLim(coord, xLim, yLim, zLim))
                if n_place == 10000:
                    print('Dead Loop located line 223 in build.py -> change zLim!!')
                    print(f'Your zLim={zLim} ... coord={coord}')

            # Might want to use other way to randomly rotate positiosn
            #fragAtms.rotate(random.random()*np.pi/2,geom.rand_direction(),center=(0,0,0))
            #fragAtms.rotate(random.random()*np.pi,growVec,center=(0,0,0))
            pos_ads = fragAtms.get_positions()
            rot = R.from_rotvec(- random.random()*np.pi * np.array([0, 0, 1]))
            pos_ads = rot.apply(pos_ads)
            #fragAtms.set_positions(geom.rotate_around_vector(fragAtms.get_positions(), np.array([0,0,1]), random.random()*np.pi))
            pos_ads = geom.rotate_around_point(pos_ads, np.array([0, 0, 0]), np.array([0, 0, 1]), growVec)
            fragAtms.set_positions(pos_ads + coord)
            # Add fragment then sort
            tmpInterfc_test.add_adsFrag(fragAtms) 
            tmpInterfc_test.sortAds_frag()
            ind_curr += 1

        # Reject if bad contacts between atoms
        if not tmpInterfc_test.has_badContact(tolerance=toler):
            badStructure = False
            if bondRejList is not None:
                mySymb = tmpInterfc_test.get_chemical_symbols()
                bps = geom.get_bondpairs(tmpInterfc_test.get_allAtoms(), min(1-toler/2, 1-toler+0.2))
                bps = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in bps]
                for br in bondRejList:
                    if br in bps or [br[1],br[0]] in bps:
                        badStructure = True
                        break
                if badStructure:
                    print('b', end='')
            if not badStructure:
                tmpInterfc = tmpInterfc_test.copy()
            else:
                print('c', end='')

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
    toler_BLmin=-0.5,
    toler_CNmax=8,
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
        testInterfc = tmpInterfc.copy()
        # get the current bond pairs
#        mySymb = testInterfc.get_chemical_symbols()
        if bondRejList is not None or bondMustList is not None:
            myBPs_old = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
            myBPs_old = [[i[0], i[1]] for i in myBPs_old]
#        myBPs_old = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
        # add adsorbate
        newAdsCoord = geom.rand_point_box(xyzLims)

        a_probe = Atoms('X', [newAdsCoord])
        surf_probe = testInterfc.get_allAtoms()
        surf_probe.extend(a_probe)
        if not (2.5 > min(surf_probe.get_distances(len(surf_probe)-1, range(len(surf_probe)-1), mic=True)) > 0.5):
            continue
        
        testInterfc.merge_adsorbate(
            Atoms(addElemList[ind_curr], [newAdsCoord]))
        myCovRad = testInterfc.get_covalRadii()
        myContact = np.array([myCovRad[i]+myCovRad[-1] for i in range(len(testInterfc)-1)])
        myDist = testInterfc.get_distances(len(testInterfc)-1, range(len(testInterfc)-1))
        bondDiff = myDist - myContact  # actual distance - ideal covalent BL
        myBond = bondDiff[bondDiff > toler_BLmin]
        myBond = myBond[myBond < toler_BLmax]
       # print(len(bondDiff[bondDiff < toler_BLmin]) == 0 , len(myBond))
        if len(bondDiff[bondDiff < toler_BLmin]) == 0 and toler_CNmin <= len(myBond) <= toler_CNmax:
            goodStruc = True
            if bondRejList is not None or bondMustList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs_now = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[i[0], i[1]] for i in myBPs_now]
                myBPs = [i for i in myBPs if i not in myBPs_old]
                myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
               # print(myBPs)
            if bondRejList is not None:
                #myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs_now]
                
                for rj in bondRejList:
                    if rj in myBPs or [rj[1], rj[0]] in myBPs:
                        #print('!has rej bond', myBPs)
                        goodStruc = False
                        break
            if bondMustList is not None and goodStruc:
                if len(myBPs) == 0:
                 #   print('!has no bond')
                    goodStruc = False
                    #break
                else:
                    # print(myBPs)
                    #print(any(x in myBPs for x in bondMustList))
                    if not any(x in myBPs or [x[1], x[0]] in myBPs for x in bondMustList):
                   #     print('!has no must bond')
                        goodStruc = False
                        #break
            if constrainTop:
                pos = testInterfc.get_pos()
                for i in testInterfc.get_bufList():
                    if np.linalg.norm(newAdsCoord - pos[i]) < 2:
                        if pos[i][2] > newAdsCoord[2]:
                            goodStruc = False
            if goodStruc:
                print(f'# Progress: {ind_curr+1}/{numAds}\t@attempt {n_attempts}')
                tmpInterfc = testInterfc
                ind_curr += 1
        # prevent dead loop
        if ind_curr == 0 and n_attempts >= 500:
            print(
                'BAD STARTING STRUCTURE! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
        if n_attempts >= 250 * numAds:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    #tmpInterfc.sort() #DO NOT SORT BY DEFAULT
    #print('%i\tplacements' % (n_attempts))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc

def boxSample_frag(
    interfc,
    growList,
    xyzLims,
    toler_BLmax=1.0,
    toler_BLmin=0.45,
    toler_CNmax=6,
    toler_CNmin=1,
    bondRejList=None,
    bondMustList=None,
    doShuffle=True,
    constrainTop=False,
    rattle=False, rattleStdev=0.05, rattleZEnhance=False,
    ljopt=False, ljstepsize=0.01, ljnsteps=400,
    bondCutoff = 0.80
):
    # PROBLEMATIC 
    
    # Translate fragments from strings into atoms objects 
    # Atoms object MUST list atom connecting fragment to substrate at origin first (called the 'bridle' atom as steers whole fragment)
    frags_to_add = []
    for fragName in growList:
        if fragName == 'CO':
            frags_to_add.append(Atoms('CO',[(0, 0, 0),(0, 0, 1.15034)]))
        elif fragName == 'H':
            frags_to_add.append(Atoms('H',[(0,0,0)]))
        elif fragName == 'O':
            frags_to_add.append(Atoms('O',[(0,0,0)]))
        elif fragName == 'Rh':
            frags_to_add.append(Atoms('Rh',[(0,0,0)]))
        elif fragName == 'Ti':
            frags_to_add.append(Atoms('Ti',[(0,0,0)]))
        elif fragName =='H2O':
            frags_to_add.append(Atoms('OH2',[(0,0,0),(0.758602,0,0.504284),(-0.758602,0,0.504284)]))
        elif fragName == 'HO':
            frags_to_add.append(Atoms('OH',[(0,0,0), (0,0,0.9)]))
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
        bondDiff = myDist / myContact  # actual distance - ideal covalent BL

        # Considered bonded if within bond limits
        myBond = bondDiff[bondDiff > 1-toler_BLmin]
        myBond = myBond[myBond < toler_BLmax]

        print(len(bondDiff[bondDiff < 1-toler_BLmin]), len(myBond), end=' | ')


        # Considered good structure if no atoms are too close and if proper number of bonds
        if len(bondDiff[bondDiff < 1-toler_BLmin]) == 0 and toler_CNmin <= len(myBond) <= toler_CNmax:
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
                if not goodStruc:
                    print('b', end='')
            
            # Or required to bond with certain atoms
            if bondMustList is not None and goodStruc:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = get_bondpairs(testInterfc.get_allAtoms(), bondCutoff)
                myBPs = [[i[0], i[1]] for i in myBPs]
                myBPs = [i for i in myBPs if i not in myBPs_old]
                if len(myBPs) == 0:
                    goodStruc = False
                else:
                    myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                    print(myBPs)
                    if not any(x in myBPs or [x[1], x[0]] in myBPs for x in bondMustList):
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
        if n_attempts >= 1000:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    tmpInterfc.sortAds_frag()
    print('%i\tplacements' % (n_attempts))
    if ljopt:
        tmpInterfc.preopt_lj(stepsize=ljstepsize, nsteps=ljnsteps)
    return tmpInterfc

def boxSample_mol(interfc,
    addMolList,
    xyzLims,
    solvMode = False,
    bindDict = None,
    bondRejList=None,
    bondMustList=None,
    doShuffle=True,
    bondCutoff = 0.8
):
    numAds = len(addMolList)
    if doShuffle:
        random.shuffle(addMolList)
    n_added = 0
    n_attempts = 0
    tmpInterfc = interfc.copy()

    trj = []

    while n_added < len(addMolList):
        n_attempts += 1
        testInterfc = tmpInterfc.copy()

        centreCoord = geom.rand_point_box(xyzLims)

        a_probe = Atoms('X', [centreCoord])
        surf_probe = testInterfc.get_allAtoms()
        surf_probe.extend(a_probe)
        if not (2.5 > min(surf_probe.get_distances(len(surf_probe)-1, range(len(surf_probe)-1), mic=True)) > 0.5):
            continue
        
        tmpAds = addMolList[n_added].copy()
        pos_ads = tmpAds.get_positions()
        # rotate along z
        rot = R.from_rotvec(- random.random()*np.pi * np.array([0, 0, 1]))
        pos_ads = rot.apply(pos_ads)
        # rotate randomly a bit
        rot_axis = geom.rand_direction()
        rot_angle = random.random()*np.pi/2
        rot = R.from_rotvec(- rot_angle * rot_axis)
        pos_ads = rot.apply(pos_ads)


        pos_ads_test = pos_ads + centreCoord

        tmpAds_test = tmpAds.copy()
        tmpAds_test.set_positions(pos_ads_test)
        testInterfc.merge_adsorbate(tmpAds_test)

        if testInterfc.has_badContact(2*(1-bondCutoff)):
            continue

        bps_all = geom.get_bondpairs(testInterfc.get_allAtoms(), scale=bondCutoff)
        list_coords = []
        for ind in range(len(testInterfc)-len(tmpAds_test), len(testInterfc)):
            list_tmp = []
            for j in geom.get_coordList(ind, bps_all):
                if j not in range(len(testInterfc)-len(tmpAds_test), len(testInterfc)):
                    list_tmp.append(j)
            list_coords.append(list_tmp)
     

        trj.append(testInterfc.get_allAtoms())


        if sum([len(c) for c in list_coords]) > 0:
            goodStruc = True
            if bondRejList is not None or bondMustList is not None:
                mySymb = testInterfc.get_chemical_symbols()
                myBPs = []
                for j in range(len(list_coords)):
                    c = list_coords[j]
                    myBPs += [[len(testInterfc)-len(tmpAds_test)+j, a] for a in c]
                myBPs = [[mySymb[bp[0]], mySymb[bp[1]]] for bp in myBPs]
                

            if bondRejList is not None:                
                for rj in bondRejList:
                    if rj in myBPs or [rj[1], rj[0]] in myBPs:
                        #print('!has rej bond', myBPs)
                        goodStruc = False
                        break
            if bondMustList is not None and goodStruc:
                if len(myBPs) == 0:
                 #   print('!has no bond')
                    goodStruc = False
                    #break
                else:
                    # print(myBPs)
                    #print(any(x in myBPs for x in bondMustList))
                    if not any(x in myBPs or [x[1], x[0]] in myBPs for x in bondMustList):
                        #print('!has no must bond')
                        goodStruc = False
                        #break
            if goodStruc:
                print('Bonds:', myBPs)
                print(f'# Progress: {n_added+1}/{numAds}\t@attempt {n_attempts}')
                tmpInterfc = testInterfc
                n_added += 1

        write('test.xyz', trj)

        # prevent dead loop
        if n_added == 0 and n_attempts >= 1000:
            print(
                'BAD STARTING STRUCTURE! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
        if n_attempts >= 500 * numAds:
            print(
                'DEAD LOOP! RESTARTING...\n(if you see this too often, try adjusting the params!)')
            return None
    #tmpInterfc.sort() #DO NOT SORT BY DEFAULT
    #print('%i\tplacements' % (n_attempts))
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
