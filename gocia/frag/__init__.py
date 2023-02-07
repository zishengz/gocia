### Helper functions for manipulating fragment lists (as list of lists)

    # Flattens list of lists
def flatten(fragList):
    return [index for frag in fragList for index in frag]

    # Refills list of lists by filling frame with content in order
def refill(fragFrame, fragContent):
    newFragList = []
    lastPos = 0
    for frag in fragFrame:
        keeper = []
        for pos in range(len(frag)):
            keeper.append(pos + lastPos)
        newFrag = []
        for pos in keeper:
            newFrag.append(fragContent[pos])
        newFragList.append(newFrag)
        lastPos += len(frag)
    return newFragList

    # Remakes list of lists by filling frame according to ruler with content
def remake(fragFrame, fragRuler, fragContent):
    newFragList = []
    if len(fragContent) != len(fragRuler) or len(fragRuler) != sum([len(frag) for frag in fragFrame]) or len(fragContent) != sum([len(frag) for frag in fragFrame]):
        raise IndexError('Oops, somehow misplaced an adsorbate index while sorting with FragList.remake()!')
    else:
        for frag in fragFrame:
                posKeeper = []
                for atmId in frag:
                    pos = fragRuler.index(atmId)
                    posKeeper.append(pos)
                posKeeper.sort()
                newFrag = []
                for pos in posKeeper:
                    newFrag.append(fragContent[pos])
                newFragList.append(newFrag)
    return newFragList   

## DO I WANT TO COMBINE THE ABOVE TWO FUNCTIONS? Making it depend on whether fragRuler=None or not 

    # Flattens list of lists then sorts resulting list
    # !!! fragList flattened and sorted should always equal adsList from Interface.get_adsList()
def flatsort(fragList):
    flatsort = flatten(fragList)
    flatsort.sort()
    return flatsort 

    # Transpose indices down from All to Ads: [[144,146],[145,147]] to [[0, 2], [1, 3]] 
    # Assumes flatsort(fragList) same as adsList
def transposeDown(fragList, surf):
    return remake(fragList,flatsort(fragList),[a - len(surf.get_subAtoms()) for a in flatsort(fragList)])

    # Transpose indices up from Ads to All: [[0, 2], [1, 3]] to [[144,146],[145,147]] 
    # Assumes flatsort(fragList) same as adsList
def transposeUp(fragList, surf):
    return remake(fragList,flatsort(fragList),[a + len(surf.get_subAtoms()) for a in flatsort(fragList)])