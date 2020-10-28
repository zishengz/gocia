import numpy as np


def get_eneFactor(eneArr):
    if type(eneArr) is list:
        eneArr = np.array(eneArr)
    eneArr = eneArr
    eneNorm = (eneArr - min(eneArr)) / (max(eneArr) - min(eneArr))
    return (1 - np.tanh(2*eneNorm - 1))/2

def get_matedFactor(matedList):
    if type(matedList) is list:
        matedList = np.array(matedList)
    return 1 / np.sqrt(1 + matedList)



