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

def get_matedFactor2(matedList):
    if type(matedList) is list:
        matedList = np.array(matedList)
    return np.power(0.5, matedList/4)

def get_matedFactor3(matedList):
    if type(matedList) is list:
        matedList = np.array(matedList)
    return 1 / np.power(1+matedList, 2/3)
