import numpy as np

def normalize_mat(matrix):
    return (matrix - matrix.min())\
        / (matrix.max() - matrix.min())

def normalize_fgp(fgpArray):
    return (fgpArray - fgpArray.min(axis=0))\
        / (fgpArray.max(axis=0) - fgpArray.min(axis=0))

def cosSimMatrix(fgpArray):
    from sklearn.metrics.pairwise import cosine_similarity
    if fgpArray.max() == 1 and fgpArray.min() == 0:
        return cosine_similarity(fgpArray)
    else:
        return cosine_similarity(normalize_fgp(fgpArray))

def diffMatrix(d1Array):
    return d1Array[:, np.newaxis] - d1Array