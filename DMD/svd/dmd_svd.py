import numpy as np
import os
from matplotlib import pyplot as plt
def getInputData(inputDirPath):
    oin = []
    datalist = os.listdir(inputDirPath)
    for i in range(1,len(datalist)+1):
        oin.append(np.loadtxt('./data/predata/data_{0:0>5}.dat'.format(i)))
    return oin

def dmdSvd(oin):
    oin = np.mat(oin).T
    V1 = oin[:,:-1]
    U,Sigma,WT = np.linalg.svd(V1,full_matrices=False)
    V2 = oin[:,1:]
    Sigma = np.diag(Sigma)
    S = U.T.dot(V2).dot(WT.T).dot(np.mat(Sigma).I)

    eigens, eigenVectors = np.linalg.eig(S)

    modes = U.dot(eigenVectors)
    energy = []
    for mode in range(1,modes.shape[1]+1):
        energy.append(np.sum(np.sqrt(np.power(modes.T[:,mode].imag,2)+np.power(modes.T[:,mode].real,2))))
    energy = np.array(energy)
    np.savetxt('./data/results/svd/DMD_eigenvalue_{0:0>4}.dat'.format(oin.shape[1]),
               np.array([eigens.imag, eigens.real, energy]).T)
    modesEigens = np.column_stack((modes.T,eigens.imag))
    sortedModes = modes.T[np.lexsort(modesEigens.T)[0]].T
    for mode in range(0,sortedModes.shape[1]):
        np.savetxt('./data/results/svd/DMD_Mode_{0:0>4}{1:0>4}_Real.dat'.format(sortedModes.shape[1]+1, mode + 1),
                   sortedModes[:, mode].real)
        np.savetxt('./data/results/svd/DMD_Mode_{0:0>4}{1:0>4}_Imag.dat'.format(sortedModes.shape[1]+1, mode + 1),
                   sortedModes[:, mode].imag)

if __name__ == '__main__':
    inData = getInputData('./data/predata/')
    dmdSvd(inData)
