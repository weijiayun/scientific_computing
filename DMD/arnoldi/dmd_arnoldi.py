import numpy as np,os
import matplotlib.pyplot as plt
def getInputData(inputDirPath):
    oin = []
    datalist = os.listdir(inputDirPath)
    for i in range(1,len(datalist)+1):
	print 'reading data_{0:0>5}'.format(i)
        oin.append(np.loadtxt('./data/predata/data_{0:0>5}.dat'.format(i)))
    return oin

def dmd_arnoldi(oin):
    oin = np.mat(oin).T
    V1 = oin[:, :-1]
    vn = oin[:, -1]

    Q, R = np.linalg.qr(V1)

    a = R.I.dot(Q.T).dot(vn)
    np.savetxt('./data/results/svd/DMD_arnoldi_a_coeff_{0:0>4}.dat'.format(oin.shape[1]), a)
    S = np.eye(N=V1.shape[1], k=-1)
    S[:, -1] = np.array(a.T)[0]
    eigens, eigenVectors = np.linalg.eig(S)
    modes = V1.dot(eigenVectors)

    energy = []
    for mode in range(1, modes.shape[1] + 1):
        energy.append(np.sum(np.sqrt(np.power(modes.T[:, mode].imag, 2) + np.power(modes.T[:, mode].real, 2))))

    energy = np.array(energy)
    np.savetxt('./data/results/svd/DMD_arnoldi_eigenvalue_{0:0>4}.dat'.format(oin.shape[1]),
               np.array([eigens.imag, eigens.real, energy]).T)

    modesEigens = np.column_stack((modes.T, eigens.imag))
    sortedModes = modes.T[np.lexsort(modesEigens.T)[0]].T
    print sortedModes.shape
    for mode in range(0, sortedModes.shape[1]):
        np.savetxt('./data/results/svd/DMD_Mode_{0:0>4}{1:0>4}_Real.dat'.format(sortedModes.shape[1] + 1, mode + 1),
                   sortedModes[:, mode].real)

        np.savetxt('./data/results/svd/DMD_Mode_{0:0>4}{1:0>4}_Imag.dat'.format(sortedModes.shape[1] + 1, mode + 1),
                   sortedModes[:, mode].imag)

if __name__ == '__main__':
    inData = getInputData('./data/predata/')
    dmd_arnoldi(inData)

