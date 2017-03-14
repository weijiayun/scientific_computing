import numpy as np
import matplotlib.pyplot as plt

def dmd_arnoldi(Oin):
    Oin = np.mat(Oin).T
    V1 = Oin[:, :-1]
    vn = Oin[:, -1]
    Q, R = np.linalg.qr(V1)
    a = R.I.dot(Q.T).dot(vn)
    S = np.eye(N=V1.shape[0], k=-1)
    S[:, -1] = a
    eigens, eigenVectors = np.linalg.eig(S)
    modes = V1.dot(eigenVectors)

