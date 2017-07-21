from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from emd import bi_emd

if __name__ == '__main__':

    no = np.load('noiseMat.npy')
    xhalf = 0.5
    yhalf = 0.6
    x = np.arange(-xhalf, xhalf, 0.01)
    y = np.arange(-xhalf, yhalf, 0.01)
    X, Y = np.meshgrid(x, y)
    lena = no[1, :, :]+ no[1, :, :] +np.random.standard_normal(no[1,:,:].shape)*0.1
    biIMFs = bi_emd(lena,
                     2,
                     2,
                     fittingType='cs',
                     SDLowerLimt=0,
                     SDUpperLimit=0.3)



    fig = plt.figure()
    ax = fig.add_subplot(221, projection='3d')
    Z = no[1,: , :]
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    ax = fig.add_subplot(222, projection='3d')
    Z = lena
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    ax = fig.add_subplot(223, projection='3d')

    print biIMFs.shape

    Z = biIMFs[0, :, :]
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    ax = fig.add_subplot(224, projection='3d')

    Z = biIMFs[1, :, :]
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    plt.tight_layout()
    plt.show()




