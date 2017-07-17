from matplotlib import pyplot as plt
import numpy as np, cv2
from emd import bi_emd
from mpl_toolkits.mplot3d import Axes3D

lena = cv2.imread('lena.jpg')
lena = cv2.cvtColor(lena, cv2.COLOR_BGR2GRAY)

biIMFs = bi_emd(lena,
                mComponentH=3,
                nComponentV=3,
                fittingType='cs',
                SDLowerLimt=0,
                SDUpperLimit=0.3)

X = np.arange(0,lena.shape[0],1)
Y = np.arange(0, lena.shape[1],1)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')
Z = lena
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

ax1 = fig.add_subplot(222, projection='3d')
lena = biIMFs[0, : , :]
Z = lena
ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

ax1 = fig.add_subplot(223, projection='3d')
lena = biIMFs[1, : , :]
Z = lena
ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

ax1 = fig.add_subplot(224, projection='3d')
lena = biIMFs[2, : , :]
Z = lena
ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
plt.savefig('3dlena.png',dpi=200)
plt.show()