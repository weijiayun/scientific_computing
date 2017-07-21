from matplotlib import pyplot as plt
import numpy as np, cv2
from emd import bi_emd
from mpl_toolkits.mplot3d import Axes3D


X = np.arange(0,134,1)
Y = np.arange(0, 134,1)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
for i in range(6):
    Z = np.loadtxt('c{0}.dat'.format(i))
    ax = fig.add_subplot(231+i, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')

plt.savefig('matlabLena.png',dpi=200)
plt.show()