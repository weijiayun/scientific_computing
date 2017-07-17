from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def sin_noise(f, phase, u):
    col = 100
    row = 150
    xhalf = 0.5
    yhalf = 0.6
    x = np.arange(-xhalf, xhalf, 0.01)
    y = np.arange(-xhalf, yhalf, 0.01)
    cx = 0
    cy = 0
    A = 0.1
    T = 1/f
    ddt = 0.01*T
    t = np.arange(0, T, ddt)
    noise = np.zeros((len(t), len(y), len(x)))

    for i in range(len(y)):
        for j in range(len(x)):
            for ti, tt in enumerate(t):
                r = np.sqrt(np.power(x[j]-cx,2)+np.power(y[i]-cy,2))
                dt = r/u
                noise[ti, i, j] = A*np.sin(2*np.pi*f*(tt-dt)+phase)

    return noise

if __name__ == '__main__':
    no = sin_noise(f=30,phase=0, u=5)
    plt.figure()
    plt.contour(no[1,:,:],colors='k')
    np.save('noiseMat', no)
    print 'show successfully'
    plt.show()




