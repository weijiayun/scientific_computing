from __future__ import division
import cv2, numpy as np
from IMF import SoloIMF, BiIMF
import matplotlib.pyplot as plt

fs = 1000
delta = 1/fs
data = np.loadtxt('instant-span/r-010th-300')
data = data[:int(len(data)/10)]
print data.shape
signal = data[:,3]
t = np.arange(0,len(signal)*delta, delta)
sData = np.zeros((len(t),2))
sData[:,0] = t
sData[:,1] = signal

IMFInst = SoloIMF(sData,'cs', 3, 0, 0.3)
soloIMFs = IMFInst.getIMF()

fig = plt.figure(figsize=(10,8))
ax0 = fig.add_subplot(221)
ax0.plot(t, signal)

for i, elem in enumerate(soloIMFs):
    component = elem.getCi()
    ax1 = fig.add_subplot(222+i)
    ax1.plot(t, component[:,1], '-',label=r'IMF: $C_{0}$'.format(i+1),c='black',lw=2)
    ax1.set_xlabel("time in seconds")
    ax1.set_ylabel(r"$f/Hz$")
    plt.tight_layout()

fig.savefig('png/testProbe.png')
plt.show()