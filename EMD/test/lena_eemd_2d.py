from __future__ import division
import cv2, numpy as np
from emd import bi_emd, bi_eemd
from cubicSpline import CubicSpline


import matplotlib.pyplot as plt

lena = cv2.imread('lenagray.png',0)

biIMFs = bi_eemd(lena,
                5,
                5,
                numOfEnsemble=1,
                fittingType='cs',
                SDLowerLimt=0,
                SDUpperLimit=0.3)

print biIMFs.shape
lena = biIMFs[0, : , :]
cv2.namedWindow('Lena')
cv2.imshow('Lena', lena)
cv2.namedWindow('Lena2')
lena = biIMFs[1, : , :]
cv2.imshow('Lena2', lena)
cv2.namedWindow('Lena3')
lena = biIMFs[2, : , :]
cv2.imshow('Lena3', lena)
cv2.waitKey(0)