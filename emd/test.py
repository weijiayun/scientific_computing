from __future__ import division
import cv2, numpy as np
from emd import bi_emd, bi_eemd
import matplotlib.pyplot as plt

lena = cv2.imread('lena.jpg')
lena = cv2.cvtColor(lena, cv2.COLOR_BGR2GRAY)

biIMFs = bi_eemd(lena,
                3,
                3,
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