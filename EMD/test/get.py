import cv2

lena = cv2.imread('lena.png')

# lena = cv2.cvtColor(lena, cv2.COLOR_BGR2GRAY)
lena = lena[::6,::6,:]
cv2.imwrite('lenagray.png',lena)