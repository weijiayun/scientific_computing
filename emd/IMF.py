from __future__ import division
from cubicSpline import CubicSpline as CS
import numpy as np
import scipy.interpolate as itp
class Ci(object):
    def __init__(self, ci, ri):
        self.ci = np.array(ci)
        self.ri = np.array(ri)

    def getCi(self):
        return  self.ci

    def getRi(self):
        return self.ri


class SoloIMF(object):
    def __init__(self, Yt, numOfIMF, fittingType='cs', SDLowerLimt=0, SDUpperLimit=0.3):
        self.Yt = Yt
        self.fittingType = fittingType
        self.numOfIMF = numOfIMF
        self.SDLowerLimit = SDLowerLimt
        self.SDUpperLimit = SDUpperLimit

    def getData(self):
        return self.Yt

    def getExtrama(self, Yt):
        maxima = [[0, Yt[0]]]
        minima = [[0, Yt[0]]]
        try:
            for k in range(1, len(Yt)-1):
                if (Yt[k - 1] <= Yt[k]) and (Yt[k] >= Yt[k + 1]):
                    maxima.append([k, Yt[k]])
                if (Yt[k-1] >= Yt[k]) and (Yt[k] <= Yt[k+1]):
                    minima.append([k, Yt[k]])

            maxima.append([len(Yt) - 1, Yt[-1]])
            minima.append([len(Yt) - 1, Yt[-1]])

            maxima = np.array(maxima)
            minima = np.array(minima)

            if len(maxima) >= 4:
                slope1 = (maxima[1, 1] - maxima[2, 1]) / (maxima[1, 0] - maxima[2, 0])
                tmp1 = slope1 * (maxima[0, 0] - maxima[1, 0]) + maxima[1, 1]
                if tmp1 > maxima[0, 1]:
                    maxima[0, 1] = tmp1

                slope2 = (maxima[-2, 1] - maxima[-3, 1]) / (maxima[-2, 0] - maxima[-3, 0])
                tmp2 = slope2 * (maxima[-1, 0] - maxima[-2, 0]) + maxima[-2, 1]
                if tmp2 > maxima[-1, 1]:
                    maxima[-1, 1] = tmp2
            if len(minima) >= 4:
                slope1 = (minima[1, 1] - minima[2, 1]) / (minima[1, 0] - minima[2, 0])
                tmp1 = slope1 * (minima[0, 0] - minima[1, 0]) + minima[1, 1]
                if tmp1 < minima[0, 1]:
                    minima[0, 1] = tmp1

                slope2 = (minima[-2, 1] - minima[-3, 1]) / (minima[-2, 0] - minima[-3, 0])
                tmp2 = slope2 * (minima[-1, 0] - minima[-2, 0]) + minima[-2, 1]
                if tmp2 < minima[-1, 1]:
                    minima[-1, 1] = tmp2

            return minima, maxima
        except Exception as e:
            raise e

    def getEnvelop(self, Yt):
        section = np.arange(0,len(Yt), 1)
        minima, maxima = self.getExtrama(Yt)

        if self.fittingType == 'cs':
            upperEnvelop = np.array(CS(maxima).cs(section))
            lowerEnvelop = np.array(CS(minima).cs(section))
        elif self.fittingType == 'cubic-zero':
            upperEnvelop = np.array(CS(maxima).cs(section, endCondition=2))
            lowerEnvelop = np.array(CS(minima).cs(section, endCondition=2))
        else:
            raise TypeError

        meanYt = (upperEnvelop[:,1]+lowerEnvelop[:,1])/2
        newYt = Yt - meanYt
        return newYt

           # def getIMF(self):
        #     it = 0
        #     Yt = self.getData()
        #     n = self.numOfIMF
        #     h_i = [self.getData()]
        #     c_i = []
        #     while n > 0:
        #         it += 1
        #         h_i.append(self.getEnvelop(Yt))
        #         SDl = SDu = 0
        #         for k in range(1,len(h_i[-1])):
        #             SDu += np.power(h_i[-2][k]-h_i[-1][k],2)
        #             SDl += np.power(h_i[-2][k],2)
        #         SD = SDu/SDl
        #         # log = ['iteration: {0:0>5}'.format(it),
        #         #        'fitting type: {0}'.format(self.fittingType),
        #         #        'imf length: {0}'.format(len(c_i)),
        #         #        'SD: {0}'.format(SD)]
        #         #print  '\t'.join(log)
        #         if self.SDLowerLimit <= SD <= self.SDUpperLimit:
        #             ci = h_i[-1]
        #             ri = Yt-ci
        #             c_i.append(Ci(ci, ri))
        #             n -= 1
        #             h_i = [ri]
        #             Yt = ri
        #         else:
        #             Yt = h_i[-1]
        #
        #     return c_i

    def getIMF(self):
        xend = self.getData()
        tmp = []
        for i in range(self.numOfIMF):
            xstart = xend
            for n in range(3):
                xstart = self.getEnvelop(xstart)

            xend = xend - xstart

            tmp.append(Ci(xstart, xend))

        return tmp


class BiIMF(object):
    def __init__(self, BiData, mComponentH, nComponentV, fittingType='cs', SDLowerLimt=0, SDUpperLimit=0.3):
        self.data = BiData
        self.fittingType = fittingType
        self.m = mComponentH
        self.n = nComponentV
        self.i = self.data.shape[0]
        self.j = self.data.shape[1]
        self.SDLowerLimt = SDLowerLimt
        self.SDUpperLimit = SDUpperLimit

    def getVerticalIMFs(self):
        print 'vertical imf...'
        RowsIMFComponentMat = np.zeros((self.m, self.i, self.j), dtype=float)
        for i in range(self.i):
            soloDimensionalIMF = SoloIMF(self.data[i,:], self.m, self.fittingType, self.SDLowerLimt, self.SDUpperLimit)
            rowComponents = soloDimensionalIMF.getIMF()
            for m in range(self.m):
                RowsIMFComponentMat[m,i,:] = rowComponents[m].getCi()

        return RowsIMFComponentMat

    def getHerisonIMFs(self):
        RowsIMFComponentMat = self.getVerticalIMFs()
        print 'herison imf...'
        BiIMFComponentMat = np.zeros((self.m, self.n, self.i, self.j), dtype=float)
        for j in range(self.j):
            for m in range(self.m):
                soloDimensionalIMF = SoloIMF(RowsIMFComponentMat[m,:,j], self.n, self.fittingType,  self.SDLowerLimt, self.SDUpperLimit)
                colComponents = soloDimensionalIMF.getIMF()
                for n in range(self.n):
                    BiIMFComponentMat[m, n, :, j] = colComponents[n].getCi()
        print 'done'
        return BiIMFComponentMat

    def getIMF(self):
        BiIMFComponentMat = self.getHerisonIMFs()
        l = min([self.m, self.n])
        pseudoBiIMFComponentMat = np.zeros((l, self.i, self.j), dtype=float)

        for ll in range(l):
            sum1 = np.zeros((self.i, self.j))
            for k in range(self.m):
                sum1 += BiIMFComponentMat[k, ll, :, :]
            sum2 = np.zeros((self.i, self.j))
            for k in range(ll+1, self.n):
                sum2 += BiIMFComponentMat[ll, k, :, :]

            pseudoBiIMFComponentMat[ll, :, :] = sum1 + sum2

        return pseudoBiIMFComponentMat















