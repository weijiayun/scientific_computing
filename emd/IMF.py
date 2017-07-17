from __future__ import division
from cubicSpline import CubicSpline as CS
import numpy as np
class Ci(object):
    def __init__(self, ci, ri, xData,lowerEnvelop, upperEnvelop, newxData, meanEnvelop):
        self.ci = np.array(ci)
        self.ri = np.array(ri)
        self.lowerEnvelop = np.array(lowerEnvelop)
        self.upperEnvelop = np.array(upperEnvelop)
        self.xData = np.array(xData)
        self.newxData = np.array(newxData)
        self.meanEnvelop = np.array(meanEnvelop)

    def getCi(self):
        return  self.ci

    def getRi(self):
        return self.ri

    def getUpper(self):
        return self.upperEnvelop

    def getLower(self):
        return self.lowerEnvelop

    def getXData(self):
        return self.xData

    def getMeanEnvelop(self):
        return self.meanEnvelop

    def getNewxData(self):
        return self.newxData


class SoloIMF(object):
    def __init__(self, Yt, numOfIMF, fittingType='cs', SDLowerLimt=0, SDUpperLimit=0.3):
        self.t = map(lambda x: x[0], Yt)
        self.y = map(lambda x: x[1], Yt)
        self.Xt = Yt
        self.fittingType = fittingType
        self.numOfIMF = numOfIMF
        self.SDLowerLimit = SDLowerLimt
        self.SDUpperLimit = SDUpperLimit

    def getData(self):
        return self.Xt

    def extendCS(self, p1, p2, section):
        try:
            k = (p2[1]-p1[1])/(p2[0]-p1[0])
            tmp = []
            for x in section:
                y = k*(x - p1[0]) + p1[1]
                tmp.append([x, y])
            return tmp
        except ZeroDivisionError as e:
            print e

    def getEnvelop(self, Yt):

        maxima = []
        minima = []

        for k in range(1,len(Yt)-1):
            if Yt[k - 1][1] <= Yt[k][1] >= Yt[k + 1][1]:
                maxima.append(Yt[k])
            elif Yt[k-1][1] >= Yt[k][1] <= Yt[k+1][1]:
                minima.append(Yt[k])
        if len(maxima) < 5 or len(minima) <5:
            raise Exception('The number of extrima is not enough, try to improve the revolution')

        minxIndex = self.t.index(min(map(lambda x:x[0],maxima)))
        maxIndex = self.t.index(max(map(lambda x:x[0],maxima)))
        section = self.t[minxIndex:maxIndex+1]

        if self.fittingType == 'cs':
            upperEnvelop = CS(maxima).cs(section)
        elif self.fittingType == 'cubic-zero':
            upperEnvelop = CS(maxima).cs(section, endCondition=2)
        else:
            raise TypeError

        p1 = upperEnvelop[0]
        p2 = upperEnvelop[1]
        p_1 = upperEnvelop[-2]
        p_2 = upperEnvelop[-1]

        upperEnvelop = self.extendCS(p1, p2, self.t[:minxIndex]) \
                       + upperEnvelop \
                       + self.extendCS(p_1, p_2,self.t[maxIndex + 1:])

        minxIndex = self.t.index(min(map(lambda x: x[0], minima)))
        maxIndex = self.t.index(max(map(lambda x: x[0], minima)))
        section = self.t[minxIndex:maxIndex + 1]

        if self.fittingType == 'cs':
            lowerEnvelop = CS(minima).cs(section)
        elif self.fittingType == 'cubic-zero':
            lowerEnvelop = CS(minima).cs(section,endCondition=2)
        else:
            raise TypeError
        
        p1 = lowerEnvelop[0]
        p2 = lowerEnvelop[1]
        p_1 = lowerEnvelop[-2]
        p_2 = lowerEnvelop[-1]

        lowerEnvelop = self.extendCS(p1, p2, self.t[:minxIndex]) \
                       + lowerEnvelop \
                       + self.extendCS(p_1, p_2,self.t[maxIndex + 1:])

        meanCurve = [[upperEnvelop[i][0],(upperEnvelop[i][1]+lowerEnvelop[i][1])/2] for i in range(len(self.t))]

        newYt = []
        for k in range(len(self.t)):
            newYt.append([Yt[k][0], Yt[k][1] - meanCurve[k][1]])
        return {'xData':Yt,'lower':lowerEnvelop,'upper': upperEnvelop,'mean': meanCurve, 'newxData':newYt}

    def getIMF(self):
        it = 0
        Xt = self.getData()
        n = self.numOfIMF
        h_i = [self.getData()]
        c_i = []
        while n > 0:
            it += 1
            envelopInfo = self.getEnvelop(Xt)
            h_i.append(envelopInfo['newxData'])
            SDu = 0
            SDl = 0
            for k in range(1,len(h_i[-1])):
                SDu += np.power(h_i[-2][k][1]-h_i[-1][k][1],2)
                SDl += np.power(h_i[-2][k][1],2)
            SD = SDu/SDl
            log = ['iteration: {0:0>5}'.format(it),
                   'fitting type: {0}'.format(self.fittingType),
                   'imf length: {0}'.format(len(c_i)),
                   'SD: {0}'.format(SD)]
            print  '\t'.join(log)
            if self.SDLowerLimit <= SD <= self.SDUpperLimit:
                ci = h_i[-1]
                ri = [[ci[k][0], Xt[k][1]-ci[k][1]] for k in range(len(self.t))]
                c_i.append(Ci(ci, ri, envelopInfo['xData'], envelopInfo['lower'], envelopInfo['upper'], envelopInfo['newxData'], envelopInfo['mean']))
                n -= 1
                h_i = [ri]
                Xt = ri
            else:
                Xt = h_i[-1]

        return c_i


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
        RowsIMFComponentMat = np.zeros((self.m, self.i, self.j), dtype=float)
        for i in range(self.i):
            rowData = np.zeros((self.j, 2))
            rowData[:,0] = np.arange(0, self.j, 1)
            rowData[:,1] = self.data[i,:]
            soloDimensionalIMF = SoloIMF(rowData, self.m, self.fittingType, self.SDLowerLimt, self.SDUpperLimit)
            rowComponents = soloDimensionalIMF.getIMF()
            for m in range(self.m):
                RowsIMFComponentMat[m,i,:] = rowComponents[m].getCi()[:, 1]

        return RowsIMFComponentMat

    def getHerisonIMFs(self):
        RowsIMFComponentMat = self.getVerticalIMFs()
        BiIMFComponentMat = np.zeros((self.m, self.n, self.i, self.j), dtype=float)
        for j in range(self.j):
            for m in range(self.m):
                colData = np.zeros((self.i, 2))
                colData[:,0] = np.arange(0, self.i, 1)
                colData[:,1] = RowsIMFComponentMat[m,:,j]
                soloDimensionalIMF = SoloIMF(colData, self.n, self.fittingType,  self.SDLowerLimt, self.SDUpperLimit)
                colComponents = soloDimensionalIMF.getIMF()
                for n in range(self.n):
                    BiIMFComponentMat[m, n, :, j] = colComponents[n].getCi()[:, 1]

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















