import numpy as np
class CubicSpline(object):
    '''y=a_i(x-x_i)^3+b_i(x-x_i)^2+c_i(x-x_i)+d_i'''
    def __init__(self, successiveData):
        self.x = map(lambda x:x[0], successiveData)
        self.y = map(lambda y:y[1], successiveData)
        self.h_i = [self.x[i + 1] - self.x[i] for i in range(0, len(self.x) - 1)]

    def cs(self, discreteArr, endCondition=1):
        if len(self.x) == 3:
            return self.secOrderSpline(discreteArr)
        elif len(self.x) == 2:
            return self.firstOrderSpline(discreteArr)
        else:
            if endCondition == 1 or endCondition == 'natural':
                return self.naturalCubicSpline(discreteArr)
            elif endCondition == 2 or endCondition == 'zero':
                return self.zeroCubicSpline(discreteArr)
            else:
                raise TypeError('wrong end condition')

    def secOrderSpline(self, section):
        matA = np.array([[self.x[0] ** 2, self.x[0], 1],
                         [self.x[1] ** 2, self.x[1], 1],
                         [self.x[2] ** 2, self.x[2], 1]])

        y = np.array([self.y[0], self.y[1], self.y[2]])
        coeff = np.linalg.inv(matA).dot(y)

        return map(lambda x:[x, coeff[0]*x**2+coeff[1]*x+coeff[2]], section)

    def firstOrderSpline(self, section):
        slope = (self.y[1] - self.y[0])/(self.x[1] - self.x[0])
        return map(lambda x:[x, slope*(x - self.x[0])+self.y[0]], section)







    def getSMat(self, discreteX):
        if min(discreteX) >= min(self.x) and max(discreteX) <= max(self.x):
            delta = [(self.y[i + 1] - self.y[i]) / self.h_i[i] - (self.y[i] - self.y[i - 1]) / self.h_i[i - 1]
                     for i in range(1, len(self.y) - 1)]
            delta = np.mat(delta).T

            tdMat = [[2 * (self.h_i[0] + self.h_i[1]), self.h_i[1]] + [0] * (len(self.y) - 4)]
            for i in range(1, len(self.y) - 3):
                tdMat.append(
                    [0] * (i - 1) + [self.h_i[i], 2 * (self.h_i[i] + self.h_i[i + 1]), self.h_i[i + 2]] + [0] * (len(self.y) - 5 - (i - 1)))
            tdMat.append([0] * (len(self.y) - 4) + [self.h_i[-2], 2 * (self.h_i[-2] + self.h_i[-1])])
            tdMat = np.mat(tdMat)

            SMat = 6 * np.linalg.inv(tdMat).dot(delta)
            return SMat
        else:
            raise Exception('input discrete x is not in the range of self.x')

    def naturalCubicSpline(self,discreteX):
        SMat = self.getSMat(discreteX)
        ddy1 = self.secDifference(self.x[0], self.y[0], self.x[1], self.y[1], self.x[2], self.y[2])
        ddyn = self.secDifference(self.x[-3], self.y[-3], self.x[-2], self.y[-2], self.x[-1], self.y[-1])
        SMat = np.row_stack(([[ddy1]], SMat))
        SMat = np.row_stack((SMat,[ddyn]))
        coeff = []
        for i in range(SMat.shape[0]-1):
            a = (SMat[i+1,0]-SMat[i,0])/(6*self.h_i[i])
            b = SMat[i,0]/2
            c = (self.y[i+1]-self.y[i])/self.h_i[i] - self.h_i[i]*(2*SMat[i,0]+SMat[i+1,0])/6
            d = self.y[i]
            coeff.append([a,b,c,d])
        discreteXY = []
        for x in discreteX:
            sIndex = self.getSectionIndex(x)
            y = coeff[sIndex][0]*(x-self.x[sIndex]) ** 3 \
                + coeff[sIndex][1]*(x - self.x[sIndex]) ** 2 \
                + coeff[sIndex][2]*(x - self.x[sIndex])\
                + coeff[sIndex][3]
            discreteXY.append([x,y])
        return discreteXY

    def zeroCubicSpline(self, discreteX):
        coeff = []

        for i in range(len(self.y)-1):
            coeff.append([3*(self.y[i+1]-self.y[i]),
                          -2*(self.y[i+1]-self.y[i]),
                          0,
                          self.y[i]])
        discreteXY = []
        for x in discreteX:
            sIndex = self.getSectionIndex(x)
            y = coeff[sIndex][0] * ((x - self.x[sIndex])/self.h_i[sIndex]) ** 3 \
                + coeff[sIndex][1] * ((x - self.x[sIndex])/self.h_i[sIndex]) ** 2 \
                + coeff[sIndex][2] * ((x - self.x[sIndex])/self.h_i[sIndex]) \
                + coeff[sIndex][3]
            discreteXY.append([x, y])
        return discreteXY
        
    def getSectionIndex(self,x):
        for i in range(len(self.x)-1):
            if self.x[i] <= x <= self.x[i+1]:
                return i
        raise IndexError

    def fistDifference(self,x1,y1,x2,y2):
        try:
            return (y2-y1)/(x2-x1)
        except Exception:
            raise ZeroDivisionError

    def secDifference(self, x1, y1, x2, y2, x3, y3):
        try:
            return (self.fistDifference(x1,y1,x2,y2)-self.fistDifference(x2,y2,x3,y3)) / (x3 - x1)
        except Exception:
            raise ZeroDivisionError




if __name__ == '__main__':


    print CubicSpline([[1,2],[2,3],[3,5]]).cs([1,1.5,2,2.5,3])
    def testFunc1(n=5):
        import random
        a = [1 * i + random.random() for i in range(n + 1)]
        return map(lambda x: [x, np.power(x, 3) - 8], a)

    import matplotlib.pyplot as plt
    a = np.array(testFunc1(5))
    cs = CubicSpline(a)
    a = np.array(a)
    csdata = cs.cs(a[:,0])
    csdata = np.array(csdata)
    print csdata.shape
    plt.figure(figsize=(9,6))
    plt.scatter(a[:,0],a[:,1],s=40, marker='s',alpha=0.5,c='crimson')
    plt.plot(csdata[:,0],csdata[:,1],lw=2)
    plt.show()


