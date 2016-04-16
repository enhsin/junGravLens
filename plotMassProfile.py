
import matplotlib.pyplot as plt
import numpy as np

def main():

    xRange = 461
    yRange = 361

    crit1 = 4.38
    crit2 =1.63
    centerX1 = xRange/2
    centerY1 = yRange/2
    centerX2 = xRange/2 -152
    centerY2 = yRange/2 - 0
    e1 = 0.32
    e2 = 0.44
    PA1 = 176
    PA2 = 14

    core1 = 1.0
    core2 = 1.0

    x = np.arange(0, xRange, 1.0)
    y = np.arange(0, yRange, 1.0)
    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    print X.shape
    for i in range(xRange):
        for j in range(yRange):
            dR1 = np.sqrt((i-centerX1)**2 + (j-centerY1)**2 + core1**2)
            coeff1 = np.sqrt(1-e1*np.cos(2*PA1*np.pi/180))

            dR2 = np.sqrt((i-centerX2)**2 + (j-centerY2)**2 + core2**2)
            coeff2 = np.sqrt(1-e2*np.cos(2*PA2*np.pi/180))

            Z[j][i] = crit1/(dR1*coeff1) + crit2/(dR2*coeff2)
           # print dR, coeff, Z1[j][i], 4.38/(dR*coeff)


    max =  np.max(Z)
    min = np.min(Z)
    diff = max - min
    print np.min(Z)
    nlayer = 10
    levels = []

    levels = np.arange(0.02, 5.0, 0.02)
    plt.contour(Z, levels)
    plt.show()


if __name__=="__main__":
    main()