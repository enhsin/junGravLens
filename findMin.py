

import matplotlib.pyplot as plt



def readTextFile(fileName):
    f = open(fileName, 'r')
    flag = 1


    for line in f.readlines():
        temp = line.split()
        if flag==1:
            numCols = len(temp)
            cols = [[] for i in range(numCols)]
            flag =0

        for i in range(numCols):

            cols[i].append(float(temp[i]))
    return cols




def main():
    fileName = "output.txt"
    cols = readTextFile(fileName)



    X = cols[0]
    Y = cols[1]
    Z = cols[2]

    uX = list(set(X))
    uY = list(set(Y))
    uZ = list(set(Z))

    #plt.contour(X, Y, Z)
    #plt.show()

    plt.plot(X, Y, '-o')
    plt.show()
    n = 3 # vegetti regularization:
    n = 4 # zeroth:
    n = 5 # gradient;
    n = 6 # curvature;

    n = 6

    index = cols[n].index(max(cols[n]))
    for i in range(len(cols)):
        print cols[i][index]


    #for y in uY:
    #    for z in uZ:
    #        for i in range(len(X)):
    #            if y==Y[i] and z==Z[i]:
    #                #plt.plot( X[i], cols[3][i], '-o')







if __name__=='__main__':
    main()
