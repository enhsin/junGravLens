import numpy as np
import pprint

def createConf(confName): 
    header= """
## Setting primary parameters: 
omega     0.28600
lambda	  0.70000
weos	  -1.000000
hubble	  0.696000
srcZ	2.379
lenZ	0.4457


## Setting secondary parameters: 
srcRes	0.1
imgRes	0.04
potRes	0.1

## Files
imageFileName	horseshoe_test/HorseShoe_large.fits
regionFileName	horseshoe_test/mask_horseshoe_large_4.reg
varFileName	  jun_var.fits
psfFileName		horseshoe_test/HorseShoe_psf.fits
criticalName	horseshoe_test/critical.reg
causticName	horseshoe_test/caustic.reg


## Source plane size
srcX	100
srcY	100
"""

    f = open(confName, 'w')
    f.write(header)
    xRange = np.arange(-150, 150, 3)
    yRange = np.arange(-150, 150, 3)
    critRange = np.arange(0.1, 5.5, 0.1)

    PTMASS_critFrom = 0.0
    PTMASS_critTo   = 5.5
    PTMASS_critInc  = 0.1

    for x in xRange:
        for y in yRange:
            f.write("PTMASS "+ str(x)       + " " + str(x)       + " " + "0.1 "
                             + str(y)       + " " + str(y)       + " " + "0.1 "
                             + str(PTMASS_critFrom) + " " + str(PTMASS_critTo) + " " + str(PTMASS_critInc)
                             + "\n")


    f.close()
    print len(xRange)*len(yRange)*len(critRange)


def readOriginalConf(confName):
    models = {"PTMASS","SIE", "NFW", "SPEMD", "SERSIC" }
    header = ""
    mCombs = []
    f = open(confName, 'r')
    for line in f.readlines():
        items = line.split()
        if len(items)>=1:
            if items[0] not in models and line[0]!='#':
                header += line
            if items[0] in models:
                mCombs.append(line)

    f.close()
    return header, mCombs


def expand(paraFrom, paraTo, paraInc) :
    ret = []

    n = int(np.floor((paraTo - paraFrom)/paraInc))
    if n==0:
        ret.append((paraFrom, paraFrom, paraInc))
    else:
        for i in range(n):
            ret.append((paraFrom +i*paraInc, paraFrom + (i+1)*paraInc, paraInc))

    #print n, ret
    return ret

def singleMix(para):
    # for SIE model:  6 levels depth for 6 parameters of SIE model;
    singlePara = []
    for p0  in para[0]:
        for p1 in para[1]:
            for p2 in para[2]:
                for p3 in para[3]:
                    for p4 in para[4]:
                        for p5 in para[5]:
                            singlePara.append([p0[0], p0[1], p0[2],
                                          p1[0], p1[1], p1[2],
                                          p2[0], p2[1], p2[2],
                                          p3[0], p3[1], p3[2],
                                          p4[0], p4[1], p4[2],
                                          p5[0], p5[1], p5[2]])
    return singlePara


def finalMix(mix):
    # depth = 3 ;  for 3 models;
    fMix = []
    for p0 in mix[0]:
        for p1 in mix[1]:
            for p2 in mix[2]:
                fMix.append((p0, p1,p2))


    return fMix

def convertStr(fMix):
    print len(fMix)
    fMixConf = []
    for models in fMix:
        oneComb = ""
        #print models
        for singleM in models:
            tempStr = "SIE "
            for num in singleM:
                tempStr += str(num) + " "
            tempStr += "\n"
            oneComb += tempStr
        fMixConf.append(oneComb)
    print fMixConf[1:3]
    return fMixConf

def writeConf(path, header, fMixConf):
    i=0
    for conf in fMixConf:
        fileName = path + str(i)+".txt"
        f = open(fileName, 'w')
        f.write(header + "\n" + conf )
        f.close()
        i+=1











def splitCombs (mCombs):
    mix = []
    # assume all of them are SIE models;

    for model in mCombs:
        items = model.split()
        i = 1
        para = []
        while i<len(items):
            para.append(expand(float(items[i]), float(items[i+1]), float(items[i+2])))
            i += 3
        singlePara = singleMix(para)
        mix.append(singlePara)
    fMix = finalMix(mix)

    fMixConf = convertStr(fMix)
    return fMixConf


def main():


    nSplits = 100
    confName="horseshoe_test/conf.txt"

    header, mCombs = readOriginalConf(confName)
    fMixConf = splitCombs(mCombs)



    writeConf("run_file/", header, fMixConf)
    #createConf(confName)





if __name__=="__main__": 
    main()