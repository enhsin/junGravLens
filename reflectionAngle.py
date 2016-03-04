__author__ = 'cheng109'


from astropy.io import fits
#from fastell4py import fastell4py
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pyregion
import scipy.sparse
import math
import numpy.random



def pt_model(imgX, imgY, critR):
    Xcenter = 100
    Ycenter = 100

    srcCenterX = 100
    srcCenterY = 100

    fX = imgX - Xcenter
    fY = imgY - Ycenter

    fDenom = fX*fX + fY*fY


    pDeltaX = fX * critR*critR/fDenom
    pDeltaY = fY * critR*critR/fDenom

    return srcCenterX + pDeltaX, srcCenterY + pDeltaY


def filterImage(maskFileName, imageFileName, filterType):
    imgList = []
    imgData, _ = readFitsImage(imageFileName)
    xlim, ylim = imgData.shape
    if filterType=="REG":
        reg = open(maskFileName, 'r').read()

        hdulist = pyfits.open(imageFileName)
        maskData = pyregion.parse(reg).get_mask(hdu=hdulist[0])
        #maskData = writeMaskFile(maskData, 'mask.fits')
    if filterType=="FITS":
        maskData, _= readFitsImage(maskFileName)
        assert (maskData.shape==imgData.shape), "Mask shape does not match image shape!"

    if filterType=="NONE":
        maskData = np.ones((xlim, ylim))

    for i in range(xlim):
        for j in range(ylim):
            if maskData[i][j]!=0:
                if (i+j)%2 ==1:
                    type = 'v'
                else:
                    type = 'o'
                #type = 'v'
                imgList.append((j,i, imgData[i][j], type))

    return imgList


def readFitsImage(imageName):
    hdulist = fits.open(imageName)
    data =  hdulist[0].data
    hdulist.close()
    vector = np.reshape(data, data.shape[0]*data.shape[1])
    return data, vector


def main():
    critR = 5.85
    imgX = 104
    imgY = 159
    srcX, srcY = pt_model(imgX, imgY, critR)

    imageFileName = "pt_test/jun_image.fits"
    maskFileName = "pt_test/mask_0.1.reg"
    filterType = "REG"

    data, vector = readFitsImage("pt_test/jun_image.fits")

    imgList = filterImage(maskFileName, imageFileName, filterType)





    print data
    print vector

if __name__=='__main__':
    main()
