/*
 * parafit2.cpp
 *
 *  Created on: Apr 22, 2016
 *      Author: En-Hsin Peng
 */


#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
#include <limits>
#include <iomanip>
#include "mc.h"

using namespace std;



void mcFit(Conf* conf, MultModelParam param_old, Image* dataImage, vec d, string dir, string outputFileName) {
    MC mc(123);
    size_t nLoops(1000), iter(0);
    double cfac(1.), weight(5e-2), zerothOrder, L0;
    double L(std::numeric_limits<double>::min());
    double LMax(L);

	Model *model = new Model(conf, param_old, 0.1);
	vector<vector<double> > critical;
	vector<double> sBright;

	ofstream output;
	output.open(outputFileName);

    for (size_t loop=0; loop<nLoops; ++loop) {
        iter+=1;
        L0 = L;
        cfac = mc.stepPar(model->param, cfac, iter);
        model->copyParam(3);
        sBright = dataImage->dataList;
        model->updatePosMapping(dataImage, conf);
        zerothOrder = model->getZerothOrderReg(conf, sBright);
        L = zerothOrder;

        cout<< loop<< " " << iter << ": " << L << " " << L0 << " " << LMax << " cfac "<<cfac<<endl;
        if (std::isnan(L) || (L<=L0 && mc.random() > exp((L-L0)*weight))) {
            L = L0;
        } else {
            model->copyParam(3, 4);
            if (L> LMax) {
                model->copyParam(3, 5);
                LMax = L;
            }
           output << model->param.printCurrentModels(4).at(0)
                  << std::scientific << std::setprecision(3) << L << endl;
        }
    }
    model->copyParam(5);
    sBright = dataImage->dataList;
    model->updatePosMapping(dataImage, conf);

	if(conf->outputSrcImg) {
		Image* srcImg 	= new Image(model->srcPosXListPixel, model->srcPosYListPixel, &sBright, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
		string outputSrcName = dir + "img_src_mc.fits";
		if (conf->srcBackground) {
			srcImg -> writeToFile(outputSrcName, conf->back_mean, conf->back_std ) ;
		}
		else
			srcImg -> writeToFile(outputSrcName) ;
		delete srcImg;
	}
	if(conf->outputModImg) {
		//model->updateLensAndRegularMatrix(dataImage, conf);
		//model->updateGradient(dataImage);
		//model->updatePenalty(&dataImage->invC, d);
		Image fullResidualImage = model->getFullResidual(dataImage);
		Image* modelImg = new Image(dataImage->xList, dataImage->yList, &model->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
		fullResidualImage.writeToFile(dir + "img_res_mc.fits");
		modelImg->writeToFile	(dir + "img_mod_mc.fits");
		delete modelImg;
	}
	if(conf->outputCritImg) {
		vector<Image* > curve =  getCritCaustic(conf, &model->param);
		Image* critImg = curve[0];
		Image* causImg = curve[1];
		critImg->writeToFile(dir + "img_crit_mc.fits");
		causImg->writeToFile(dir + "img_caus_mc.fits");
		delete critImg;
		delete causImg;
	}
	if(conf->outputLensImg) {
		Image* lensImg = createLensImage(conf, &model->param);
		lensImg->writeToFile(dir + "img_lens_mc.fits");
		delete lensImg;
	}
	output.close();

	// Print out the best model :
	cout << "************************\nThe best models : "<< endl;
	cout << model->param.printCurrentModels(5).at(1);
	cout << "************************\n" << endl;
	cout << model->param.printCurrentModels(6).at(1);
	cout << model->param.printCurrentModels(7).at(1);
	cout << model->param.printCurrentModels(8).at(1);
	delete model;

}



