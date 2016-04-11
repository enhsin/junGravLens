/*
 * parafit.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: cheng109
 */


#include "gsl/gsl_multimin.h"
#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
using namespace std;
//
void gridSearch(Conf* conf, MultModelParam param, Image* dataImage, vec d, string dir, double lambdaS) {
	Model *model1 = new Model(conf, param, lambdaS);
	int i=0;

	
	//vector<vector<mixModels> > mixAllModels;
	ofstream output; 
	output.open("output.txt"); 
	for(int i=0 ; i< param.nComb ; ++i) {
		MultModelParam newParam (param);  
	 
		newParam.parameter.clear(); 
		newParam.parameter.resize(param.nLens); 

		for(int j=0; j<param.nLens; ++j) {   // max of j is 3; 
			SingleModelParam s; 
			s.name = param.mixAllModels[i][j].name; 

			if(s.name=="SIE") {
				s.critRad = param.mixAllModels[i][j].paraList[0]; 
				s.centerX = param.mixAllModels[i][j].paraList[1]; 
				s.centerY = param.mixAllModels[i][j].paraList[2]; 
				s.e       = param.mixAllModels[i][j].paraList[3]; 
				s.PA 	  = param.mixAllModels[i][j].paraList[4];
				s.core 	  = param.mixAllModels[i][j].paraList[5];  
				newParam.parameter.push_back(s); 
			}
		}
		// Do work; 
		param.printCurrentModels(i); 
		model1 = new Model(conf, newParam, lambdaS);
		model1->updatePosMapping(dataImage, conf);
		//model1->updateLensAndRegularMatrix(dataImage, conf);
		//model1->updateGradient(dataImage);
		//model1->updatePenalty(&dataImage->invC, d);

							// Write residual image:
		//model1->updateReducedResidual(dataImage);
		//vector<double> sBright = eigenV_to_cV(model1->s);
		vector<double> sBright = dataImage->dataList; 
		//double srcReg = model1->getRegularizationSrcValue(d);
		//double vegetiiReg = d.transpose()*model1->HtH*d;
		double zerothOrder 		= model1->getZerothOrderReg   (conf, sBright);
	  	double gradientOrder 	= model1->getGradientOrderReg (conf, sBright); 
		double curvatureOrder 	= model1->getCurvatureOrderReg(conf, sBright);

		++i;
		int j=0;
		Image* resImg 	= new Image(dataImage->xList, dataImage->yList, &model1->simple_res_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
		Image* modelImg = new Image(dataImage->xList, dataImage->yList, &model1->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
		Image* srcImg 	= new Image(model1->srcPosXListPixel, model1->srcPosYListPixel, &sBright, conf->srcSize[0], conf->srcSize[1], conf->bitpix);

		if(1) {
			model1->updateCritCaustic(dataImage, conf);
			Image* critImg = new Image(dataImage->xList, dataImage->yList, &model1->critical, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
								

		 	critImg->writeToFile(dir + "img_crit.fits");
			Image* causticImg = new Image(model1->srcPosXListPixel, model1->srcPosYListPixel, &model1->critical, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
			causticImg->writeToFile( dir + "img_caustic.fits");
		}	


		//resImg	->writeToFile  	(dir + "img_res_" + to_string(i) + "_" + to_string(j) +".fits");
		//modelImg->writeToFile	(dir + "img_mod_" + to_string(i) + "_" + to_string(j) +".fits");
		//model1	->writeSrcImage	(dir + "img_src_" + to_string(i) + "_" + to_string(j) +".fits", conf);
		srcImg -> writeToFile (dir + "img_src_" + to_string(i)) ; 
		cout << zerothOrder << "\t" ;
		cout << gradientOrder << "\t" ;
		cout << curvatureOrder << "\t"; 
							//cout << model1->chi2 << "\t" << model1->srcR << "\t" << model1->penalty ; 
		cout << endl;

		output << zerothOrder << "\t" ;
		output << gradientOrder << "\t" ;
		output << curvatureOrder << "\t"; 
		output << endl;
					
		delete model1;
	}
	output.close(); 

}



