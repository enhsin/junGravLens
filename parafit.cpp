/*
 * parafit.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: cheng109
 */


//#include "gsl/gsl_multimin.h"
#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
using namespace std;
//
void gridSearch(Conf* conf, MultModelParam param, Image* dataImage, vec d, string dir, string outputFileName) {
	Model *model1 = new Model(conf, param, 0.1);
	

	vector<int> maxIndex(3,0); 
	vector<double> maxObjFunc(3, -1.0); //maximum of a double value; 
	//vector<vector<mixModels> > mixAllModels;
	ofstream output; 
	output.open(outputFileName); 
	MultModelParam newParam (param);  
	for(int i=0 ; i< param.nComb ; ++i) {
		

		newParam.parameter.clear(); 
		//newParam.parameter.resize(param.nLens); 

		for(int j=0; j<param.nLens; ++j) {   // max of j is 3; 
			SingleModelParam s; 
			s.name = param.mixAllModels[i][j].name; 
			if(s.name=="PTMASS") {

				
				s.critRad = param.mixAllModels[i][j].paraList[0]; 
				s.centerX = param.mixAllModels[i][j].paraList[1]; 
				s.centerY = param.mixAllModels[i][j].paraList[2]; 
				newParam.parameter.push_back(s); 
			}


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
		// Do serious work; 
		//param.printCurrentModels(i); 
		vector<double> sBright = dataImage->dataList; 

		Image* srcImg; 
		model1 = new Model(conf, newParam, 0.1);
		model1->updatePosMapping(dataImage, conf);

		if(1) {
			// Output src image: 
			srcImg 	= new Image(model1->srcPosXListPixel, model1->srcPosYListPixel, &sBright, conf->srcSize[0], conf->srcSize[1], conf->bitpix);		
			srcImg -> writeToFile(dir + "img_src_" + to_string(i) +".fits" ) ; 
			delete srcImg; 
		}
		

		if (0) {
			//  Output model and residual image: 
			model1->updateLensAndRegularMatrix(dataImage, conf);
			model1->updateGradient(dataImage);
			model1->updatePenalty(&dataImage->invC, d);
			Image fullResidualImage = model1->getFullResidual(dataImage);
			Image* modelImg = new Image(dataImage->xList, dataImage->yList, &model1->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			fullResidualImage.writeToFile(dir + "img_res_" + to_string(i) + ".fits");
			modelImg->writeToFile	(dir + "img_mod_" + to_string(i) +".fits");
			delete modelImg; 

			// Output critical and caustic image: 
			model1->updateCritCaustic(dataImage, conf);
			Image* critImg = new Image(dataImage->xList, dataImage->yList, &model1->critical, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
		 	critImg->writeToFile(dir + "img_crit.fits");
			Image* causticImg = new Image(model1->srcPosXListPixel, model1->srcPosYListPixel, &model1->critical, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
			causticImg->writeToFile( dir + "img_caustic.fits");
			delete critImg; 
			
		}

	

		//vector<double> sBright = eigenV_to_cV(model1->s);
		
		//double srcReg = model1->getRegularizationSrcValue(d);
		//double vegetiiReg = d.transpose()*model1->HtH*d;
		double zerothOrder 		= model1->getZerothOrderReg   (conf, sBright);
	  	double gradientOrder 	= model1->getGradientOrderReg (conf, sBright); 
		double curvatureOrder 	= model1->getCurvatureOrderReg(conf, sBright);
		double scatter 			= model1->getScatterReg(); 
		if(zerothOrder > maxObjFunc[0])  {
			maxObjFunc[0] = zerothOrder; 
			maxIndex[0] = i; 
		}
		if(gradientOrder > maxObjFunc[1])  {
			maxObjFunc[1] = gradientOrder; 
			maxIndex[1] = i; 
		}
		if(curvatureOrder > maxObjFunc[2])  {
			maxObjFunc[2] = curvatureOrder; 
			maxIndex[2] = i; 
		}
		string pStatus = "[" + to_string(i+1) + "/" + to_string(param.nComb) + "]\t" ; 
		string resultStatus =  to_string(newParam.parameter[2].centerY)  + "\t"
				+ to_string(scatter*1000) + "\t"
				+ to_string(zerothOrder) + "\t" 
				+ to_string(gradientOrder) + "\t" 
				+ to_string(curvatureOrder) + "\n"; 


		cout 	<< pStatus << resultStatus ; 
		output  << pStatus << param.printCurrentModels(i).at(0) << resultStatus ; 		
		delete model1;
	}
	output.close(); 

	// Print out the best model : 
	cout << "************************\nThe best models : " << endl;
	cout << param.printCurrentModels(maxIndex[0]).at(1); 
	cout << param.printCurrentModels(maxIndex[1]).at(1);
	cout << param.printCurrentModels(maxIndex[2]).at(1); 
	cout << "************************\n" << endl;
 

}



