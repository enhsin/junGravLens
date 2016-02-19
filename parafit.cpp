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
typedef struct minimiser_params {
	Model* model;
	Image* dataImage;
	Conf* conf;
} minimiser_params;




static double penalty_func(const gsl_vector *v, void *voidparams) {

	minimiser_params *params = new minimiser_params();
	params = (minimiser_params *) voidparams;




	Image* dataImage = params->dataImage;
	Conf* conf = params->conf;
	Model* model=params->model;

	Model* localModel = new Model(params->conf, model->param, model->lambdaS);

	int iParam = 0;
	for(int i=0; i<localModel->nLens; ++i) {
		for(int j=0; j<localModel->param.parameter[i].nParam; ++j) {
			localModel->param.parameter[i].critRad = gsl_vector_get(v,iParam);
			iParam ++;
		}
	}

	vec d =dataImage->getMatrixD();
	localModel->updatePosMapping(dataImage, conf);
	localModel->updateLensAndRegularMatrix(dataImage, conf);
	localModel->updateGradient(dataImage);
	localModel->updatePenalty(&dataImage->invC, d);

	
	cout << "iteration" << endl;  
	delete localModel;
	return params->model->penalty;
}


/***************************
Function:   gsl_min_wrap
Description:    Wrapper function that sets up the GSL amoeba style minimiser then calls the minimiser.
Arguments:
Returns:    0: success.
****************************/
int	gsl_min_wrap(minimiser_params *min_params) {

	cout << "hello0 " << endl;
	size_t iter = 0;
	int status=0,iParam=0,iCount=0,iOffset=0; 
	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *x=NULL;
	gsl_vector *step=NULL;
	gsl_multimin_function my_func;
	double fSize=0;
	double fTemp=0;


	cout << "hello1 " << endl;

	my_func.f = &penalty_func;
	my_func.n = 1; //params->model->totalParam;
	my_func.params = min_params;

	int nParam = 1; 
	iParam=0;


	//params->model = new Model(params->conf, model->param, model->lambdaS);

	/* Set up initial guess. Just take the middle of the param range */
	x = gsl_vector_alloc(nParam);


	//for (int j=0; j< params->model->nLens; j++) {
	//	for (int i=0; i< params->model->param.parameter[j].nParam ; i++) {
	fTemp =  7.0; // (params->model->params.parameter[0].critRadFrom + params->model->params.parameter[0].critRadTo)*0.5;
	gsl_vector_set (x, iParam,fTemp);
	//			cout << "number of param: " <<  iParam << endl;
	//		iParam++;
	//	}
	//}

	/* set step sizes */
	//step = gsl_vector_alloc(params->model->totalParam);
	step = gsl_vector_alloc(nParam);
	double stepSize = 0.1 ; //params->model->parameter[0].critRadInc; 
	gsl_vector_set (step, iParam,stepSize);

	//for (int j=0; j< params->model->nLens; j++) {
	//	for (int i=0; i< params->model->param.parameter[j].nParam ; i++) {
	//		fTemp = 0.5; //(params->pLens->pComponent[j].fParamto[i] - params->pLens->pComponent[j].fParamfrom[i])*0.5;
	
	//		iParam++;
	//	}
	//}

	T = gsl_multimin_fminimizer_nmsimplex;
	s = gsl_multimin_fminimizer_alloc(T, min_params->model->totalParam);
	s = gsl_multimin_fminimizer_alloc(T, nParam);	

	gsl_multimin_fminimizer_set(s, &my_func, x, step);

	cout << "hello2 " << endl; 

	do {
		iter++;

		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;
		if (status == GSL_SUCCESS)
			printf ("\nMinimum found at:\n");



		fSize = gsl_multimin_fminimizer_size (s);
		cout << iter << "\t" << gsl_vector_get(s->x, 0) << "\t" << s->fval << "\t" << fSize <<  endl;
		status = gsl_multimin_test_size (fSize, 1e-6);
	} while (status == GSL_CONTINUE && iter < 1000);


	if (s!=NULL) gsl_multimin_fminimizer_free(s);
	if (x!=NULL) gsl_vector_free(x);
	if (step!=NULL) gsl_vector_free(step);

	return 0;
}



void gridSearch(Conf* conf, MultModelParam param, Image* dataImage, vec d, string dir, double lambdaS) {



  //string dir("pt_test/");

	//double lambdaS = 0.01;
	Model *model1 = new Model(conf, param, lambdaS);

	int i=0;

	// PTMASS model:
	if(param.parameter[0].name=="PTMASS") {
		for (param.parameter[0].critRad=param.parameter[0].critRadFrom ;
				param.parameter[0].critRad <= param.parameter[0].critRadTo;
				param.parameter[0].critRad += param.parameter[0].critRadInc) {   // elliptictiy // Critical radius
			model1 = new Model(conf, param, lambdaS);
			model1->updatePosMapping(dataImage, conf);
			model1->updateLensAndRegularMatrix(dataImage, conf);
			model1->updateGradient(dataImage);
			model1->updatePenalty(&dataImage->invC, d);
			// Write residual image:
			model1->updateReducedResidual(dataImage);

			//vector<double> sBright = eigenV_to_cV(model1->s);
			vector<double> sBright = dataImage->dataList; 
			//double srcReg = model1->getRegularizationSrcValue(d);
			double zerothOrder 		= model1->getZerothOrderReg   (conf, sBright);
			double gradientOrder 	= model1->getGradientOrderReg (conf, sBright); 
			double curvatureOrder 	= model1->getCurvatureOrderReg(conf, sBright);
			++i;


			
			Image* resImg   	= 	new Image(dataImage->xList, dataImage->yList, &model1->res_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			//Image* simple_resImg= 	new Image(dataImage->xList, dataImage->yList, &model1->simple_res_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			Image* modelImg 	= 	new Image(dataImage->xList, dataImage->yList, &model1->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			//Image* srcImg 		= 	new Image(model1->srcPosXList, model1->srcPosYList, &dataImage->dataList, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
			Image* srcImg 		= 	new Image(model1->srcPosXListPixel, model1->srcPosYListPixel, &sBright, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
			resImg		 ->writeToFile  (dir + "img_res_"		+to_string(i) +".fits");
			//simple_resImg->writeToFile  (dir + "img_simple_res_"+to_string(i) +".fits");
			modelImg	 ->writeToFile	(dir + "img_mod_"		+to_string(i) +".fits");
			srcImg 		 -> writeToFile (dir + "img_src_" 		+to_string(i) +".fits");

			//model1	->writeSrcImage	(dir + "img_src_" + to_string(i) +".fits", conf);

			cout << model1->param.parameter[0].critRad  << "\t" ;
			cout << zerothOrder << "\t" ;
			cout << gradientOrder << "\t" ;
			cout << curvatureOrder << "\t"; 
			cout << model1->chi2 << "\t" << model1->srcR << "\t" << model1->penalty ; 
			cout << endl;
			delete model1;

		}
	}

	// SIE model:

	if(param.parameter[0].name=="SIE") {
		ofstream output; 
		output.open("output.txt"); 



		for (param.parameter[0].critRad=param.parameter[0].critRadFrom ;
			param.parameter[0].critRad <=param.parameter[0].critRadTo;
			param.parameter[0].critRad += param.parameter[0].critRadInc) {   // elliptictiy // Critical radius

		 	for (param.parameter[0].e=param.parameter[0].eFrom ; 
		  		param.parameter[0].e <=param.parameter[0].eTo;
		  		param.parameter[0].e += param.parameter[0].eInc) {
			  	
			  	for (param.parameter[0].PA=param.parameter[0].PAFrom ;
			  		param.parameter[0].PA <=param.parameter[0].PATo;
			  		param.parameter[0].PA += param.parameter[0].PAInc) {   

			  		for (param.parameter[0].centerX=param.parameter[0].centerXFrom ;
			  			param.parameter[0].centerX <=param.parameter[0].centerXTo;
			  			param.parameter[0].centerX += param.parameter[0].centerXInc) {  

			  			for (param.parameter[0].centerY=param.parameter[0].centerYFrom ;
			  				param.parameter[0].centerY <=param.parameter[0].centerYTo;
			  				param.parameter[0].centerY += param.parameter[0].centerYInc) {   
				  			//param.parameter[0].PA = 90; 
					  		model1 = new Model(conf, param, lambdaS);
							model1->updatePosMapping(dataImage, conf);
							model1->updateLensAndRegularMatrix(dataImage, conf);
							model1->updateGradient(dataImage);
							model1->updatePenalty(&dataImage->invC, d);

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


							resImg	->writeToFile  	(dir + "img_res_" + to_string(i) + "_" + to_string(j) +".fits");
							modelImg->writeToFile	(dir + "img_mod_" + to_string(i) + "_" + to_string(j) +".fits");
							//model1	->writeSrcImage	(dir + "img_src_" + to_string(i) + "_" + to_string(j) +".fits", conf);
							srcImg -> writeToFile (dir + "img_src_" + to_string(i) + "_" + to_string(j) +".fits");

							// output to console; 
							cout << model1->param.parameter[0].centerX  << "\t" ;
							cout << model1->param.parameter[0].centerY  << "\t" ;
							cout << model1->param.parameter[0].critRad  << "\t" ;
							cout << model1->param.parameter[0].e  << "\t" ;
							cout << model1->param.parameter[0].PA  << "\t" ;
					

							//cout << vegetiiReg << "\t" ; 
							cout << zerothOrder << "\t" ;
							cout << gradientOrder << "\t" ;
							cout << curvatureOrder << "\t"; 
							//cout << model1->chi2 << "\t" << model1->srcR << "\t" << model1->penalty ; 
							cout << endl;

							// output to file; 
							output << model1->param.parameter[0].centerX  << "\t" ;
							output << model1->param.parameter[0].centerY  << "\t" ;
							output << model1->param.parameter[0].critRad  << "\t" ;
							output << model1->param.parameter[0].e  << "\t" ;
							output << model1->param.parameter[0].PA  << "\t" ;
					

							//output << vegetiiReg << "\t" ; 
							output << zerothOrder << "\t" ;
							output << gradientOrder << "\t" ;
							output << curvatureOrder << "\t"; 
							//output << model1->chi2 << "\t" << model1->srcR << "\t" << model1->penalty ; 
							output << endl;
					
							delete model1;
							}
						}
					}
				}
			}


		output.close(); 
		
	}





}



