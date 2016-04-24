/*
 * parafit.h
 *
 *  Created on: Dec 24, 2015
 *      Author: cheng109
 */

#ifndef PARAFIT_H_
#define PARAFIT_H_

#include "gsl/gsl_multimin.h"
#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
using namespace std;
//
typedef struct minimiser_params {
	Model* model;
	Image* dataImage;
	Conf* conf;
} minimiser_params;



//static double penalty_func(const gsl_vector *v, void *voidparams);
int	gsl_min_wrap(minimiser_params *params);
void gridSearch(Conf* conf, MultModelParam param, Image* dataImage, vec d, string dir, string outputFileName);
void mcFit(Conf* conf, MultModelParam param, Image* dataImage, vec d, string dir, string outputFileName);



#endif /* PARAFIT_H_ */



