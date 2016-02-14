#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <ctime>
#include "Image.h"
#include "commons.h"
#include "fitsio.h"
#include "Model.h"
#include "parafit.h"
#include "gsl/gsl_multimin.h"
//#include <armadillo>
#include <tuple>
#include <map>
//#include <boost/python.hpp>
#include "gl_crit.h"
//#include "gl_crit.h"

#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

using namespace std;
//using namespace arma;





/*
int test() {

	**** prepare****
	string dir("pt_test/");
	//string dir("horseshoe/");
	//string dir("sie_test/");
	//string dir("blind_test/");
	string confFileName = dir+ "conf.txt";
	map<string, string> mapConf = parseConfigure(confFileName);
	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], 1);
	dataImage->updateGridPointType();
	dataImage->updateVarList(1.5, 0.25); // (threshold, var);
	dataImage->invC = dataImage->getVarMatrix();
	Conf *conf = new Conf(dataImage, mapConf);


	vec d =dataImage->getMatrixD();
	conf->printConfList();


	vector<ParaList> p;
	p.push_back(ParaList("SIE", 0, 0, 16.41, 0.4, 167, 0.0, 0.0));
	//p.push_back(ParaList("PTMASS", 30, 0, 2.6212, 0.0, 0.0, 0.0, 0.0));
	Model *model = new Model(conf, p, 1.0);
	model->updatePosMapping(dataImage, conf);  // update pDeltaX, pDeltaY;
	model->updateCritCaustic(dataImage, conf);

	Image* critImg = new Image(dataImage->xList, dataImage->yList, &model->critical, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
	critImg->writeToFile(dir + "critCurve.fits");
	Image* causticImg = new Image(model->srcPosXList, model->srcPosYList, &model->critical, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
	causticImg->writeToFile( dir + "causticCurve.fits");

	******************************************************************************
	**** routine to manipulate an image ****
//	Image* inputImage = new Image("PSF_test/cut_psf.fits");
//	inputImage->updateFilterImage("none", 0);
//	for(int i=0; i<inputImage->length; ++i) {
//		if (inputImage->xList[i]==10 && inputImage->yList[i] == 10) {
//			inputImage->dataList[i] = 1;
//		}
//		else inputImage->dataList[i] = 0;
//	}
//	inputImage->writeFilterImage("outputImage.fits");
	******************************************************************************
	**** routine to generate critical and caustic curve  ****
//	vector<ParaList> p;
//	p.push_back(ParaList("PTMASS", 0, 30, 4.3070, 0.0, 0.0, 0.0, 0.0));
//	p.push_back(ParaList("PTMASS", 30, 0, 2.6212, 0.0, 0.0, 0.0, 0.0));
//	Model *model = new Model(conf, p);
//	model->updatePosMapping(dataImage, conf);  // update pDeltaX, pDeltaY;
//
//	model->updateCritCaustic(dataImage, conf);
//
//	Image* critImg = new Image(dataImage->xList, dataImage->yList, &model->critical, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
//	critImg->writeToFile("critCurve.fits");
//
//
//	Image* causticImg = new Image(model->srcPosXList, model->srcPosYList, &model->critical, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
//	causticImg->writeToFile("causticCurve.fits");
	******************************************************************************
//	Image* psfTestImage= new Image("outputImage.fits");
//	psfTestImage->updateFilterImage("none", 0);
//
//	sp_mat blurOperator = generatePSFoperator(mapConf["psfFileName"],psfTestImage);
//	psfTestImage->printImageInfo(0, 10, 0, 10);
//
//	vec ss = cV_to_eigenV(psfTestImage->data);
//	vec afterBlur = blurOperator*ss;
//	vector<double> ss2= eigenV_to_cV(afterBlur);
//
//	Eigen::MatrixXd dMat = Eigen::MatrixXd(blurOperator);
//	// plot all ss;
//
//	double sum = 0 ;
//	for(int i=0; i<psfTestImage->data.size(); ++i) {
//		sum += ss2[i];
//		//sum += psfTestImage->data[i];
//	}
//
//	cout << "After blur sum: " << sum << endl ;
//	Image* bluredImage = new Image(psfTestImage->xList, psfTestImage->yList, &ss2, psfTestImage->naxis1, psfTestImage->naxis2, conf->bitpix);
//	bluredImage->writeToFile("PSF_test/image_after_psf.fits");
	******************************************************************************



	//ParaList para("PTMASS", 0, 0, 5.6826, 0.0, 0.0, 0.0, 0.0);

	//vector<double> srcBriList(conf->length, 0);



	Model *model = new Model(conf, para);
	model->updatePosMapping(dataImage, conf);
	model->updateLensAndRegularMatrix(dataImage, conf);
	model->updateGradient(dataImage);
	model->updateSrc(&invC, d);
	vector<double> src= model-> getDeflectionAngle(conf, 35, 52);
	cout << src[0] << endl;
	cout << src[1] << endl;


	//Image* srcImg = new Image(model->srcPosXList, model->srcPosYList, &srcBriList, conf->srcSize[0], conf->srcSize[1], conf->bitpix);
	//srcImg->writeToFile();

	model->writeSrcImage( "src.fits", conf);



	double minPenalty = 10e9;
	double minR = 10e9;
	//chrono::time_point<std::chrono::system_clock> start, end;

	vector<ParaList> paraList;
	paraList.push_back(ParaList("SIE", 0, 0, 0, 0.0, 0.0, 0.0, 0.0));
	double lambdaS = 1.0;
	Model *model1 = new Model(conf, paraList, lambdaS); // = new Model(conf, paraList); //


	for (int i=0; i<20;  ++i) {   // elliptictiy // Critical radius
		for (int j=0; j<1; ++j) {   // PA;
			//start = std::chrono::system_clock::now();
			//paraList[0].critRad = 5.61+i*0.05;  //5.41 + 0.05*i;
			paraList[0].critRad = 5.46 + 0.05*i ;
			//paraList[0].e = 0.0 + i*0.002;
			//paraList[0].PA = 135; // 0 + j*5.0;
			model1 = new Model(conf, paraList, lambdaS);
			model1->updatePosMapping(dataImage, conf);
			model1->updateLensAndRegularMatrix(dataImage, conf);
			model1->updateGradient(dataImage);
			model1->updatePenalty(&dataImage->invC, d);
			// Write residual image:
			model1->updateReducedResidual(dataImage);


			Image* resImg = new Image(dataImage->xList, dataImage->yList, &model1->res_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			Image* modelImg = new Image(dataImage->xList, dataImage->yList, &model1->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			resImg->writeToFile  (dir + "img_res_"+to_string(i) + "_" +to_string(j) +".fits");
			modelImg->writeToFile(dir + "img_mod_"+to_string(i) + "_" +to_string(j) +".fits");
			model1->writeSrcImage(dir + "img_src_"+to_string(i) + "_" +to_string(j) +".fits", conf);

			if(model1->penalty < minPenalty)  {
				minPenalty = model1->penalty;
				minR = model1->param[0].critRad;
			}
			cout << model1->param[0].critRad  << "\t" ;
			cout << model1->param[0].e  << "\t" ;
			cout << model1->param[0].PA << "\t" <<  model1->chi2 << "\t" << model1->srcR << "\t" << model1->penalty << endl;
			delete model1;
		}
	}
	return 0;
}

*/

// Ameoba test:

/* Paraboloid centered on (p[0],p[1]), with
   scale factors (p[2],p[3]) and minimum p[4] */

double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *p = (double *)params;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return p[2] * (x - p[0]) * (x - p[0]) +
           p[3] * (y - p[1]) * (y - p[1]) + p[4];
}

/* The gradient of f, df = (df/dx, df/dy). */


/* Compute both f and df together. */


int Ameoba_test() {
	double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};

	  const gsl_multimin_fminimizer_type *T =
	    gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc (2);
	  gsl_vector_set (x, 0, 5.0);
	  gsl_vector_set (x, 1, 7.0);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (2);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 2;
	  minex_func.f = my_f;
	  minex_func.params = par;

	  s = gsl_multimin_fminimizer_alloc (T, 2);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate(s);

	      if (status)
	        break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 1e-2);

	      if (status == GSL_SUCCESS)
	        {
	          printf ("converged to minimum at\n");
	        }

	      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
	              iter,
	              gsl_vector_get (s->x, 0),
	              gsl_vector_get (s->x, 1),
	              s->fval, size);
	    }
	  while (status == GSL_CONTINUE && iter < 100);

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free (s);

	  return status;

}


int main() {
	/***** prepare*****/
	string dir("pt_test/");
	//string dir("horseshoe/");
	//string dir("sie_test/");
	//string dir("blind_test/");
	string confFileName = dir+ "conf.txt";
	map<string, string> mapConf = parseConfigure(confFileName);


	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], 1);
	dataImage->updateGridPointType();
	dataImage->updateVarList(1, 0.1); // (threshold, var);
	dataImage->invC = dataImage->getVarMatrix();
	Conf *conf = new Conf(dataImage, mapConf);


	vec d =dataImage->getMatrixD();
	conf->printConfList();

	MultModelParam param = MultModelParam(mapConf);



	minimiser_params *params = new minimiser_params();

	//params->model = model1;
	params->dataImage = dataImage;
	params->conf = conf;


	double lambdaS = 1.0;

	gridSearch(conf, param,  dataImage, d, dir, lambdaS);





	//int status = gsl_min_wrap(params);

	//int test = Ameoba_test();

	return 0;
}




