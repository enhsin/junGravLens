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


int main() {
	/***** prepare*****/
	//string dir("pt_test/");
	//string dir("nfw_test/"); 
	//string dir("sersic_test/"); 
	string dir("sie_test/");
	//string dir("horseshoe_test/");cd
	//string dir("spemd_test/"); 
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

	//cout << "start " << endl; 
	MultModelParam param = MultModelParam(mapConf);
	param.printModels();

	
	cout << param.nLens << endl;


	


	double lambdaS = 0.001;

	minimiser_params *min_params = new minimiser_params();

	//params->model = model1;
	min_params->dataImage = dataImage;
	min_params->conf = conf;
	min_params->model = new Model(conf, param, lambdaS); 


	gridSearch(conf, param,  dataImage, d, dir, lambdaS);	


	//int status = gsl_min_wrap(min_params);

	//int test = Ameoba_test();
	delete conf; 
	delete dataImage; 
	delete min_params->model; 
	delete min_params; 
	return 0;
}





