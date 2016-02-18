/*
 * Model.h
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#ifndef MODEL_H_
#define MODEL_H_
#include <string>
#include <vector>
#include "commons.h"
#include "Image.h"
#include <map>
#include <Eigen/Sparse>
#include <fstream>
#include <sstream>
#define NUM_PTMASS_PARAM 3
#define NUM_SIE_PARAM 5
#define NUM_NFW_PARAM 6
#define NUM_SPEMD_PARAM 6

typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::VectorXd vec;

using namespace std;

struct SingleModelParam {
	string name;
	double centerX, centerXFrom, centerXTo, centerXInc;
	double centerY, centerYFrom, centerYTo, centerYInc;
	double critRad, critRadFrom, critRadTo, critRadInc;
	double massScale, massScaleFrom, massScaleTo, massScaleInc; 
	double radScale, radScaleFrom, radScaleTo, radScaleInc; 

	double e, eFrom, eTo, eInc;
	double q, qFrom, qTo, qInc;
	double PA, PAFrom, PATo, PAInc;
	double power, powerFrom, powerTo, powerInc; 
	double mass;
	double core, coreFrom, coreTo, coreInc;
	double nParam;
};

class MultModelParam{


public:
	vector<SingleModelParam> parameter;
	int nLens ;
	vector<int> nParam;

	MultModelParam(map<string,string> confMap) {
		nLens = 0;
		map<string, string>::iterator itPTMASS = confMap.find("PTMASS") ;
		map<string, string>::iterator itSIE   = confMap.find("SIE");
		map<string, string>::iterator itNFW   = confMap.find("NFW");
		map<string, string>::iterator itSPEMD   = confMap.find("SPEMD");

		if(itPTMASS != confMap.end()) {
			vector<string> items = splitString(itPTMASS->second);
			SingleModelParam tempParam;
			tempParam.name = "PTMASS";
			tempParam.centerXFrom   = stof(items[0]);
			tempParam.centerXTo 	= stof(items[1]);
			tempParam.centerXInc 	= stof(items[2]);
			tempParam.centerYFrom 	= stof(items[3]);
			tempParam.centerYTo 	= stof(items[4]);
			tempParam.centerYInc 	= stof(items[5]);
			tempParam.critRadFrom	= stof(items[6]);
			tempParam.critRadTo   	= stof(items[7]);
			tempParam.critRadInc  	= stof(items[8]);

			parameter.push_back(tempParam);
			nParam.push_back(NUM_PTMASS_PARAM);
			nLens +=1;


		}
		if(itSIE != confMap.end()) {
			vector<string> items = splitString(itSIE->second);

			SingleModelParam tempParam;

			tempParam.name = "SIE";

			tempParam.centerXFrom 	= stof(items[0]);
			tempParam.centerXTo 	= stof(items[1]);
			tempParam.centerXInc 	= stof(items[2]);
			tempParam.centerYFrom 	= stof(items[3]);
			tempParam.centerYTo 	= stof(items[4]);
			tempParam.centerYInc 	= stof(items[5]);
			tempParam.critRadFrom 	= stof(items[6]);
			tempParam.critRadTo   	= stof(items[7]);
			tempParam.critRadInc  	= stof(items[8]);
			tempParam.eFrom			= stof(items[9]);
			tempParam.eTo 			= stof(items[10]);
			tempParam.eInc 			= stof(items[11]);
			tempParam.PAFrom 		= stof(items[12]);
			tempParam.PATo			= stof(items[13]);
			tempParam.PAInc			= stof(items[14]);
			parameter.push_back(tempParam);
			nParam.push_back(NUM_SIE_PARAM);

			nLens +=1;
		}
		if(itNFW != confMap.end()) {
			vector<string> items = splitString(itNFW->second);

			SingleModelParam tempParam;

			tempParam.name = "NFW";

			tempParam.centerXFrom 	= stof(items[0]);
			tempParam.centerXTo 	= stof(items[1]);
			tempParam.centerXInc 	= stof(items[2]);
			tempParam.centerYFrom 	= stof(items[3]);
			tempParam.centerYTo 	= stof(items[4]);
			tempParam.centerYInc 	= stof(items[5]);
			tempParam.massScaleFrom = stof(items[6]);
			tempParam.massScaleTo   = stof(items[7]);
			tempParam.massScaleInc  = stof(items[8]);
			tempParam.radScaleFrom 	= stof(items[9]);
			tempParam.radScaleTo   	= stof(items[10]);
			tempParam.radScaleInc  	= stof(items[11]);
			tempParam.eFrom			= stof(items[12]);
			tempParam.eTo 			= stof(items[13]);
			tempParam.eInc 			= stof(items[14]);
			tempParam.PAFrom 		= stof(items[15]);
			tempParam.PATo			= stof(items[16]);
			tempParam.PAInc			= stof(items[17]);
			parameter.push_back(tempParam);
			nParam.push_back(NUM_NFW_PARAM);

			nLens +=1;
		}

		if(itSPEMD != confMap.end()) {
			vector<string> items = splitString(itSPEMD->second);

			SingleModelParam tempParam;

			tempParam.name = "SPEMD";

			tempParam.centerXFrom 	= stof(items[0]);
			tempParam.centerXTo 	= stof(items[1]);
			tempParam.centerXInc 	= stof(items[2]);
			tempParam.centerYFrom 	= stof(items[3]);
			tempParam.centerYTo 	= stof(items[4]);
			tempParam.centerYInc 	= stof(items[5]);
			tempParam.critRadFrom 	= stof(items[6]);
			tempParam.critRadTo   	= stof(items[7]);
			tempParam.critRadInc  	= stof(items[8]);
			tempParam.eFrom			= stof(items[9]);
			tempParam.eTo 			= stof(items[10]);
			tempParam.eInc 			= stof(items[11]);
			tempParam.PAFrom 		= stof(items[12]);
			tempParam.PATo			= stof(items[13]);
			tempParam.PAInc			= stof(items[14]);
			tempParam.powerFrom		= stof(items[15]);
			tempParam.powerTo 		= stof(items[16]);
			tempParam.powerInc 		= stof(items[17]);
			tempParam.coreFrom 		= stof(items[18]);
			tempParam.coreTo		= stof(items[19]);
			tempParam.coreInc		= stof(items[20]);

			parameter.push_back(tempParam);
			nParam.push_back(NUM_SPEMD_PARAM);
			nLens +=1;
		}



	}

};

class Model {
	int length;


public:

	// For regularization;
	double chi2;
	double srcR;
	double phiR;
	double penalty;
	// Number of lens models;
	double nLens;
	double totalParam;
	MultModelParam param;

	vector<double> srcPosXList;	  // Source position after deflection in X direction, in arcsec;
	vector<double> srcPosYList;	  // Source position after deflection in Y direction, in arcsec;

	vector<double> srcPosXListPixel;	  // Source position after deflection in X direction, in pixel;
	vector<double> srcPosYListPixel;	  // Source position after deflection in Y direction, in pixel;

	vector<double> pDeltaX;  	// Deflection angle in X direction;
	vector<double> pDeltaY; 	// Deflection angle in Y direction;
	vector<double> invMag;
	vector<double> dSy1; 
	vector<double> dSy2; 

	vector<double> res_img;    // Residual brightness.
	vector<double> simple_res_img;


	vector<double> mod_img;
	vector<double> red_res_img; //
	vector<double> critical;
	vector<double> caustic;

	//vector<double> s;

	vector<double> something;
	sp_mat L;
	sp_mat M;
	vec r;
	vec new_r;
	vec s;
	vec phi;
	vec square_s;
	//
	sp_mat Ds;
	sp_mat Dphi;
	sp_mat Hs1;
	sp_mat Hs2;
	sp_mat Hphi;
	sp_mat HtH;
	sp_mat HphiH;
	sp_mat T;

	sp_mat RtR;

	sp_mat H0;  // zeroth-order regularization;
	sp_mat H1;  // gradient-order regularization;
	sp_mat H2;  // curvature-order regularization;
	double lambdaS;
	double lambdaPhi;


	vector<vector<normVec> > normV;
	vector<normVec> meanNormV;
	map<pair<int, int>,int> posMap;



public:
	Model();
	Model(Conf* conList, MultModelParam multModelParam, double lambdaS);
	void updateMatrixT(Conf* conf);
	vector<double> getDeflectionAngle(Conf* conList, int imgX, int imgY, double *pDeltaX, double *pDeltaY);
	void updatePosMapping(Image* image,  Conf* conList);
	void updateLensAndRegularMatrix(Image* dataImage,  Conf* constList);
	void updateGradient(Image* dataImage);
	void updateSrc(sp_mat* invC, vec d);
	void Logging(Image* dataImage, Conf* conList, string outFileName);
	void updateRegularMatrix();
	void updatePenalty( sp_mat* invC, vec d);
	void writeSrcImage(string outFileName, Conf* conList);
	void updateCritCaustic(Image* dataImage,  Conf* constList);
	virtual ~Model();
	void updateReducedResidual(Image* dataImage);
	double getRegularizationSrcValue (vec d);

	double getZerothOrderReg  	(Conf* conf, vector<double> briList);
	double getGradientOrderReg	(Conf* conf, vector<double> briList); 
	double getCurvatureOrderReg	(Conf* conf, vector<double> briList); 

	void clearMatrix();
};

#endif /* MODEL_H_ */
