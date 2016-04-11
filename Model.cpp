/*
 * Model.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: juncheng
 */

#include "Model.h"
#include "commons.h"
//#include <armadillo>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <cmath>
#include <fstream>
#include "fastell.h"
//#include "libiomp/omp.h"
//#include "gsl/gsl_multimin.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

using namespace std;
//using namespace arma;
using namespace Eigen;


/* Model::Model() {
	// TODO Auto-generated constructor stub
} */

Model::Model(Conf* conf, MultModelParam param, double lambdaS):
		length(conf->length),
		param(param),
		mod_img(length),
		red_res_img(length),
		L(length,length),
		r(2*length),
		new_r(2*length),
		s(length),
		phi(length),
		square_s(conf->srcSize[0]*conf->srcSize[1]),
		Ds(length, 2*length),
		Dphi(2*length,length),
		Hs1(length, length),
		Hs2(length, length),
		Hphi(length, length),
		HtH(length, length),
		HphiH(length, length),
		T(conf->srcSize[0]*conf->srcSize[1], length),
		RtR(2*length, 2*length),
		H0(conf->srcSize[0]*conf->srcSize[1], conf->srcSize[0]*conf->srcSize[1]),
		lambdaS(lambdaS)
		{
	// initial s;
	nLens = param.parameter.size();

	for(int i=0; i<nLens; ++i) {

		totalParam += param.nParam[i];
	}

	L.reserve(3*length);
	Ds.reserve(2*length);
	Dphi.reserve(2*length);
	Hs1.reserve(6*length);
	Hs2.reserve(6*length);
	RtR.reserve(200*length);
	H0.reserve(conf->srcSize[0]*conf->srcSize[1]);

	T.reserve(length);

	normVec n(0, 0, 0);
	vector<normVec> temp;
	temp.push_back(n);

	//lambdaS = 0.001; //0.01; //0.01;
	lambdaPhi = 10;
	for(int i=0; i<conf->length; ++i) {
		//s.push_back(0);
		normV.push_back(temp);
		dSy1.push_back(0);
		dSy2.push_back(0);
	}

	for(int i=0; i<conf->srcSize[0]*conf->srcSize[1]; ++i) {
			H0.insert(i, i)= 1;
	}

}

/***************************
Function:   	updateMatrixT
Description:  	(1) Member function of class Model;
				(2) T is matrix convert an adaptive grid to a regular grid;
				(3) T*S_adpative = S_regular ;
				(4) It will use the position of adaptive grid nodes:  srcPosXList, srcPosYList.
Arguments:		Conf*
Returns:		None
****************************/
void Model::updateMatrixT(Conf* conf) {

	// T size:  [srcSize[0]*srcSize[1],  2*length]
	int x, y, iList;
	for(int i=0; i<length; ++i) {
		x = nearbyint(srcPosXList[i]);
		y = nearbyint(srcPosYList[i]);

		if(x>0 && x< conf->srcSize[0] && y>0 && y<conf->srcSize[1]) {
			iList = conf->srcSize[0]*y+x;
			cout << iList << "\t ";
			T.insert(iList, i) = 1;
		}
	}

	cout << endl;
}

vector<double> Model::getDeflectionAngle(Conf* conf, int imgX, int imgY, double *pDeltaX, double *pDeltaY) {
	double fDenom = 0 ; 
	double srcX   = 0;  
	double srcY   = 0 ; 
	double fX  = 0; 
	double fY  = 0; 
	double pfX = 0; 
	double pfY = 0;
	vector<double> srcPos;

	

	pfX = (imgX - conf->imgXCenter ) * conf->imgRes; 			// image center frame;
	pfY = (imgY - conf->imgYCenter ) * conf->imgRes;


	//cout << "nLens: " << nLens << endl; 
	for(int i=0; i<nLens; ++i) {
		fX =  (imgX - conf->imgXCenter - param.parameter[i].centerX ) * conf->imgRes;   // lens center frame, 
		fY =  (imgY - conf->imgYCenter - param.parameter[i].centerY ) * conf->imgRes;
		// Unit:  aresecond.

		if(param.parameter[i].name.compare("PTMASS")==0) {
			fDenom = fX*fX+fY*fY;
			double fMult = param.parameter[i].critRad*param.parameter[i].critRad/fDenom;
			*pDeltaX +=  fX*fMult;
			*pDeltaY +=  fY*fMult;
		}


		if(param.parameter[i].name.compare("SIE")==0) {

			double phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;
			if (fX == 0 && fY == 0)  {
				*pDeltaX += param.parameter[i].critRad; 
				*pDeltaY += param.parameter[i].critRad;

			}
			

			//pre-calculate constants
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fq = 1-param.parameter[i].e;
			if (fq>1.0) cout << "Axis ratio should be smaller than 1. " << endl;
			if (fq==1.0) fq = 0.999;

			//rotate reference frame to x-axis
			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;

			root1mq = sqrt(1.0-fq*fq);
			phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
			fac = param.parameter[i].critRad*sqrt(fq)/root1mq;
			deltax1 = fac*atan(root1mq*x1/(phi + fCore));
			deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));

			//cout << root1mq << "\t" << y1 << "\t " << phi << "\t" << fq << "\t" << deltay1 << endl; 
			*pDeltaX += (deltax1*fCosTheta - deltay1*fSinTheta);
			*pDeltaY += (deltay1*fCosTheta + deltax1*fSinTheta);
			


		}
		

		if(param.parameter[i].name.compare("NFW")==0) {
		
			double fEllip,fCosTheta,fSinTheta,x1,y1,fPhi,fAngRadius,fTempResult,fCosPhi,fSinPhi,fScale;
  			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);

			fEllip = param.parameter[i].e; 
			fScale = param.parameter[i].radScale;

            if (fEllip >= 1.0 || fEllip < 0 || fScale < 0) 
                cout << "Bad parameters of 'e' or 'scale'. " << endl; 
                    
      		// create elliptical co-ords still in angle units from rotated frame sky coords 
			x1 = sqrt(1.0 - fEllip)*(fX*fCosTheta + fY*fSinTheta);
			y1 = sqrt(1.0 + fEllip)*(-fX*fSinTheta + fCosTheta*fY);
			fPhi = atan2(y1,x1);

			// angular radius is in dimensionless units 
			fAngRadius = sqrt(x1*x1 + y1*y1)/fScale;

			if (fAngRadius > 0.0) {
				double	deflx,defly;

				fCosPhi = cos(fPhi);
				fSinPhi = sin(fPhi);
					
				fTempResult = param.parameter[i].massScale * fScale * lm_nfw_mass(fAngRadius)/(fAngRadius);
				//cout << param.parameter[i].radScale << "\t" << fScale << "\t" << fAngRadius << "\t" << fTempResult  << endl; 
				deflx = sqrt(1.-fEllip)*fTempResult*fCosPhi;
				defly = sqrt(1.+fEllip)*fTempResult*fSinPhi;
				*pDeltaX += (deflx*fCosTheta - defly*fSinTheta);
				*pDeltaY += (deflx*fSinTheta + defly*fCosTheta);
			}
			else {
				*pDeltaX += 0.0;
				*pDeltaY += 0.0;
			} 
		}
		

		if(param.parameter[i].name.compare("SERSIC")==0)  {  
			double fCosTheta,fSinTheta,x1,y1,deflx,defly;
			double fEllip = param.parameter[i].e; 
            if (fEllip > 1.0 || fEllip <= 0) 
            	  cout << "Bad parameters of 'e' . " << endl; 
			
			/* rotate so that major axis of mass in x direction */
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);

			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;
			
			//iStatus = lm_deflCacheLookup(LM_SERSIC,&deVaucCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
			/*
			if (iStatus == LM_CACHE_MISS) {
					iStatus = lm_CalcSersicDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
					if (iStatus != 0) {
						goto EXIT;
					}
				}
			*/
				/* rotate back again */
			*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*param.parameter[i].kap;
			*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*param.parameter[i].kap;

			

		}



		if(param.parameter[i].name.compare("SPEMD")==0)  {  

			double	fTempKappa = 0.0,fTempCoreSqu=0.0, fTempAxratio=1.0, fTempDefl[2], fTempGamma = 0.0, fTempCenterX=0, fTempCenterY=0;
			double	x1, y1;
			double 	fCosTheta, fSinTheta;

            
			fTempAxratio = 1.0 - param.parameter[i].critRad;
			fTempCoreSqu = param.parameter[i].core*param.parameter[i].core;
			fTempGamma = param.parameter[i].power;
            if (fTempAxratio > 1.0 || fTempAxratio<=0 ) {
                cout << "Ellipticity is out of range; " << endl; 
            }
			if (fX == 0 && fY == 0 && fTempGamma >= 0.5) {
				*pDeltaX = *pDeltaY = param.parameter[i].critRad;
				//iStatus =LM_IGNORE_PROJ;
				//break;
			}
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);

                
            x1 = (fX*fCosTheta + fY*fSinTheta);
            y1 = (-fX*fSinTheta + fY*fCosTheta);
			fTempKappa = 0.5 * param.parameter[i].critRad * pow((2.0-2.0*fTempGamma)/fTempAxratio,fTempGamma);
			//extern  void fastelldefl_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2]);
			
			fastelldefl_(&x1,&y1,&fTempKappa,&fTempGamma,&fTempAxratio,&fTempCoreSqu,fTempDefl);

			*pDeltaX = fTempDefl[0]*fCosTheta - fTempDefl[1]*fSinTheta;
			*pDeltaY = fTempDefl[1]*fCosTheta + fTempDefl[0]*fSinTheta;  
			

			}
		}

	// SrcX and SrcY are in unit aresec;   
	// (SrcX, SrcY) = (0, 0) is the center point of source plane; 
	srcX = pfX - (*pDeltaX);     
	srcY = pfY - (*pDeltaY);


	/*   
	ofstream debug("debug.txt" , ios::out | ios::app  );
	
	debug << imgX << "\t" << imgY  << "\t" <<*pDeltaX<< "\t" << *pDeltaY << "\t"<<  srcX << "\t" << srcY << endl;
	debug.close();
	*/

	// In arcsecond
	srcPos.push_back(srcX);
	srcPos.push_back(srcY);
	return srcPos;

}


Model::~Model() {
	// TODO Auto-generated destructor stub
}



void Model::updatePosMapping(Image* image, Conf* conf) {

	length = conf->length;
	vector<double> srcPos;

	vector<double> test_srcPos; 
	//getDeflectionAngle(Conf* conf, int imgX, int imgY, double *pDeltaX, double *pDeltaY)

	double testDx =0; 
	double testDy = 0; 
	double testImgX = 136; 
	double testImgY = 136; 
	test_srcPos = getDeflectionAngle(conf, testImgX, testImgY, &testDx, &testDy); 

	//cout << testImgX << "\t" << testImgY << "\t" << test_srcPos[0]/conf->srcRes+conf->srcXCenter << "\t" 
	//<< test_srcPos[1]/conf->srcRes+conf->srcYCenter << "\t" << testDx << "\t" << testDy << "\t" <<endl; 

	for(int i=0; i<length; ++i) {
		int imgX = image->xList[i];
		int imgY = image->yList[i];
		//cout << param.critRad << endl;
		double defX = 0 ;
		double defY = 0 ;
		srcPos = getDeflectionAngle(conf,imgX, imgY, &defX, &defY);
		pDeltaX.push_back(defX);
		pDeltaY.push_back(defY);

		srcPosXList.push_back(srcPos[0]);
		srcPosYList.push_back(srcPos[1]);
		

		//cout << "dx: \t" << srcPos[0]/conf->srcRes
		 //<< "\tdy: " << srcPos[1]/conf->srcRes << endl; 


		srcPosXListPixel.push_back(srcPos[0]/conf->srcRes+conf->srcXCenter);
		srcPosYListPixel.push_back(srcPos[1]/conf->srcRes+conf->srcYCenter);

		posMap[make_pair(imgX, imgY)] = i;
	}

}

void Model::updateLensAndRegularMatrix(Image* dataImage,  Conf* constList) {

	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w, w5;


	for (int i=0; i<constList->length; ++i) {

	//	if(dataImage->type[i]==1) {// || dataImage->type[i]==0) {


			//L.insert(i,i)=1;

			left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
			right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
			up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
			down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

//			left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]-1));
//			right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]+1));
//			up    = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]+1));
//			down  = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]-1));

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {

				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;
				int iRight= right->second;
				Point A(srcPosXList[iLeft ], srcPosYList[iLeft ], s(iLeft ));
				Point B(srcPosXList[iUp   ], srcPosYList[iUp   ], s(iUp   ));
				Point C(srcPosXList[i     ], srcPosYList[i     ], s(i     ));
				Point D(srcPosXList[iDown ], srcPosYList[iDown ], s(iDown ));
				Point E(srcPosXList[iRight], srcPosYList[iRight], s(iRight));
				w5 = getPentWeigth(A, B, C, D, E);
				Hs1.insert(i, iLeft	) 	= w5[0];
				Hs1.insert(i, iUp	) 	= w5[1];
				Hs1.insert(i, i		) 	= w5[2];
				Hs1.insert(i, iDown	) 	= w5[3];
				Hs1.insert(i, iRight) 	= w5[4];

				Hs2.insert(i, iLeft	) 	= w5[5];
				Hs2.insert(i, iUp	) 	= w5[6];
				Hs2.insert(i, i		) 	= w5[7];
				Hs2.insert(i, iDown	) 	= w5[8];
				Hs2.insert(i, iRight) 	= w5[9];
			}


		//}
		if(dataImage->type[i]==1) {
			L.insert(i,i)=1;

		}
		if(dataImage->type[i]==0) {
			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;

				Point A(srcPosXList[iLeft], srcPosYList[iLeft], s(iLeft));
				Point B(srcPosXList[iUp  ], srcPosYList[iUp	 ], s(iUp  ));
				Point C(srcPosXList[iDown], srcPosYList[iDown], s(iDown));
				Point P(srcPosXList[i	 ], srcPosYList[i	 ], s(i    ));

				w = getTriWeight( A, B, C, P);
				L.insert(i, iLeft) 	= w[0];
				L.insert(i, iUp	 )  = w[1];
				L.insert(i, iDown) 	= w[2];

				normVec n = getNormVector(A, B, C);
				normV[iLeft].push_back(n);
				normV[iUp  ].push_back(n);
				normV[iDown].push_back(n);
			}
			else L.insert(i, i) = 1;
		}
	}



	HtH = Hs1.transpose()*Hs1 + Hs2.transpose()*Hs2;

	// test a new HtH
	//updateMatrixT(constList);
	//sp_mat uH = H0*T;
	//cout << H0.rows() << " " << H0.cols() << endl;
	//cout << T.rows() << " " << T.cols() << endl;

	//HtH = uH.transpose()*uH;
	//cout << "H0: " << H0 << endl;
	// end test



	//start = std::chrono::system_clock::now();
	//for (int k=0; k<HtH.outerSize(); ++k)
	//  for (SparseMatrix<double>::InnerIterator it(HtH,k); it; ++it)
	//	  RtR.insert(it.row(), it.col()) = it.value()*lambdaS*lambdaS;
	RtR = lambdaS*lambdaS*HtH;
	RtR.conservativeResize(2*constList->length, 2*constList->length);

}

void Model::updateGradient(Image* dataImage) {
	vec vSy1(length);  vSy1.fill(0);
	vec vSy2(length);  vSy2.fill(0);

	for (int i=0; i<length; ++i) {
		if(dataImage->type[i]==1) {
			normVec mean =  meanNormVector(normV[i]);
			if (mean.n2==0) {
				vSy1[i] = 0;
				vSy2[i] = 0;
			}
			else {
				vSy1[i] = -mean.n0/mean.n2;
				vSy2[i] = -mean.n1/mean.n2;
			}
		}
	}

	vSy1 = L*vSy1;
	vSy2 = L*vSy2;
	//update Ds and Dphi

	for(int i=0; i<length; ++i){
		Ds.coeffRef(i, 2*i+0) = vSy1(i);
		Ds.coeffRef(i, 2*i+1) = vSy2(i);
	}

	sp_mat extension = -1*Ds*Dphi;
	M.resize(L.rows(), L.cols()+extension.cols());
	M.middleCols(0,L.cols()) = L;
	M.middleCols(L.cols(), extension.cols()) = extension;


}



void Model::Logging(Image* dataImage, Conf* conList, string outFileName) {
	ofstream f(outFileName);
	string tab = "\t";
	string entry;
	f << "#1 index\n" << "#2 imgX\n" << "#3 imgY\n" << "#4 imgBri\n";
	f << "#5 srcX\n"  << "#6 srcY\n" << "#7 mapIndex" << "#8 ";
	for(int i=0; i<conList->length; ++i) {
		entry =  to_string(dataImage->iList[i]) + "\t"
				+to_string(dataImage->xList[i]) + "\t"
				+to_string(dataImage->yList[i]) + "\t"
				+to_string(dataImage->dataList[i]) + "\t"
				+to_string(srcPosXList[i]) + "\t"
				+to_string(srcPosYList[i]) + "\t"
				+to_string(posMap[make_pair(dataImage->xList[i], dataImage->yList[i])]) + "\t" ;  //  Matched position
				//+to_string()
		f << entry << endl ;
	}
	f.close();
}

void Model::updateSrc(sp_mat* invC, vec d) {

}

void Model::updatePenalty(sp_mat* invC, vec d) {

	// problem:  Ax=b
	sp_mat A = M.transpose()*(*invC)*M + RtR;

	vec b = M.transpose()*(*invC)*d;
	// sp_mat A1 = A.submat(0, 0, length-1, length-1);
	sp_mat A1 = A.block(0,0, length, length);
	vec b1 = b.head(length);

	// Using Eigen methods;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(A1);
	s = chol.solve(b1);

	// get Penalty Function;
	for(int i=0; i<length; ++i) {
		r(i) = s(i);
		new_r(i) = d(i);
	}

	for(int i=length; i<2*length; ++i) {
		r(i) = 0;
		new_r(i) = 0;
	}
	vec res = M*r-d;
	res_img = eigenV_to_cV(res);
	simple_res_img = eigenV_to_cV(M*new_r-d);

	mod_img = eigenV_to_cV(M*new_r);

	vec temp1 =  res.transpose()*(*invC)*res;
	//cout << "HtH: " << HtH << endl ;
	double temp2 = s.transpose()*HtH*s;
	//cout << "source regularization: " << temp2 << endl;
	srcR =  lambdaS*lambdaS*temp2;
	chi2 = temp1(0,0);
	penalty = chi2 + srcR;


}



void Model::updateCritCaustic(Image* dataImage,  Conf* conf) {

	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w, w5;

	vector<double> a11, a12, a21, a22;
	double h = conf->imgRes;
	for (int i=0; i<conf->length; ++i) {
		left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
		right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
		up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
		down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

		if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {

			int iLeft = left->second;
			int iUp   = up  ->second;
			int iDown = down->second;
			int iRight= right->second;

			a11.push_back((pDeltaX[iRight]-pDeltaX[iLeft])/(2*h));
			a12.push_back((pDeltaX[iUp]-pDeltaX[iDown])/(2*h));
			a21.push_back((pDeltaY[iRight]-pDeltaY[iLeft])/(2*h));
			a22.push_back((pDeltaY[iUp]-pDeltaY[iDown])/(2*h));
		}
		else {
			a11.push_back(0);
			a12.push_back(0);
			a21.push_back(0);
			a22.push_back(0);
		}
		invMag.push_back(1.0-(a11[i]+a22[i])+a11[i]*a22[i]-a12[i]*a21[i]);
	}
	// update inverse magnification
	int sign_t = 0;

	string header = "# Region file format: DS9 version 4.1\n\
					# Filename: horseshoe_test/HorseShoe_new.fits\n\
					global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n\
					image\n"; 
	ofstream criticalFile; 
	criticalFile.open(conf->criticalName); 
	criticalFile << header; 

	for (int i=0; i<conf->length; ++i) {
			left  = posMap.find(make_pair(dataImage->xList[i]-1, dataImage->yList[i]));
			right = posMap.find(make_pair(dataImage->xList[i]+1, dataImage->yList[i]));
			up    = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]+1));
			down  = posMap.find(make_pair(dataImage->xList[i], dataImage->yList[i]-1));

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {

				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;
				int iRight= right->second;
				sign_t = sign(invMag[i])*(sign(invMag[iLeft]) + sign(invMag[iRight]) + sign(invMag[iUp]) + sign(invMag[iDown]));
			}
			else
				sign_t = 5;  // Assign a value bigger than 4;


			double distSqure = (dataImage->xList[i]-150)*(dataImage->xList[i]-150) + (dataImage->yList[i]-150)*(dataImage->yList[i]-150); 

			if(sign_t<4 && distSqure > 50*50) {
				// Distance from center; 
				
				critical.push_back(1);
				criticalFile << "point(" << to_string(dataImage->xList[i]) << "," << to_string(dataImage->yList[i]) << ")\n" ; 
			}
			else
				critical.push_back(0);
	}

	criticalFile.close(); 
}

void Model::updateReducedResidual(Image* dataImage) {

	for(int i=0; i<length; ++i) {
		red_res_img[i] = res_img[i]/dataImage->varList[i];
	}




}

void Model::writeSrcImage(string outFileName, Conf* conList) {
	vector<double> sBright = eigenV_to_cV(s);

	Image* srcImg = new Image(srcPosXListPixel, srcPosYListPixel, &sBright, conList->srcSize[0], conList->srcSize[1], conList->bitpix);
	srcImg->writeToFile(outFileName);
	delete srcImg;

}


double Model::getRegularizationSrcValue (vec d) {
	// HtH and d is known
	return  d.transpose()*HtH*d;
}


double Model::getZerothOrderReg (Conf* conf, vector<double> briList) {
	//s is known;
	double sum = 0;
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];
	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		for (int j=0; j< naxis2 ; ++j) {
			int index = i+j*naxis1;
			sum += srcImg->data[index] * srcImg->data[index] ;
		}
	}
	delete srcImg;
	return sum ;
}


double Model::getGradientOrderReg(Conf* conf, vector<double> briList) {
	double sum = 0;
	double diff = 0;  
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];
	int index = 0 ; 
	int index_edge = 0; 
	int index_next = 0; 
	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		index_edge = (naxis2 - 1)* naxis1 + i; 
		sum +=  srcImg->data[index_edge] * srcImg->data[index_edge]; 		
		for (int j=0; j< naxis2-1 ; ++j) {
			index  		= i + j 	* naxis1;
			index_next 	= i + (j+1) * naxis1; 
			diff = srcImg->data[index] - srcImg->data[index_next]; 
			sum += diff*diff; 
		}
	}

	for (int j=0; j< naxis2 ; ++j) {
		index_edge = j * naxis2 + (naxis1-1) ; 
		sum +=  srcImg->data[index_edge] *  srcImg->data[index_edge]; 
		for (int i=0; i< naxis1-1 ; ++i) {
			index  		= i 	+ j * naxis1;
			index_next 	= (i+1) + j * naxis1; 
			diff = srcImg->data[index] - srcImg->data[index_next]; 
			sum += diff*diff; 
		}
	}
	delete srcImg; 
	return sum ; 
}


double Model::getCurvatureOrderReg(Conf* conf, vector<double> briList) {

	double sum = 0;
	double diff = 0;  
	double diff_edge = 0; 
	long naxis1 = conf->srcSize[0];
	long naxis2 = conf->srcSize[1];

	int index1 = 0 ;
	int index2 = 0 ;
	int index3 = 0 ;

	int index_edge1 = 0; 
	int index_edge2 = 0;	 

	Image* srcImg =new Image(srcPosXListPixel, srcPosYListPixel, &briList, naxis1, naxis2, conf->bitpix);

	for (int i=0; i< naxis1 ; ++i) {
		index_edge1 = (naxis2 - 2)* naxis1 + i;
		index_edge2 = (naxis2 - 1)* naxis1 + i; 
		diff_edge = srcImg->data[index_edge1] -  srcImg->data[index_edge2]; 
		sum +=  diff_edge * diff_edge +  srcImg->data[index_edge2]*  srcImg->data[index_edge2]; 

		for ( int j=0 ; j< naxis2-2; ++j) {
			index1 = i + j  * naxis1; 
			index2 = index1 + naxis1; 
			index3 = index2 + naxis1; 
			diff = srcImg->data[index1] - 2* srcImg->data[index2] + srcImg->data[index3]; 
			sum += diff * diff; 

		}

	}	
	
	for ( int j=0 ; j< naxis2; ++j) {
		index_edge1 = j * naxis1 + naxis1-2;
		index_edge2 = j * naxis1 + naxis1-1; 
		diff_edge = srcImg->data[index_edge1] -  srcImg->data[index_edge2]; 
		sum +=  diff_edge * diff_edge +  srcImg->data[index_edge2]*  srcImg->data[index_edge2];

		for (int i=0; i< naxis1-2 ; ++i) {
			index1 = i + j  * naxis1; 
			index2 = index1 + 1; 
			index3 = index2 + 1; 
			diff = srcImg->data[index1] - 2* srcImg->data[index2] + srcImg->data[index3]; 
			sum += diff * diff; 

		}

	}	

	delete srcImg; 
	return sum ; 
}

MultModelParam::MultModelParam(map<string,string> confMap) {
		

		nLens = 0;
		map<string, string>::iterator itPTMASS	= confMap.find("PTMASS") ;
		map<string, string>::iterator itSIE    	= confMap.find("SIE");
		map<string, string>::iterator itNFW   	= confMap.find("NFW");
		map<string, string>::iterator itSPEMD 	= confMap.find("SPEMD");
		map<string, string>::iterator itSERSIC 	= confMap.find("SERSIC");

		if(itSERSIC != confMap.end()) {
			vector<string> items = splitString(itSERSIC->second);

			SingleModelParam tempParam;

			tempParam.name = "SERSIC";

			tempParam.centerXFrom 		= stof(items[0]);
			tempParam.centerXTo 		= stof(items[1]);
			tempParam.centerXInc 		= stof(items[2]);
			tempParam.centerYFrom 		= stof(items[3]);
			tempParam.centerYTo 		= stof(items[4]);
			tempParam.centerYInc 		= stof(items[5]);
			tempParam.kapFrom 			= stof(items[6]);
			tempParam.kapTo   			= stof(items[7]);
			tempParam.kapInc  			= stof(items[8]);
			tempParam.eFrom 			= stof(items[9]);
			tempParam.eTo   			= stof(items[10]);
			tempParam.eInc  			= stof(items[11]);
			tempParam.PAFrom			= stof(items[12]);
			tempParam.PATo 				= stof(items[13]);
			tempParam.PAInc 			= stof(items[14]);
			tempParam.sersicScaleFrom 	= stof(items[15]);
			tempParam.sersicScaleTo		= stof(items[16]);
			tempParam.sersicScaleInc	= stof(items[17]);
			tempParam.mFrom 			= stof(items[18]);
			tempParam.mTo   			= stof(items[19]);
			tempParam.mInc  			= stof(items[20]);
			
			parameter.push_back(tempParam);
			nParam.push_back(NUM_SERSIC_PARAM);

			nLens +=1;
		}



		if(itPTMASS != confMap.end()) {

			vector<string> strs;

			std::string s = itPTMASS->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
				break; 
			}; 
			strs.push_back(s); 

			for(int i=0; i<strs.size(); ++i) {
				vector<string> items = splitString(strs[i]);
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

		}
		if(itSIE != confMap.end()) {
			
			vector<string> strs;
			std::string s = itSIE->second; 
			string delimiter = "_&&_";

			size_t pos = s.find(delimiter) ; 
			while( pos!=std::string::npos) {
				strs.push_back(s.substr(0, pos)); 
				s = s.substr(pos+4); 
				pos = s.find(delimiter);
			}; 
			strs.push_back(s); 
			for(int i=0; i<strs.size(); ++i) {

				vector<string> items = splitString(strs[i]);
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


void MultModelParam::printModels() {
	/*
	string name;
	double mass;
	double nParam;
	// All shared parameters: 
	double centerX, centerXFrom, centerXTo, centerXInc;
	double centerY, centerYFrom, centerYTo, centerYInc;
	
	double e, eFrom, eTo, eInc;
	double q, qFrom, qTo, qInc;
	double PA, PAFrom, PATo, PAInc;

	// For PTMASS model and SIE model: 
	double critRad, critRadFrom, critRadTo, critRadInc;
	double power, powerFrom, powerTo, powerInc; 
	double core, coreFrom, coreTo, coreInc;
	// For NFW model
	double massScale, massScaleFrom, massScaleTo, massScaleInc;   
	double radScale, radScaleFrom, radScaleTo, radScaleInc; 	
	// For sersic model: 
	double kap, kapFrom, kapTo, kapInc; 
	double sersicScale, sersicScaleFrom, sersicScaleTo, sersicScaleInc; 
	double m, mFrom, mTo, mInc; 

	int nLens ;
	vector<int> nParam;

	*/

	cout << "\n*********** Models *********" << endl;
	cout << "nLens: 	 " << nLens << endl; 
	for(int i=0; i<nLens; ++i) {
		if (parameter[i].name=="PTMASS") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].critRadFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< endl; 
				
		}

		if (parameter[i].name=="SIE") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].critRadFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t"
				<< endl; 
				
		}

		if (parameter[i].name=="NFW") {

			cout << "[" << parameter[i].name  << "]:" 
				<< "\n_From\t"
				<< parameter[i].centerXFrom << "\t" 
				<< parameter[i].centerYFrom << "\t"
				<< parameter[i].massScaleFrom << "\t"
				<< parameter[i].radScaleFrom << "\t"
				<< parameter[i].eFrom << "\t" 
				<< parameter[i].PAFrom << "\t" 
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].massScaleTo << "\t"
				<< parameter[i].radScaleTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].massScaleInc << "\t"
				<< parameter[i].radScaleInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t" 
				<< endl; 
		}
	}
	cout << "*******************************" << endl;
}


void MultModelParam::mix() { 
	
	
	/*	mixModel structure:    {name, paraList[0, 0, 0, 0, 0, 0, 0, 0]}
	*
	*/

	vector<vector<mixModels> > mix; 
	for(int i=0; i<nLens; ++i) {
		if (parameter[i].name=="SIE") {
			vector<mixModels> v1; 
			for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
				for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
					for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
						for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
							for (double PA = parameter[i].PAFrom;	PA <= parameter[i].eTo; PA += parameter[i].PAInc) {
								mixModels sModel("SIE");
								
								sModel.paraList[0] = critRad; 
								sModel.paraList[1] = centerX; 
								sModel.paraList[2] = centerY; 
								sModel.paraList[3] = e;
								sModel.paraList[4] = PA;
								sModel.paraList[5] = 0;  //parameter[i].core;
								v1.push_back(sModel); 
							}
						}
					}
				}
			}
		
			mix.push_back(v1); 
		}
	}
	mix.resize(3); 
	/* for maximum 3 models:  j, k, m */

		for(int j=0; j<mix[0].size(); ++j) {
			for(int k=0; k<mix[1].size(); ++k ) {
				for(int m=0; m<mix[2].size(); ++m) {
					vector<mixModels> v2; 
					v2.push_back(mix[0][j]); 
					v2.push_back(mix[1][k]); 
					v2.push_back(mix[2][m]); 
					mixAllModels.push_back(v2); 
				}

			}

		}
	nComb = mixAllModels.size(); 
}


void MultModelParam::printCurrentModels(int curr) {
		
	//cout <<  mixAllModels[curr].size() << endl; 
	cout << "[" << curr+1 << "/" <<  nComb << "]" << endl; 
	for(int i=0; i<mixAllModels[curr].size(); ++i) {
		cout << mixAllModels[curr][i].name << ":\t" ; 
		for (int j=0; j< 8; ++j) {
			cout << mixAllModels[curr][i].paraList[j] << "\t" ; 
		}
		cout << "\n"; 
	} 
	cout << endl; 
}


#if 0
#endif 

