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
	occupation  = 0 ; 

	/*

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
	*/
	

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

vector<double> Model::getDeflectionAngle(Conf* conf, double pfX, double pfY, double *pDeltaX, double *pDeltaY, MultModelParam * param) {
	double fDenom = 0 ; 
	double srcX   = 0;  
	double srcY   = 0 ; 
	double fX  = 0; 
	double fY  = 0; 
	//double pfX = 0; 
	//double pfY = 0;
	vector<double> srcPos;
	int nLens = param->parameter.size(); 
	
	//cout << "nLens: " << nLens << endl; 
	for(int i=0; i<nLens; ++i) {
		fX =  pfX - (param->parameter[i].centerX * conf->imgRes);   // lens center frame, 
		fY =  pfY - (param->parameter[i].centerY * conf->imgRes);
		// Unit:  aresecond.

		if(param->parameter[i].name.compare("PTMASS")==0) {
			fDenom = fX*fX+fY*fY;
			double fMult = param->parameter[i].critRad*param->parameter[i].critRad/fDenom;
			*pDeltaX +=  fX*fMult;
			*pDeltaY +=  fY*fMult;
		}


		if(param->parameter[i].name.compare("SIE")==0) {

			double phi,root1mq,fq,fac,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;

			double fCore= param->parameter[i].core; 
			if (fX == 0 && fY == 0)  {
				*pDeltaX += param->parameter[i].critRad; 
				*pDeltaY += param->parameter[i].critRad;

			}
			

			//pre-calculate constants
			fCosTheta = cos(param->parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fSinTheta = sin(param->parameter[i].PA*M_PI/180 + 0.5* M_PI);
			fq = 1-param->parameter[i].e;
			if (fq>1.0) cout << "Axis ratio should be smaller than 1. " << endl;
			if (fq==1.0) fq = 0.999;

			//rotate reference frame to x-axis
			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;

			root1mq = sqrt(1.0-fq*fq);
			phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
			fac = param->parameter[i].critRad*sqrt(fq)/root1mq;
			deltax1 = fac*atan(root1mq*x1/(phi + fCore));
			deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));

			//cout << root1mq << "\t" << y1 << "\t " << phi << "\t" << fq << "\t" << deltay1 << endl; 
			*pDeltaX += (deltax1*fCosTheta - deltay1*fSinTheta);
			*pDeltaY += (deltay1*fCosTheta + deltax1*fSinTheta);
			


		}
		
		/*
		if(param->parameter[i].name.compare("NFW")==0) {
		
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
			
			// rotate so that major axis of mass in x direction 
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);

			x1 = fX*fCosTheta + fY*fSinTheta;
			y1 = -fX*fSinTheta + fY*fCosTheta;
			
			//iStatus = lm_deflCacheLookup(LM_SERSIC,&deVaucCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
			
			if (iStatus == LM_CACHE_MISS) {
					iStatus = lm_CalcSersicDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
					if (iStatus != 0) {
						goto EXIT;
					}
				}
			
				
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
		*/
		}

		
	// SrcX and SrcY are in unit aresec;   
	// (SrcX, SrcY) = (0, 0) is the center point of source plane; 
	srcX = pfX - (*pDeltaX);     
	srcY = pfY - (*pDeltaY);
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

	
	for(int i=0; i<length; ++i) {
		int imgX = image->xList[i];
		int imgY = image->yList[i];
		//cout << param.critRad << endl;
		double defX = 0 ;
		double defY = 0 ;

		double pfX = (imgX - conf->imgXCenter ) * conf->imgRes; 			// image center frame;
		double pfY = (imgY - conf->imgYCenter ) * conf->imgRes;

		srcPos = getDeflectionAngle(conf,pfX, pfY, &defX, &defY, &param);


		pDeltaX.push_back(defX);
		pDeltaY.push_back(defY);

		srcPosXList.push_back(srcPos[0]);
		srcPosYList.push_back(srcPos[1]);
		
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





Image Model::getFullResidual(Image* dataImage) {
	// Assume "mod_image" is known; 
	mod_img = dataImage->dataList; 
	Image fullResidualImage(* dataImage); 
	map<pair<int, int>, double> modMap; 
	for(size_t i=0; i< mod_img.size(); ++i)  {
		int pos = dataImage->yList[i] * dataImage->naxis1 + dataImage->xList[i]; 
		fullResidualImage.data[pos] -= mod_img[i]; 
	}

	return fullResidualImage; 	

	/*
	for(int i=0; i<length; ++i) {
		red_res_img[i] = res_img[i]/dataImage->varList[i];
	} */




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


double Model::getScatterReg() {

	// Now we know the 'srcPosXListPixel' and 'srcPosYListPixel': 
	double scatter = 0;   // STD of position 'x' and 'y'
	double sumX = 0; 
	double sumY = 0;  
    int counter(0);
	for (size_t i=0; i<srcPosXList.size(); ++i) {
        if (std::isnan(srcPosXList[i]) || std::isnan(srcPosYList[i])) continue;
		sumX += srcPosXList[i]; 
		sumY += srcPosYList[i];
        counter++;
	}
    cout<<"sum "<<sumX<<" "<<sumY<<" "<<counter<<endl;

	double xPosMean = sumX / counter; 
	double yPosMean = sumY / counter;
	for (size_t i=0; i<srcPosXList.size(); ++i) {
        if (std::isnan(srcPosXList[i]) || std::isnan(srcPosYList[i])) continue;
		scatter += (srcPosXList[i]-xPosMean) * (srcPosXList[i]-xPosMean) ;
		scatter += (srcPosYList[i]-yPosMean) * (srcPosYList[i]-yPosMean) ;
	}
	scatter = scatter / counter;  
	return scatter; 
	
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

			double val = srcImg->data[index] * srcImg->data[index]; 
			if(val!=0) {
				sum += val;
				occupation += 1; 
			}

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

			for(size_t i=0; i<strs.size(); ++i) {
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
			for(size_t i=0; i<strs.size(); ++i) {

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
				tempParam.coreFrom		= stof(items[15]); 
				tempParam.coreTo		= stof(items[16]); 
				tempParam.coreInc		= stof(items[17]); 

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
				<< parameter[i].coreFrom << "\t"
				<< "\n_To\t"
				<< parameter[i].centerXTo<< "\t" 
				<< parameter[i].centerYTo << "\t"
				<< parameter[i].critRadTo << "\t"
				<< parameter[i].eTo << "\t" 
				<< parameter[i].PATo << "\t" 
				<< parameter[i].coreTo << "\t"
				<< "\n_Inc\t"
				<< parameter[i].centerXInc << "\t" 
				<< parameter[i].centerYInc << "\t"
				<< parameter[i].critRadInc << "\t"
				<< parameter[i].eInc << "\t" 
				<< parameter[i].PAInc << "\t"
				<< parameter[i].coreInc << "\t"
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

		if (parameter[i].name=="PTMASS") {
			vector<mixModels> v1; 
			for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
				for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
					for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
						
								mixModels sModel("PTMASS");
								
								sModel.paraList[0] = critRad; 
								sModel.paraList[1] = centerX; 
								sModel.paraList[2] = centerY; 
								sModel.paraList[3] = 0;
								sModel.paraList[4] = 0;
								sModel.paraList[5] = 0;  //parameter[i].core;
								v1.push_back(sModel); 
					}
				}
			}
		
			mix.push_back(v1); 
		}

		if (parameter[i].name=="SIE") {

			vector<mixModels> v1; 
			for (double critRad = parameter[i].critRadFrom;	critRad <= parameter[i].critRadTo; critRad += parameter[i].critRadInc) {
				for (double centerX = parameter[i].centerXFrom;	centerX <= parameter[i].centerXTo; centerX += parameter[i].centerXInc) {
					for (double centerY = parameter[i].centerYFrom;	centerY <= parameter[i].centerYTo; centerY += parameter[i].centerYInc) {
						for (double e = parameter[i].eFrom;	e <= parameter[i].eTo; e += parameter[i].eInc) {
							for (double PA = parameter[i].PAFrom; PA <= parameter[i].PATo; PA += parameter[i].PAInc) {
								for(double core = parameter[i].coreFrom; core <= parameter[i].coreTo; core += parameter[i].coreInc) {
									mixModels sModel("SIE");
									sModel.paraList[0] = critRad; 
									sModel.paraList[1] = centerX; 
									sModel.paraList[2] = centerY; 
									sModel.paraList[3] = e;
									sModel.paraList[4] = PA;
									sModel.paraList[5] = core;  //parameter[i].core;
									v1.push_back(sModel); 
								}
							}
						}
					}
				}
			}
		
			mix.push_back(v1); 
		}



	}
    cout << "mix size: " << mix.size() << endl;
	if (mix.size() > 3) mix.resize(3);

	/* for maximum 3 models:  j, k, m */

		for(size_t j=0; j<mix[0].size(); ++j) {
			for(size_t k=0; k<mix[1].size(); ++k ) {
				for(size_t m=0; m<mix[2].size(); ++m) {
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


vector<string> MultModelParam::printCurrentModels(int curr) {
		
	//cout <<  mixAllModels[curr].size() << endl; 
	vector<string> ret; 
	string modelsInRow; 
	string modelsInCol; 

	// return a string for "output.txt"; 
	modelsInCol += ( "[" + to_string(curr+1)+"/"+to_string(nComb)  + "]\n" ); 
	for(size_t i=0; i<mixAllModels[curr].size(); ++i) {
		modelsInCol += (mixAllModels[curr][i].name  + ":\t") ; 
		for (int j=0; j< 8; ++j) {
			modelsInCol += (to_string(mixAllModels[curr][i].paraList[j]) + "\t") ; 
			modelsInRow += (to_string(mixAllModels[curr][i].paraList[j]) + "\t")  ; 

		}
		modelsInCol += "\n";  
	} 
	modelsInCol += "\n";
	//modelsInRow += "\n"; 
	ret.push_back(modelsInRow); 
	ret.push_back(modelsInCol); 

	return ret; 
}


void Model::clearVectors() {

	/*vector<double> srcPosXListPixel;	  // Source position after deflection in X direction, in pixel;
	vector<double> srcPosYListPixel;	  // Source position after deflection in Y direction, in pixel;

	vector<double> pDeltaX;  	// Deflection angle in X direction;
	vector<double> pDeltaY; 	// Deflection angle in Y direction;
	vector<double> invMag;
	vector<double> dSy1; 
	vector<double> dSy2; 

	vector<double> res_img;    // Residual brightness
	vector<double> res_full_img;  //  the full residual map; 
	vector<double> simple_res_img;


	vector<double> mod_img;
	vector<double> critical;
	*/



	param.parameter.clear(); 
	srcPosXListPixel.clear(); 
	srcPosYListPixel.clear(); 
	srcPosXList.clear(); 
	srcPosYList.clear(); 
	pDeltaX.clear(); 
	pDeltaY.clear(); 
	critical.clear();

	res_img.clear(); 
	res_full_img.clear(); 
	simple_res_img.clear(); 
	mod_img.clear(); 


}


vector<vector<double> > getCritCausticFine(vector<double> xPosListArc, vector<double> yPosListArc, Conf* conf, MultModelParam * param, int level) {
	// level = 5; 
	vector<double> invMag;
	map<pair<int, int>,int> posMap;
	vector<double> srcPos, pDeltaX, pDeltaY; 
	vector<double> srcPosXListArc; 
	vector<double> srcPosYListArc; 

	

	double curr_res = conf->imgRes/level; 

	vector<vector<double> > ret; 

	vector<double> imgXList; 
	vector<double> imgYList; 

	double side = 1.0/level;   // unit is pixel; 

	int index = 0 ; 
	for (size_t i=0; i<xPosListArc.size(); ++i) {
		int old_imgX = int(round(xPosListArc[i] / conf->imgRes + conf->imgXCenter )); 
		int old_imgY = int(round(yPosListArc[i] / conf->imgRes + conf->imgYCenter )); 

		for(int j=0; j<level; ++j) {		
			for(int k=0; k<level; ++k) {
			double new_posX = (j * side) * curr_res + xPosListArc[i];   // in arcsecond
			double new_posY = (k * side) * curr_res + yPosListArc[i];   // in arcsecond

			double defX = 0 ;
			double defY = 0 ;  
			srcPos = Model::getDeflectionAngle(conf, new_posX, new_posY, &defX, &defY, param);
			pDeltaX.push_back(defX);    // in arcsec; 
			pDeltaY.push_back(defY);		
			srcPosXListArc.push_back(srcPos[0]); 
			srcPosYListArc.push_back(srcPos[1]);

			int imgX = old_imgX * level + k; 
			int imgY = old_imgY * level + j; 
			imgXList.push_back(imgX); 
			imgYList.push_back(imgY); 
			posMap[make_pair(imgX, imgY)] = index;
			++index; 

			}
		}
	}		

	map<pair<int, int>,int>::iterator left, right, up, down;
	vector<double> w, w5;

	vector<double> a11, a12, a21, a22;
	vector<double> new_xPosListArc; 
	vector<double> new_yPosListArc; 
	vector<double> new_xSrcPosListArc; 
	vector<double> new_ySrcPosListArc; 


	double h = curr_res;
	for (int i=0; i< index; ++i) {
		left  = posMap.find(make_pair(imgXList[i]-1, imgYList[i]));
		right = posMap.find(make_pair(imgXList[i]+1, imgYList[i]));
		up    = posMap.find(make_pair(imgXList[i], imgYList[i]+1));
		down  = posMap.find(make_pair(imgXList[i], imgYList[i]-1));

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
	

	for (int i=0; i<index; ++i) {
			
			left  = posMap.find(make_pair(imgXList[i]-1, imgYList[i]));
			right = posMap.find(make_pair(imgXList[i]+1, imgYList[i]));
			up    = posMap.find(make_pair(imgXList[i], imgYList[i]+1));
			down  = posMap.find(make_pair(imgXList[i], imgYList[i]-1));

			if(left!=posMap.end() && up!=posMap.end() && down!=posMap.end() && right!=posMap.end()) {
				int iLeft = left->second;
				int iUp   = up  ->second;
				int iDown = down->second;
				int iRight= right->second;
				sign_t = sign(invMag[i])*(sign(invMag[iLeft]) + sign(invMag[iRight]) + sign(invMag[iUp]) + sign(invMag[iDown]));
	
			}
			else
				sign_t = 5;  // Assign a value bigger than 4;

			if(sign_t<4 ) { //} && distSqure > 50*50) {
				new_xPosListArc.push_back((imgXList[i]-conf->imgXCenter*level)* curr_res); 
				new_yPosListArc.push_back((imgYList[i]-conf->imgYCenter*level)* curr_res);
				new_xSrcPosListArc.push_back(srcPosXListArc[i]); 
				new_ySrcPosListArc.push_back(srcPosYListArc[i]); 


			}
	}

	// Center:  (0, 0)  in arcsecond
	ret.push_back(new_xPosListArc); 
	ret.push_back(new_yPosListArc); 
	ret.push_back(new_xSrcPosListArc); 
	ret.push_back(new_ySrcPosListArc); 
	return ret; 
}


vector<Image* > getCritCaustic(Conf* conf, MultModelParam * param) {
	// automatic using  full Image ( no region files ); 
	
	vector<Image*> ret; 
	vector<double> critical, xList, yList;
	vector<double> caustic, srcXList, srcYList; 
	vector<double> srcPosXListPixel; 
	vector<double> srcPosYListPixel; 
	int level =  conf->causticLevel; 

	Image* dataImage = new Image(conf->imageFileName); 
	dataImage->updateFilterImage("whatever..", 0);   // without any region file --- using all data points; 
	int length = dataImage->data.size();	

	vector<double> xPosListArc; 
	vector<double> yPosListArc; 

	for(int i=0; i<length; ++i) {
		xPosListArc.push_back((dataImage->xList[i]-conf->imgXCenter)*conf->imgRes); 
		yPosListArc.push_back((dataImage->yList[i]-conf->imgYCenter)*conf->imgRes); 
	}

	// Basic critical searching; 
	vector<vector<double> >  critXY  = getCritCausticFine(xPosListArc, yPosListArc, conf, param, 1);
	// Finer critical searching with higher level; 
	if (level > 1 ) 
		critXY = getCritCausticFine(critXY[0], critXY[1], conf, param, level); 

	double newResolution = conf->imgRes/level; 
	double newSrcResolution = conf->srcRes/level; 

	for(size_t i=0; i<critXY[0].size(); ++i) {


		// Image coordinate (in pixel)
		xList.push_back(critXY[0][i]/newResolution + conf->imgXCenter*level +1); 
		yList.push_back(critXY[1][i]/newResolution + conf->imgYCenter*level +1); 

		srcXList.push_back(critXY[2][i]/newSrcResolution + conf->srcXCenter*level +1); 
		srcYList.push_back(critXY[3][i]/newSrcResolution + conf->srcYCenter*level +1); 
		critical.push_back(1); 

	}
	createDs9Contour(&xList,    &yList,    level, conf->contourCritName) ; 
	createDs9Contour(&srcXList, &srcYList, level, conf->contourCausName) ; 

	/// modify ending: 
	Image* critImg = new Image(xList, yList, &critical, conf->imgSize[0]*level, conf->imgSize[1]*level, conf->bitpix);
	Image* causImg = new Image(srcXList, srcYList, &critical, conf->srcSize[0]*level, conf->srcSize[1]*level, conf->bitpix);

	ret.push_back(critImg); 
	ret.push_back(causImg); 
 	return ret; 
}


void createDs9Contour(vector<double>* xList, vector<double>* yList, double level, string contourFileName) {
	//   contour in image coordinate; 
	ofstream contourFile(contourFileName);
	map<pair<int, int>, int> posMap;
	for (size_t i=0; i<xList->size(); ++i) {
		posMap[make_pair(int((*xList)[i]), int((*yList)[i]))] = i; 
	}
	vector <map<pair<int, int>,int>::iterator> n(9); //  n1, n2, n3, n4, n5, n6, n7, n8, n9; 
	/* find all the neighbors of each points; 
	 *  n0, || n1,  n2
	 * 	n3, || n4,  n5
	 *	n6, || n7,  n8
	 */
	for (size_t i=0; i<xList->size(); ++i) {
			int x = int((*xList)[i]); 
			int y = int((*yList)[i]); 

			n[0] = posMap.find(make_pair(x-1, y+1));
			n[1] = posMap.find(make_pair(x  , y+1));
			n[2] = posMap.find(make_pair(x+1, y+1));
			n[3] = posMap.find(make_pair(x-1, y  ));
			n[4] = posMap.find(make_pair(x  , y  ));    // self; 
			n[5] = posMap.find(make_pair(x+1, y  ));
			n[6] = posMap.find(make_pair(x-1, y-1));
			n[7] = posMap.find(make_pair(x  , y-1));
			n[8] = posMap.find(make_pair(x+1, y-1));


			for (int j=0; j<9; ++j) {
				if(n[j]!=posMap.end()) {
					int ind = n[j]->second; 
					contourFile << to_string((*xList)[i]/level)   << "\t" << to_string((*yList)[i]/level)   << endl; 
 					contourFile << to_string((*xList)[ind]/level) << "\t" << to_string((*yList)[ind]/level) << endl; 
 					contourFile << endl; 
				}

			}

		}
 	contourFile.close(); 
}



Image* createLensImage(Conf* conf, MultModelParam * param) {

	double level = 1; 
	vector<int> xList;  // in pixel, the same size as 'dataImage'; 
	vector<int> yList;
	vector<double> val;
	for(int i=0; i<conf->imgSize[1]*level; ++i) {
		for(int j=0; j<conf->imgSize[0]*level; ++j) {
			xList.push_back(j); 
			yList.push_back(i); 
			val.push_back(0); 
		}
	}


	cout << param->nLens << endl; 
	for(int i=0; i< param->nLens; ++i) {

		if(param->parameter[i].name =="PTMASS") {
		}
		if (param->parameter[i].name =="SIE") {

			cout << i << endl; 
			double centerX = param->parameter[i].centerX;
			double centerY = param->parameter[i].centerY; 
			double critRad = param->parameter[i].critRad; 
			double q = 1-param->parameter[i].e; 
			double PA = param->parameter[i].PA/180 * M_PI + 0.5 * M_PI; 
			double core = param->parameter[i].core; 
			if (q==1)  q=0.999; 
			for(size_t j=0; j<xList.size(); ++j) {
				double x = (xList[j] - (centerX + conf->imgXCenter))*conf->imgRes  ; 
				double y = (yList[j] - (centerY + conf->imgYCenter))*conf->imgRes  ; 
				double new_x = x*cos(PA) + y*sin(PA); 
				double new_y = x*sin(PA) - y*cos(PA);  
				double coeff =1 ; // sqrt(q/(1-q*q)); 
				val[j] += critRad*coeff /sqrt(new_x*new_x*q + new_y*new_y/q + core*core); 
			}


		}
		if(param->parameter[i].name =="NFW") {
		}


		

	}



	cout << "val: " << xList[0] << "\t" << val[0] << endl; 
	Image* lensImg = new Image(xList, yList, &val, conf->imgSize[0]*level, conf->imgSize[1]*level, conf->bitpix);

	return lensImg; 


}


#if 0
#endif 






