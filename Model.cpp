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
	double fDenom, srcX, srcY, fX, fY, pfX, pfY;
	vector<double> srcPos;
	for(int i=0; i<nLens; ++i) {
		// Unit:  aresecond.
		fX = (imgX -conf->imgXCenter-param.parameter[i].centerX)*conf->imgRes;   // lens center frame, 
		fY = (imgY -conf->imgYCenter-param.parameter[i].centerY)*conf->imgRes;

		pfX = (imgX-conf->imgXCenter)*conf->imgRes; 			// image center frame;
		pfY = (imgY-conf->imgYCenter)*conf->imgRes;

		//double pDeltaX = 0;
		//double pDeltaY = 0;

		if(param.parameter[i].name.compare("PTMASS")==0) {
			fDenom = fX*fX+fY*fY;
			double fMult = param.parameter[i].critRad*param.parameter[i].critRad/fDenom;
			*pDeltaX +=  fX*fMult;
			*pDeltaY +=  fY*fMult;


		}


		if(param.parameter[i].name.compare("SIE")==0) {

			double phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;
			//	fX -= g_PixelResn*pLensComp->fParameter[3];
			//	fY -= g_PixelResn*pLensComp->fParameter[4];
			
			if (fX == 0 && fY == 0)
				*pDeltaX = *pDeltaY = param.parameter[i].critRad; // pLensComp->fParameter[0];

			//pre-calculate constants
			fCosTheta = cos(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
			fSinTheta = sin(param.parameter[i].PA*M_PI/180 + 0.5*M_PI);
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


			*pDeltaX = deltax1*fCosTheta - deltay1*fSinTheta;
			*pDeltaY = deltay1*fCosTheta + deltax1*fSinTheta;


		}


	}
	// SrcX and SrcY are in unit aresec;   
	// (SrcX, SrcY) = (0, 0) is the center point of source plane; 
	srcX = pfX - (*pDeltaX);     
	srcY = pfY - (*pDeltaY);

	ofstream debug("debug.txt" , ios::out | ios::app  );

	debug << imgX << "\t" << imgY  << "\t" <<*pDeltaX<< "\t" << *pDeltaY << "\t"<<  srcX << "\t" << srcY << endl;
	debug.close();
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
		srcPos = getDeflectionAngle(conf,imgX, imgY, &defX, &defY);
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

	mod_img = eigenV_to_cV(M*r);

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
	double h = 0.2; //conf->res;
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

			if(sign_t<4)
				critical.push_back(1);
			else
				critical.push_back(0);
	}
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


	return sum ; 
}


