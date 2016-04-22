/*
 * commons.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: cheng109
 */
#include "commons.h"
#include "fitsio.h"
#include <string>
#include "Image.h"
//#include "Model.h"
//#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <Eigen/Sparse>
#include <map>
#include <ctime>
#include <chrono>


typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::VectorXd vec;

#define buffsize 1000
using namespace std;


Conf::Conf(Image* dataImage, map<string, string> confMap) {
		// Cosmological constants:

		omega = stod(confMap["omega"]);
		lambda = stod(confMap["lambda"]);
		weos = stod(confMap["weos"]);
		hubble = stod(confMap["hubble"]);
		srcZ = stod(confMap["srcZ"]);
		lenZ = stod(confMap["lenZ"]);
		imageFileName = confMap["imageFileName"]; 

		long naxis1, naxis2, len ;
		int bit;
		double res;
		dataImage->getConstants(&len, &naxis1, &naxis2, &res, &bit);
		res = 0.3;
		//cout << "len: " << len << endl;
		imgSize[0]=naxis1; imgSize[1] =naxis2;
		potSize[0]=naxis1; potSize[1] =naxis2;


		imgXCenter = naxis1/2.0;
		imgYCenter = naxis2/2.0;
		potXCenter = naxis1/2.0;
		potYCenter = naxis2/2.0 ;
		length = dataImage->length;
		bitpix = bit;

		criticalName = confMap["criticalName"]; 
		causticName = confMap["causticName"]; 
		contourCritName = confMap["contourCritName"]; 
		contourCausName = confMap["contourCausName"]; 

		srcRes = stod(confMap["srcRes"]);
		imgRes = stod(confMap["imgRes"]);
		potRes = stod(confMap["potRes"]);

		srcSize[0] =stod(confMap["srcX"]);
		srcSize[1] =stod(confMap["srcY"]);

		back_mean = stod(confMap["back_mean"]); 
		back_std = stod(confMap["back_std"]);



		usingRegion   = stoi(confMap["usingRegion"]); 
		outputSrcImg  = stoi(confMap["outputSrcImg"]); 
		outputModImg  = stoi(confMap["outputModImg"]); 
		outputCritImg = stoi(confMap["outputCritImg"]); 
		outputLensImg = stoi(confMap["outputLensImg"]); 
		srcBackground = stoi(confMap["srcBackground"]); 	

		causticLevel  = stoi(confMap["causticLevel"]); 

		srcXCenter = srcSize[0]/2.0;
		srcYCenter = srcSize[1]/2.0;
}

void Conf::printConfList(){
		cout << "*********** Constants *********" << endl;
		cout << "srcSize:    " << srcSize[0] << ",\t"<< srcSize[1] << endl;
		cout << "imgSize:    " << imgSize[0] << ",\t"<<imgSize[1]  << endl;
		cout << "potSize:    " << potSize[0] << ",\t"<<potSize[1]  << endl;
		cout << "srcRes:     " << srcRes << endl;
		cout << "imgRes:     " << imgRes << endl;
		cout << "potRes:     " << potRes << endl;
		cout << "srcCenter:  " << srcXCenter << ",\t"<<srcYCenter << endl;
		cout << "imgCenter:  " << imgXCenter << ",\t"<<imgYCenter << endl;
		cout << "potCenter:  " << potXCenter << ",\t"<<potYCenter << endl;
		cout << "length:     " << length << endl;
		cout << "*******************************" << endl;
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

bool pnpoly(size_t nvert, vector<double> *vertx, vector<double> *verty, double testx, double testy)
{
  int i, j;
  bool c = FALSE;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty->at(i)>testy) != (verty->at(j)>testy)) &&
     (testx < (vertx->at(j)-vertx->at(i)) * (testy-verty->at(i)) / (verty->at(j)-verty->at(i)) + vertx->at(i)) )
       c = !c;
  }
  return c;
}


string parseReagionFile(string regionFileName, vector<double> *xpos, vector<double> *ypos) {
	// regionType = 0:   not supported yet; 
	// regionType = 1:   Polygon region;
	// regionType = 2:   Point region; 

	ifstream regionFile(regionFileName.c_str());
	string line, token;
	size_t pos=0;
	string regionType; 
	size_t pos1 = 0; 
	size_t pos2 = 0; 
	while (getline(regionFile, line)) {
		if (line[0]!='#' && line.substr(0, 6)!="global" && line.substr(0, 5)!="image") {
			// For polygon region file
			if(line.substr(0,7)=="polygon") {
				pos1 = line.find_first_of("(", pos);
				pos2 = line.find_first_of(")", pos);

				istringstream ss(line.substr(pos1+1, pos2-pos1));
				int flag=-1;
				while(getline(ss, token, ',' )){
					if(flag<0) xpos->push_back(stod(token));
					if(flag>0) ypos->push_back(stod(token));
					//cout << token << endl;
					flag = (-1)*flag;
				};
				regionType = "polygon"; 
			}	
			if(line.substr(0,5)=="point") { 
				pos1 = line.find_first_of("(", pos);
				pos2 = line.find_first_of(")", pos);

				istringstream ss(line.substr(pos1+1, pos2-pos1));	
				getline(ss, token, ','); 
				xpos->push_back(stod(token));
				getline(ss, token, ')'); 
				ypos->push_back(stod(token)); 
				//cout << token << endl; 
				regionType = "point"; 
			}

			if(line.substr(0,3)=="box") {
				pos1 = line.find_first_of("(", pos);
				pos2 = line.find_first_of(")", pos);
				istringstream ss(line.substr(pos1+1, pos2-pos1));	
				for (int j=0;j<4 ; ++j) {
					getline(ss, token, ','); 
					xpos->push_back(stod(token));
					ypos->push_back(stod(token));
				}
				getline(ss, token, ')'); 
				xpos->push_back(stod(token)); 
				ypos->push_back(stod(token));

				regionType = "box"; 
			}


			if(line.substr(0,6)=="circle") {
				pos1 = line.find_first_of("(", pos);
				pos2 = line.find_first_of(")", pos);
				istringstream ss(line.substr(pos1+1, pos2-pos1));	
				for (int j=0;j<2 ; ++j) {
					getline(ss, token, ','); 
					xpos->push_back(stod(token));
					ypos->push_back(stod(token));
				}
				getline(ss, token, ')'); 
				xpos->push_back(stod(token)); 
				ypos->push_back(stod(token));

				regionType = "circle"; 
			}
		}

		if(xpos->size()!=ypos->size())
			cout << "Error when reading region File!" << endl;
	}
	return regionType; 

}


map<string,string> parseConfigure(string confFileName) {

	map<string, string> confMap;
	ifstream confFile(confFileName.c_str());
	string line;
	string modelString;
	while (getline(confFile, line)) {
		if(line[0]!='#') {
			vector<string> items = splitString(line);
			if (items.size()==2)
				confMap[items[0]] = items[1];
			else if (items.size() >2) {
				modelString.clear();
				for(int i=1; i<items.size(); ++i) {
					modelString += ( items[i] + "\t") ;
				}
				if (confMap.find(items[0])!= confMap.end()) {
					string newAdd = "_&&_" + modelString; 
					confMap[items[0]] += newAdd; 
				} 
				else 
					confMap[items[0]] = modelString;

			}

		}
		//cout <<items.size() << "\t" <<  items[0] << "\t => \t" << confMap[items[0]] ; //<< endl;
	}
	return confMap;
}

inline double dist(Point A, Point B) {
	return sqrt((B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y));
}

inline double area(Point A, Point B, Point C) {
	double side_a = dist(B, C);
	double side_b = dist(A, C);
	double side_c = dist(A, B);

	double s = 0.5*(side_a + side_b + side_c);
	return sqrt(s*(s-side_a)*(s-side_b)*(s-side_c));
}

double	lm_nfw_mass(double x) {
	double	logx_2=0,result=0;

	if (x < 0) {
		cout << "error running lm_nfw_mass" << endl; //fprintf(stderr,"lm_nfw_mass: invalid x: %g. Must be > 0\n",x);
		return result;
	}

	if (x==0) return 0;
	logx_2 = log(0.5 * x);
	if (x < 1.0) {
		result = lm_arccosh(1/x)/sqrt(1 - x*x);
	}
	else if (x == 1.0) {
		result =  1.0;
	}
	else {	 /* x > 1 */
		result =  acos(1.0/x)/sqrt(x*x -1.);
	}
	return 4.0*(logx_2 + result);
}


double lm_arctanh(double x) {
	if (x < -1 || x > 1.) {
		fprintf(stderr,"lm_arctanh: invalid x: %g. Must be 0 <= x <= 1\n",x);
		return 0;
	}
	return log(sqrt((1.+x)/(1.-x)));
}

double	lm_arccosh(double x) {
	if (x < 1.) {
		fprintf(stderr,"lm_arccosh: invalid x: %g. Must be >= 1.\n",x);
		return 0;
	}
	return log(sqrt(x*x-1.)+x);
}

vector<double> getTriWeight(Point A, Point B, Point C, Point P) {
	double areaA = area(P, B, C);
	double areaB = area(P, A, C);
	double areaC = area(P, A, B);
	double S = areaA + areaB + areaC;

	vector<double> w;
	w.push_back(areaA/S);
	w.push_back(areaB/S);
	w.push_back(areaC/S);
	return w;
}

double getPenalty(sp_mat* M, vec* r, vec* d, sp_mat* invC) {

	//  chi2 =  (M*r-d)^T(M*r-d)
	//cout << *M << endl;
	vec res = (*M)*(*r)-(*d);
	vec chi2 =  res.transpose()*(*invC)*res;

	return chi2(0,0);
}


void getLinearInterpolate(Point A, Point B,  Point C,  Point *P,  char direction) {

	double a=0, b=0; 

	if(direction=='x') {
		if(abs(A.x-B.x)<10e-8)  {
			P->x = 0.5*(A.x+B.x);
			P->y = C.y;
		}
		else {
			a = (A.y-B.y)/(A.x-B.x);
			b = A.y-a*A.x;
			P->x = (C.y-b)/a;
			P->y = C.y;
		}
	}
	if(direction=='y') {

			a = (A.y-B.y)/(A.x-B.x);
			b = A.y-a*A.x;
			P->x = C.x;
			P->y = a*C.x+b;

	}

}



vector<double> getPentWeigth(Point A, Point B, Point C, Point D, Point E) {
	vector<double> pentWeight;

	Point P(0, 0, 0);
	Point Q(0, 0, 0);
	Point M(0, 0, 0);
	Point N(0, 0, 0);
	getLinearInterpolate(A, B, C, &Q, 'x');
	//cout << Q.x << endl;
	//cout << Q.y << endl;
	//cout << Q.z << endl;
	getLinearInterpolate(D, E, C, &P, 'x');

	getLinearInterpolate(B, E, C, &M, 'y');
	getLinearInterpolate(A, D, C, &N, 'y');

	double dCQ_AB = dist(C, Q)*dist(A, B);

	double dCP_DE = dist(C, P)*dist(D, E);
	//double dAB = ;

	pentWeight.push_back(dist(Q, B)/dCQ_AB);    			// 	wAx
	pentWeight.push_back(dist(Q, A)/dCQ_AB); 				// 	wBx
	pentWeight.push_back(-(1/dist(C, P)+1/dist(C, Q)));  	//	wCx
	pentWeight.push_back(dist(P, E)/dCP_DE);   				//	wDx
	pentWeight.push_back(dist(P, D)/dCP_DE); 				//	wEx

	double dCN_AD = dist(C, N)*dist(A, D);
	double dCM_BE = dist(C, M)*dist(B, E);

	pentWeight.push_back(dist(N, D)/dCN_AD);   				//	wAy
	pentWeight.push_back(dist(M, E)/dCM_BE);   				//	wBy
	pentWeight.push_back(-(1/dist(C, N)+1/dist(C,M)));  	//	wCy
	pentWeight.push_back(dist(A, N)/dCN_AD);    			//	wDy
	pentWeight.push_back(dist(B, M)/dCM_BE);    			// 	wEy

	return pentWeight;

}



normVec getNormVector(Point A, Point B, Point C) {
	normVec norm; 
	norm.n0 = (C.y-B.y)*(B.z-A.z)-(C.z-B.z)*(B.y-A.y); 
	norm.n1 = (C.z-B.z)*(B.x-A.x)-(C.x-B.x)*(B.z-A.z); 
	norm.n2 = (C.x-B.x)*(B.y-A.y)-(C.y-B.y)*(B.x-A.x); 
	double s = sqrt(norm.n0*norm.n0 + norm.n1*norm.n1 + norm.n2*norm.n2);
	norm.n0 = norm.n0/s;
	norm.n1 = norm.n1/s;
	norm.n2 = norm.n2/s;

	return norm;
}

normVec meanNormVector(vector<normVec>  normList) {
	normVec meanNorm(0, 0, 0);
	for(int i=0; i<normList.size(); ++i) {
		meanNorm.n0 += normList[i].n0;
		meanNorm.n1 += normList[i].n1;
		meanNorm.n2 += normList[i].n2;
	}

	double s = meanNorm.n0*meanNorm.n0
			+ meanNorm.n1*meanNorm.n1
			+ meanNorm.n2*meanNorm.n2;

	if(s==0) {
		meanNorm.n0 = 0;
		meanNorm.n1 = 0;
		meanNorm.n2 = 0;
	}
	else {
		meanNorm.n0 = meanNorm.n0/s;
		meanNorm.n1 = meanNorm.n1/s;
		meanNorm.n2 = meanNorm.n2/s;
	}
	return meanNorm;

}

/***************************
Function:   	getAngularSizeDistance
Description:    Compuate comoving distance and angular distance
				from cosmological constant.
Arguments:		(1) Cosmological constants;
				(2) Redshift of object;
Returns:    	None
****************************/
void getAngularSizeDistance(Conf* conf, double z, double* comoveD, double* angularD) {
	// Cosmology calculator ala Ned Wright (www.astro.ucla.edu/~wright);
	// Asume a flat unverse.  Wv = 1-Omega_M;
	double WV = 1-conf->omega;
	double c = 299792.458;
	double DCMR = 0.0;      //  comoving radial distance in units of c/H0
	double WM = conf->omega;
	//double WR = 4.165E-5/(conf->hubble*conf->hubble) ;   // includes 3 massless neutrino species, T0 = 2.72528
	//double WK = 1-WM-WR-WV;
    double WR = 0.;
    double WK = 0.;
	double az = 1.0/(1+1.0*z);
	double a, adot;
	double n=1000;        // number of points in integrals
	for(int i=0; i<1000; ++i) {
		a = az+(1-az)*(i+0.5)/n;
		adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
		DCMR = DCMR + 1./(a*adot);
	}
	DCMR = (1.-az)*DCMR/n;
	double x = sqrt(abs(WK))*DCMR;
	double ratio;
	if(x > 0.1){
		if (WK > 0)
			ratio =  0.5*(exp(x)-exp(-x))/x;
		else
			ratio = sin(x)/x;
	}
	else {
		double y = x*x;
		if (WK < 0) y = -y;
		ratio = 1. + y/6. + y*y/120.;
	}
	double DCMT = ratio*DCMR;

	*angularD = (0.01*c/conf->hubble)*az*DCMT;  // in unit MPC
	*comoveD  = (0.01*c/conf->hubble)*DCMT;    // in unit MPC
}

/***************************
Function:   	getEisteinRadius
Description:    Compuate EisteinRadius based on the Mass of dark matter;
Arguments:		(1) Cosmological constants;
				(2) Mass;
Returns:    	Einstein Radius (arcsecond);
****************************/
double getEisteinRadius(Conf* conf, double Mtot) {
	// Mtot in unit M_sun;
	double G =  4.301e-9;      // in km^2 Mpc Msun^-1 s^-2
	double c = 299792.458;    // # velocity of light in km/sec
	double DMs, Ds, DMd, Dd, Dds;

	getAngularSizeDistance(conf, conf->srcZ, &DMs, &Ds);
	getAngularSizeDistance(conf, conf->lenZ, &DMd, &Dd);
	Dds = 1/(1+conf->srcZ)*(DMs-DMd);
	double R =  sqrt(4*G*Mtot*Dds/(c*c*Ds*Dd))*206265;

	return R;
}

/***************************
Function:   	sign
Description:    Sign function;
Arguments:
Returns:    	-1 for x<0;  1 for x>0;
****************************/
int sign(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

vec cV_to_eigenV(vector<double> s) {
	int n = s.size();
	vec v(n);
	for (int i=0; i<n; ++i) {
		v(i) = s[i];

	}
	return v;

}

vector<double> eigenV_to_cV(vec v) {
	int n = v.size();
	vector<double> s(n);
	for (int i=0; i<n; ++i) {
			s[i] = v(i);

	}
	return s;


}

sp_mat generatePSFoperator(string psfFileName, int naxis1, int naxis2) {
	/*Reference:   Treu-2004-Massive dark matter-Appendix B
	 * u 		: 	psfImage -> xList;
	 * v		:	psfImage -> yList;
	 * P(u,v)	:  	psfImage -> data;
	 */

	int N = naxis1;
	int M = naxis2;
	Image* psfImage = new Image(psfFileName);
	psfImage->updateFilterImage("none", 0);

	/*
	 *
	 * Normalize PSF data to be 1 (divided by the sum of FITS);
	 *
	 */
	psfImage->normalizeData();
	double psfLen = psfImage->naxis1*psfImage->naxis2;

	cout << "psfLen: " << psfLen << endl;
	int u, v;
	double p;
	sp_mat blurOperator(M*N, M*N);
	blurOperator.reserve(450*naxis1*naxis2);
	int g, h;
	int psfCenterX = (psfImage->naxis1-1)/2;
	int psfCenterY = (psfImage->naxis2-1)/2;
	cout << "centerX: " << psfCenterX << endl;
	cout << "centerY: " << psfCenterY << endl;
	int counter = 0;
	// create a psfMap;
	int Hx = (psfImage->naxis1-1)/2;
	int Hy = (psfImage->naxis2-1)/2;
	map<pair<int, int>,double> psfMap;
	for(int i=0; i<psfLen; ++i) {
		u = psfImage->xList[i]-Hx;
		v = psfImage->yList[i]-Hy;
		p = psfImage->data[i];
		//posMap[make_pair(imgX, imgY)] = i;
		psfMap[make_pair(u,v)] = p;

	}

	chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	int dx, dy;
	for(int x0=0; x0<naxis1; ++x0) {
		for(int y0=0; y0<naxis2; ++y0) {
			int limX_Low =(x0-Hx>0)?(x0-Hx):0;
			int limX_Up=(x0+Hx<naxis1)?(x0+Hx):naxis1;
			int limY_Low  =(y0-Hy>0)?(y0-Hy):0;
			int limY_Up =(y0+Hy<=naxis2)?(y0+Hy):naxis2;
			h = y0*naxis1 + x0;
			for(int iterX=limX_Low; iterX<limX_Up; ++iterX) {
				for(int iterY=limY_Low; iterY<limY_Up; ++iterY) {
					dx = x0 - iterX;
					dy = y0 - iterY;
					g = iterY*naxis1 + iterX;
					blurOperator.insert(h,g) = psfMap[make_pair(dx, dy)];
				}
			}
		}
	}

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout << "elapsed time _ 1: " << elapsed_seconds.count() << endl;

	cout << "Counter: " << counter <<endl;
	return blurOperator;
}


sp_mat generatePSFoperator(string psfFileName, Image* image) {
	/*Reference:   Treu-2004-Massive dark matter-Appendix B
	 * u 		: 	psfImage -> xList;
	 * v		:	psfImage -> yList;
	 * P(u,v)	:  	psfImage -> data;
	 */
	int naxis1 = image->naxis1;
	int naxis2 = image->naxis2;
	Image* psfImage = new Image(psfFileName);
	psfImage->updateFilterImage("none", 0);
	psfImage->normalizeData();
	int u, v, g;
	double p;
	sp_mat blurOperator(image->length,image->length);
	blurOperator.reserve(450*image->length);

	// create a psfMap;
	map<pair<int, int>,double> psfMap;
	int Hx = (psfImage->naxis1-1)/2;
	int Hy = (psfImage->naxis2-1)/2;
	for(int i=0; i<psfImage->length; ++i) {
		u = psfImage->xList[i]-Hx;
		v = psfImage->yList[i]-Hy;
		p = psfImage->data[i];
		psfMap[make_pair(u,v)] = p;
	}
	// create a posMap;
	map<pair<int, int>,int> posMap;
	for(int i=0; i<image->length; ++i)  {
		posMap[make_pair(image->xList[i], image->yList[i])] = i;
	}

	int dx, dy, x0, y0;
	for(int i=0; i<image->length; ++i) {
		x0 = image->xList[i];
		y0 = image->yList[i];

		int limX_Low =(x0-Hx>0)?(x0-Hx):0;
		int limX_Up=(x0+Hx<naxis1)?(x0+Hx):naxis1;
		int limY_Low  =(y0-Hy>0)?(y0-Hy):0;
		int limY_Up =(y0+Hy<=naxis2)?(y0+Hy):naxis2;

		for(int iterX=limX_Low; iterX<limX_Up; ++iterX) {
			for(int iterY=limY_Low; iterY<limY_Up; ++iterY) {
				dx = x0 - iterX;
				dy = y0 - iterY;
				g =  posMap[make_pair(iterX, iterY)];
				blurOperator.insert(i,g) = psfMap[make_pair(dx, dy)];
			}
		}

	}
	return blurOperator;
}




/***************************
Function:   	splitString
Description:    Split a string by whitespace;
Arguments:		(1) input string

Returns:    	string vector including all the split strings
****************************/
vector<string> splitString(string s) {
	vector<string> items;
	istringstream iss(s);
	//cout << iss << endl;
	copy(istream_iterator<string>(iss),
			istream_iterator<string>(),
			back_inserter(items));
	return items;
}





