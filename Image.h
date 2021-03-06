/*
 * Image.h
 *
 *  Created on: Oct 8, 2015
 *      Author: cheng109
 */

#ifndef IMAGE_H_
#define IMAGE_H_


#include<string>
#include<vector>
//#include <armadillo>
//#include "Model.h"
using namespace std;
#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::VectorXd vec;

class Conf;

class Image {
friend class Model;
friend class Const;
	string fileName;
public:
	int naxis;
	long naxis1;
	long naxis2;
	long npixels;

	int bitpix;
	double res;



	vector<double> data;
	vector<double> varList;
	/*** ImageList After Filter  ***/
	vector<double> dataList;
	vector<int> xList;
	vector<int> yList;
	vector<int> iList;

	vector<int> type;
	long length;

	sp_mat invC;

public:
	Image();

	Image(vector<double> xpos, vector<double> ypos,vector<double> *briList, long naxis1, long naxis2, int bitpix);
	Image(vector<int> xpos, vector<int> ypos,vector<double> *briList, long naxis1, long naxis2, int bitpix);

	Image(string imgFileName);
	void getConstants(long *filterPixelNum, long* naxis1, long* naxis2, double *res, int* bit);
	void printImageInfo(int x1, int y1, int x2, int y2);
	void updateFilterImage(string regionFileName,int flag) ;
	void updateGridPointType();
	void writeFilterImage(string imgFileName);
	void writeToFile(string imgFileName);
	void updateVarList(double threshold, double backVar);
	void updateVarList(string varFileName, string regionFileName);
	void normalizeData();
	vec getMatrixD();
	sp_mat getVarMatrix();
	sp_mat getPSFMatrix(string psfFileName, long dim);
	void erasePixel(int index);
	int sign(double x) ;
	virtual ~Image();
};


#endif /* IMAGE_H_ */



