/*
 * mc.cpp
 *
 *  Created on: Apr, 2016
 *      Author: En-Hsin Peng
 */

#include "mc.h"
#include <iostream>

MC::MC(unsigned seed) {
    rng_engine.seed(seed);
    rng = std::bind(std::uniform_real_distribution<double>(0.,1.), std::ref(rng_engine));
    makeCgauss();
}

double MC::random() {
    return rng();
}

double MC::gauss(double sigma, double mu, double x) {
    double dx = x- mu;
    return exp(-dx*dx/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}

void MC::makeCgauss() {
    cgArr.resize(100,0.);
    cgArr[0] = gauss(1., 0.,  -5.)*0.1;
    for (size_t i=1; i<100; ++i) {
        cgArr[i] = cgArr[i-1] + gauss(1., 0.,  -5. + i*0.1)*0.1;
    }
}

double MC::cgauss() {
    double xv = 0.;
    double r = rng();
    for (size_t i=0; i<100; ++i) {
        if (cgArr[i]>r) {
            xv = (-5. + i*0.1) - 0.1*rng();
            break;
        }
    }
    return xv;
}

double MC::stepPar(MultModelParam &param, double cfac, size_t &iter) {
    double minSig = 1e-6;
    double eps=0.02;
    cfac *= (1+eps);
    if (cfac > 1e10) {
        cfac = 1.0;
        for(int j=0; j<param.nLens; ++j) {
            for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
                param.mixAllModels[0][j].paraList[k] = 0.;
                param.mixAllModels[1][j].paraList[k] = 0.;
                param.mixAllModels[2][j].paraList[k] = 0.;
            }
        }
    }
    // mixAllModels[0]: weight * par^2
    // mixAllModels[1]: weight
    // mixAllModels[2]: weight * par
    // mixAllModels[3]: new parameter
    // mixAllModels[4]: original parameter
    // mixAllModels[5]: best parameter
    // mixAllModels[6]: lower bound
    // mixAllModels[7]: upper bound
    // mixAllModels[8]: step
    for(int j=0; j<param.nLens; ++j) {
        for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
            if (param.mixAllModels[6][j].paraList[k] < param.mixAllModels[7][j].paraList[k]) {
                double stepSig = param.mixAllModels[8][j].paraList[k];
                double par0 = param.mixAllModels[4][j].paraList[k];

                param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
                param.mixAllModels[1][j].paraList[k] += cfac;
                param.mixAllModels[2][j].paraList[k] += cfac*par0;
                if (iter > 3) {
                    double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                            - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
                    stepSig = 3.*2.38*sig/sqrt(param.mixAllModels[0][j].paraList.size());
                    if (std::isnan(stepSig)) {
                        iter=0;
                        param.mixAllModels[0][j].paraList[k] = 0.;
                        param.mixAllModels[1][j].paraList[k] = 0.;
                        param.mixAllModels[2][j].paraList[k] = 0.;
                        stepSig = param.mixAllModels[8][j].paraList[k];
                    }
                }
                if (stepSig < minSig) stepSig = minSig;
                param.mixAllModels[3][j].paraList[k] += cgauss()*stepSig;
                if (param.mixAllModels[3][j].paraList[k] < param.mixAllModels[6][j].paraList[k])
                    param.mixAllModels[3][j].paraList[k] = param.mixAllModels[6][j].paraList[k];
                if (param.mixAllModels[3][j].paraList[k] > param.mixAllModels[7][j].paraList[k])
                    param.mixAllModels[3][j].paraList[k] = param.mixAllModels[7][j].paraList[k];
            }
        }
    }
    return cfac;
}
