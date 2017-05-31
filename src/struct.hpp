//
//  struct.h
//  inference_test
//
//  Created by David Eyre on 27/01/2017.
//  Copyright © 2017 David Eyre. All rights reserved.
//

#ifndef struct_h
#define struct_h
#include <vector>
#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace std;

//structure for storing transmission events in
struct trans {
    int victim;
    int t;
    int srcType;
    int ward;
};

//structure for returning the details of a conditional source
struct Src {
    int srcIndex;
    int srcRoute;
    double srcP;
};

//structure for returning list of potential sources
struct SrcList {
    vector<int> sourceList; //vector of sources
    vector<int> sourceTypeList; //vector of source types
    vector<double> sourceProbabilities; //vector of probabilities for each source
};


//operator to allow comparison of transmission events
bool operator==(const trans& lhs, const trans& rhs);

//class to hold values of source types, e.g. SrcType::WARD == 1
class SrcType {
public:
    enum Src {
        BGROUND_HOSP = 0,
        WARD = 1,
        HOSP = 2,
        BGROUND_COMM = 3,
        START_POS = 4,
        SPORE = 5
    };
};

//structure to hold paramters
struct Parm {
    //"beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm", "sporeProb", "betaSpore", "recSize", "recMu"
    double betaBgroundHosp;
    double betaWard;
    double betaHosp;
    double sampleEpsilon;
    double directNe;
    double introNe;
    double mu;
    double probStartInf;
    double betaComm;
    double sporeProb; //for geometric distn
    double betaSpore;
    double recSize;
    double recMu;
};

//convert negative/positive number to 0 to 1 scale
double logistic(double x);

//convert 0 to 1 to any real number
double logit(double x);

//function to return spore prob
double getSporeP(vector<double> &parm);




#endif /* struct_h */
