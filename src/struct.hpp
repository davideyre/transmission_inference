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
    //parameters for model
    //"beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm", "sporeProb", "betaSpore", "recSize", "recMu"
    double betaBgroundHosp;
    double betaWard;
    double betaHosp;
    double sampleEpsilon;
    double directNe;
    double introNe;
    double mu;
    double probStartInfLogit;
    double betaComm;
    double sporeProbLogit; //for geometric distn
    double betaSpore;
    double recSize;
    double recMu;
    
    //likelihood values stored with parameters in chain
    double currentLL, currentLLTrans, currentLLGenetic, currentLLSample, currentLLRecovery;

    //allow retrival by numerical index
    double &operator[]( size_t idx ) {
        switch( idx ) {
            case 0 : return betaBgroundHosp;
            case 1 : return betaWard;
            case 2 : return betaHosp;
            case 3 : return sampleEpsilon;
            case 4 : return directNe;
            case 5 : return introNe;
            case 6 : return mu;
            case 7 : return probStartInfLogit;
            case 8 : return betaComm;
            case 9 : return sporeProbLogit;
            case 10 : return betaSpore;
            case 11 : return recSize;
            case 12 : return recMu;
                
            case 13 : return currentLL;
            case 14 : return currentLLTrans;
            case 15 : return currentLLGenetic;
            case 16 : return currentLLSample;
            case 17 : return currentLLRecovery;
                
            default: throw std::runtime_error( "Parameter structure: bad index\n" );
        }
    }
    
};


//convert negative/positive number to 0 to 1 scale
double logistic(double x);

//convert 0 to 1 to any real number
double logit(double x);

//function to return spore prob
double getSporeP(Parm &parm);

//function to return start infected prob
double getStartInfP(Parm &parm);


#endif /* struct_h */
