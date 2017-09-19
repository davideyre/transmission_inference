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
#include <stdexcept> //required for runtime_error
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

//structure for holding spore events
struct SporeEvent {
    int start; //date spores start being present
    int end; //date spores end, i.e. 1 time step prior to readmission
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

//class to hold paramters
class Parm {
    
public:
    //parameters for model
    double betaBgroundHosp;
    double betaWard;
    double betaHosp;
    double sampleSize, sampleMu; //neg binomial parameters for time between infection and sampling
    double directNe;
    double introNe;
    double mu;
    double probStartInfLogit;
    double betaComm;
    double sporeProbLogit; //for geometric distn
    double recSize;
    double recMu;
    double sporeMultiplier;
    
    //likelihood values stored with parameters in chain
    double currentLL, currentLLTrans, currentLLGenetic, currentLLSample, currentLLRecovery;

    //allow retrival by numerical index
    double &operator[]( size_t idx ) {
        switch( idx ) {
            case 0 : return betaBgroundHosp;
            case 1 : return betaWard;
            case 2 : return betaHosp;
            case 3 : return sampleSize;
            case 4 : return sampleMu;
            case 5 : return directNe;
            case 6 : return introNe;
            case 7 : return mu;
            case 8 : return betaComm;
            case 9 : return sporeProbLogit;
            case 10 : return probStartInfLogit;
            case 11 : return recSize;
            case 12 : return recMu;
            case 13 : return sporeMultiplier;
                
            case 14 : return currentLL;
            case 15 : return currentLLTrans;
            case 16 : return currentLLGenetic;
            case 17 : return currentLLSample;
            case 18 : return currentLLRecovery;
                
            default: throw runtime_error( "Parameter structure: bad index\n" );
        }
    }
    
private:
    double logistic(double x) {
        return 1/(1+exp(-x));
    }
    
public:
    void displayLog() {
        printf("Current Values\n beta0: %0.7f\t beta1: %0.7f\t beta2: %0.7f\n sampleSize: %0.4f\t sampleMu: %0.4f\n directNe: %0.4f\t introNe: %0.4f\tmu: %0.4f\n pStartInfLogit: %0.6f\t pStartInf: %0.6f\n betaComm: %0.7f\n parm.sporeProbLogit: %0.4f\tsporeProb: %0.4f\tsporeMultiplier: %0.4f\n recSize: %0.4f\t recMu %0.4f\n\nCurrentLL:  %0.1f\n",
               betaBgroundHosp, betaWard, betaHosp, sampleSize, sampleMu, directNe, introNe, mu, probStartInfLogit, logistic(probStartInfLogit), betaComm, sporeProbLogit, logistic(sporeProbLogit), logistic(sporeMultiplier), recSize, recMu, currentLL);
        
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

//function to report spore multiplier
double getSporeMultiplier(Parm &parm);

#endif /* struct_h */
