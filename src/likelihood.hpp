//
//  likelihood.hpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#ifndef likelihood_hpp
#define likelihood_hpp

#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <unordered_map>
#include <vector>
#include "tools.hpp"
#include "struct.hpp"

using namespace std;


//log likelihood contribution from transmission model - p(I | parm)
double llTrans(vector<int> &infTimes, vector<int> &infSourceType, vector<int> &infSources,
               vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
               vector<vector<vector<int>>> &wardLogInf, vector<vector<int>> &wardLogNeverInf,
               vector<vector<vector<int>>> &inPtDays,
               vector<vector<int>> &ptLocation,
               vector<vector<int>> &wardI, int nInfPatients, int nNeverInfPatients, int nWards, int maxTime, Parm &parm);


//log likelihood contribution from sampling times - p(S | I, parm)
double llSample(int nInfPatients, vector<int> &infTimes, vector<int> &sampleTimes, Parm &parm);

//log likelihood contribution for recover, p(R | I, parm)
double llRecover(int nInfPatients, vector<int> &sampleTimes, vector<int> &recTimes, Parm &parm);


//log likelihood contribution from the genetic distance matrix - p(G | I, S, parm)
double llGenetic(vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm);


//genetic log likelihood for single patient pair
double llGeneticSingle(vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm);

//prior
double getPrior(Parm &parm);


//target distribution, i.e. non-normalised posterior
double targetDist (vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes,
                   vector<int> &infSources, vector<int> &infSourceType,
                   vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
                   vector<vector<vector<int>>> &wardLogInf, vector<vector<int>> &wardLogNeverInf,
                   vector<vector<vector<int>>> &inPtDays,
                   vector<vector<int>> &ptLocation,
                   vector<vector<int>> &wardI, int nInfPatients, int nNeverInfPatients, int nWards, int maxTime, vector<vector<double>> &geneticDist, Parm &parm);

#endif /* likelihood_hpp */
