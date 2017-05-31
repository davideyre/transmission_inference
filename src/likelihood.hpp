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

double llTransAvoid(vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &infSourceType, vector<int> &infSources,
                    vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
                    vector<vector<vector<int>>> &wardLog,
                    vector<vector<vector<int>>> &inPtDays,
                    vector<vector<int>> &ptLocation,
                    vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, vector<double> &parm);

//log likelihood contribution from transmission model - p(I | parm)
double llTrans(vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &infSourceType, vector<int> &infSources,
               vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
               vector<vector<vector<int>>> &wardLog,
               vector<vector<vector<int>>> &inPtDays,
               vector<vector<int>> &ptLocation,
               vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, vector<double> &parm);


//log likelihood contribution from sampling times - p(S | I, parm)
double llSample(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<double> parm);
double llSampleMP(vector<int> &infTimes, vector<int> &sampleTimes, vector<double> parm);

//log likelihood contribution for recover, p(R | I, parm)
double llRecover(vector<int> &infectedPatients, vector<int> &sampleTimes, vector<int> &recTimes, vector<double> parm);


//log likelihood contribution from the genetic distance matrix - p(G | I, S, parm)
double llGenetic(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, vector<double> parm);

//log likelihood over all pairs calling the llGeneticSingle function (slower x4, but used for testing only)
double llGeneticAlt(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, vector<double> parm);


//genetic log likelihood for single patient pair
double llGeneticSingle(vector<int> &infectedPatients, vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, vector<double> parm);

//prior
double getPrior(vector<double>parm);


//target distribution, i.e. non-normalised posterior
double targetDist (vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes,
                   vector<int> &infSources, vector<int> &infSourceType,
                   vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
                   vector<vector<vector<int>>> &wardLog,
                   vector<vector<vector<int>>> &inPtDays,
                   vector<vector<int>> &ptLocation,
                   vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap, vector<double> &parm);

#endif /* likelihood_hpp */
