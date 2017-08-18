//
//  testing.hpp
//  inference_test
//
//  Created by David Eyre on 06/01/2017.
//  Copyright © 2017 David Eyre. All rights reserved.
//

#ifndef testing_hpp
#define testing_hpp

#include <stdio.h>
#include <unordered_map>
#include <vector>
#include "likelihood.hpp"

using namespace std;

void runTest(vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes, Parm &parm,
             vector<int> &infSources, vector<int> &infSourceType,
             vector<vector<vector<int>>> &sporeI, vector<vector<vector<double>>> &sporeForce, vector<vector<double>> &sporeForceSummary,
             vector<vector<double>> &geneticDist, vector<vector<vector<int>>> &wardLogInf, vector<vector<int>> &wardLogNeverInf,
             vector<vector<vector<int>>> &inPtDays,
             vector<vector<int>> &ptLocation,
             vector<vector<int>> &wardI, int nInfPatients, int nNeverInfPatients, int nWards, int maxTime);

#endif /* testing_hpp */
