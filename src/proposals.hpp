//
//  proposals.hpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

//This contains the Hasting's ratio for each of the possible proposals in the data augmentation

#ifndef proposals_hpp
#define proposals_hpp
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "tools.hpp"
#include "struct.hpp"

using namespace std;



//function to propose a new infection time (note will update recovery time at same iteration, and so not constrained by recovery)
int proposeInfectionTime(int i, int proposedPatient, int currentInfTime,
                         vector<vector<trans>> &onwardTransmission,
                         vector<int> &sampleTimes,
                         int maxTime,
                         double delta,
                         vector<vector<int>> &ptLocation);

//function for initiating infection times (at start or after move from being infected at t=0
int proposeInfectionTimeInitial(int proposedPatient, vector<int> &sampleTimes, Parm parm);

//function to propose a new recovery time
int proposeRecoveryTime(int proposedPatient, int currentRecTime,
                        vector<vector<trans>> &onwardTransmission, vector<int> &sampleTimes, int maxTime, int minTime, double delta,
                        vector<vector<int>> &ptLocation, vector<int> &infTimes);


//get list of sources and their probabilities for a given patient and infection time, and parameter set
SrcList getSourceProb(int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                      vector<vector<vector<int>>> &wardLogInf, vector<int> infSourceType,
                      vector<vector<vector<SporeEvent>>> &sporePatientI,
                      vector<vector<int>> &ptLocation,
                      vector<vector<double>> &geneticDist,
                      int nWards,
                      int nInfPatients, vector<vector<int>> &hospitalWards,  Parm &parm);



//function to determine proposed source of infection conditionally
Src proposeConditionalSource(int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                             vector<vector<vector<int>>> &wardLogInf, vector<int> infSourceType,
                             vector<vector<vector<SporeEvent>>> &sporePatientI, vector<vector<int>> &ptLocation,
                             vector<vector<double>> &geneticDist,
                             int nWards, int nInfPatients, vector<vector<int>> &hospitalWards, Parm &parm);

#endif /* proposals_hpp */
