//
//  tools.hpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#ifndef tools_hpp
#define tools_hpp
#include <vector>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <algorithm> //for removing first element from a vector
#include <unordered_map>
#include "likelihood.hpp"
#include "proposals.hpp"
#include "struct.hpp"

using namespace std;




//function to obtain the number of infectious individuals on a given ward at a specific timepoint
int getI(int t, int ward, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLogInf);

//function to return a 2d vector of the number of infectious individuals on each ward at all time points
vector<vector<int>> getWardI(int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLogInf);

//function to get sporeI - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
void getSporeI(vector<vector<vector<int>>> &sporeI, int nInfPatients, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes,
               vector<vector<int>> &ptLocation);

//udpate sporeI for a single patient and return proposed copy
void updateSporeI(vector<vector<vector<int>>> &sporeI, int updatePt, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes,
                  vector<vector<int>> &ptLocation);

void getSporeForceSummary(vector<vector<double>> &sporeForceSummary, vector<vector<vector<int>>> &sporeI,
                          int maxTime, int nWards, int nInfPatients, vector<int> &infTimes, Parm parm);


//function to get vector of days an intpatient - inPtDays[patient][ward] = {times...} (whereas wardLog[time][ward] = {patients...})
vector<vector<vector<int>>> getInPtDays(int nInfPatients, int maxTime, int nWards,
                                        vector<vector<vector<int>>> &wardLogInf);

//function to store the location of patients - ptLocation[patient][time] = wardId
vector<vector<int>> getPtLocation(int nInfPatients, int maxTime, int nWards, vector<vector<vector<int>>> &inPtDays);

// Incomplete gamma function, as defined in Maple - required in llGenetic
double igamma (double a, double z);
double logIgamma (double a, double z); //log version

//factorial function
double factorial (double x);

//function to return log factorial
double logFactorial(double x);

//remove the first element from a vector that matches
void removeFirst(int n, vector<int> &vect);

//create a normalised probability vector from vector of log likelihoods
vector<double> normaliseLL(vector<double> sourceLikelihood);

//choose item based on vector of probabilities, return the index of the chosen item
int sampleProbVector(vector<double> sourceProbability);



#endif /* tools_hpp */
