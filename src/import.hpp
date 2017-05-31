//
//  import.hpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#ifndef import_hpp
#define import_hpp
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include "csv.hpp"
#include "proposals.hpp"

using namespace std;

//function to import ward log
void importWardLog(string filePath, vector<vector<vector<int>>> &wardLog, int &maxTime, int &nWards, int nPatients,
                   vector<int> &sampleTimes);

//function to import patient infection log
void importPatientLog(string filePath, vector<int> &infTimes, vector<int> &infSources, 
                      vector<int> &infSourceType, vector<int> &sampleTimes, vector<int> &recoverTimes, int &nPatients);

//function to import genetic data
void importGeneticData(string filePathGenetic, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap);

#endif
