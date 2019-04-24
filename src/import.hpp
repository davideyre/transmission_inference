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
#include "tools.hpp"

using namespace std;


//function to import patient infection log
void importPatientLog(string filePath, unordered_map<string,int> &ptLookup, unordered_map<int,string> &ptLookupRev, vector<int> &infTimes, vector<int> &infSources,
                      vector<int> &infSourceType, vector<int> &sampleTimes, vector<int> &recoverTimes, int &nInfPatients, int &nNeverInfPatients);

//function to import ward log
void importWardLog(string filePath, unordered_map<string,int> &hospitalLookup, unordered_map<string,int> &wardLookup,
                   unordered_map<int,string> &wardLookupRev,
                   unordered_map<string,int> &ptLookup, vector<vector<vector<int>>> &wardLogInf,
                   vector<vector<int>> &wardLogNeverInf, int &maxTime, int &nWards, vector<int> &sampleTimes, vector<vector<int>> &hospitalWards, vector<int> &ward2Hospital, vector<vector<int>> &hospitalWardList);

//function to import genetic data
void importGeneticData(string filePathGenetic, vector<vector<double>> &geneticDist, unordered_map<string,int> &ptLookup, int nInfPatients);

#endif
