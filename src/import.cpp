//
//  import.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#include "import.hpp"


//function to import epi data
//function to import patient infection log
void importPatientLog(string filePath, vector<int> &infTimes, vector<int> &infSources,
                      vector<int> &infSourceType, vector<int> &sampleTimes, vector<int> &recoverTimes, int &nPatients) {
    io::CSVReader<5> in(filePath);
    in.read_header(io::ignore_extra_column, "t_inf", "source", "source_type", "t_sample", "t_recover");
    string tmp_infTimes; string tmp_infSources; string tmp_infSourceType; string tmp_sampleTimes; string tmp_recoverTimes;
    
    while(in.read_row(tmp_infTimes, tmp_infSources, tmp_infSourceType, tmp_sampleTimes, tmp_recoverTimes)){
        //save data - converting to be indexed from 0
        if(tmp_infTimes=="NA") { infTimes.push_back(-1); }
        else { infTimes.push_back(stoi(tmp_infTimes)-1); }
        
        if(tmp_infSources=="NA" | tmp_infSources=="-1") { infSources.push_back(-1); }
        else { infSources.push_back(stoi(tmp_infSources)-1); }
        
        if(tmp_infSourceType=="NA") {infSourceType.push_back(-1);}
        else { infSourceType.push_back(stoi(tmp_infSourceType)); }
        
        if(tmp_sampleTimes=="NA") { sampleTimes.push_back(-1); }
        else { sampleTimes.push_back(stoi(tmp_sampleTimes)-1); }
        
        if(tmp_recoverTimes=="NA") { recoverTimes.push_back(-1); }
        else { recoverTimes.push_back(stoi(tmp_recoverTimes)-1); }
        
    }
    nPatients = (int)infTimes.size();
}




//function to import ward admissions
void importWardLog(string filePath, vector<vector<vector<int>>> &wardLog, int &maxTime, int &nWards, int nPatients,
                   vector<int> &sampleTimes) {
    io::CSVReader<4> in(filePath);
    in.read_header(io::ignore_extra_column, "patient_id", "ward", "t_admit", "t_discharge");
    int tmp_pt, tmp_ward, tmp_admit, tmp_discharge;
    vector<int> ptList, wardList, admitList, dischargeList;
    
    while(in.read_row(tmp_pt, tmp_ward, tmp_admit, tmp_discharge)){
        //save data - converting to be indexed from 0
        ptList.push_back(tmp_pt-1);
        wardList.push_back(tmp_ward-1);
        admitList.push_back(tmp_admit-1);
        dischargeList.push_back(tmp_discharge-1);
    }
    
    
    int maxTimeWard = *max_element(begin(dischargeList), end(dischargeList));
    int maxTimeSample = *max_element(begin(sampleTimes), end(sampleTimes));
    maxTime = max({maxTimeWard, maxTimeSample});
    
    nWards = *max_element(begin(wardList), end(wardList))+1;
    
    //set up sizes of wardLog
    wardLog.resize(maxTime+1);
    for (int i = 0; i<=maxTime; i++) {
        wardLog[i].resize(nWards);
    }
    
    //populate wardLog -  // wardLog[time][ward] = {patients...}
    for (int logIndex=0; logIndex<ptList.size(); logIndex++) {
        int pt = ptList[logIndex];
        int ward = wardList[logIndex];
        int t_admit = admitList[logIndex];
        int t_discharge = dischargeList[logIndex];
        for (int t=t_admit; t<=t_discharge; t++) {
            
            wardLog[t][ward].push_back(pt);
        }
    }
 
}

//function to import genetic data
void importGeneticData(string filePathGenetic, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap) {
    string line;
    ifstream file(filePathGenetic);
    int lineNumber = 0;
    int fieldNumber = 0;
    string field;
    vector<int> header;
    
    if (file.is_open()) {
        while (getline(file, line)) {
            lineNumber ++;
            istringstream iss(line);
            if (lineNumber==1) {
                while(std::getline(iss, field, ' ')) {
                    string strRemove = "patient_";
                    string::size_type i = field.find(strRemove);
                    if (i != std::string::npos) {
                        field.erase(i, strRemove.length());
                    }
                    header.push_back(stoi(field));
                    geneticMap[stoi(field)-1] = fieldNumber;
                    //printf("patient %d, map %d\n",stoi(field)-1, fieldNumber);
                    fieldNumber ++;
                }
                //resize vector to store genetic distances
                geneticDist.resize(header.size());
            }
            else {
                while(std::getline(iss, field, ' ')) {
                    string startsWith = "patient";
                    if(field.find(startsWith)!=0) {
                        geneticDist[lineNumber-2].push_back(stod(field));
                    }
                }
            }
        }
        file.close();
    }
}
