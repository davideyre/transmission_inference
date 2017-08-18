//
//  import.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#include "import.hpp"


// *** functions to import epi data *** //

//function to import details for infected patients (allows import of known infection and recovery times, source and source type for testing purposes)
void importPatientLog(string filePath, unordered_map<string,int> &ptLookup, unordered_map<int,string> &ptLookupRev, vector<int> &infTimes, vector<int> &infSources,
                      vector<int> &infSourceType, vector<int> &sampleTimes, vector<int> &recoverTimes, int &nInfPatients, int &nNeverInfPatients) {

    //patient log has to contain all never infected patients
    //first set up ptLookup, save original patient ids and then re-number infected cases consecutively from zero
    //also save sampleTimes
    io::CSVReader<2> in(filePath);
    in.read_header(io::ignore_extra_column, "patient_id", "t_sample");
    string tmp_ptId, tmp_sampleTimes;
    int i = 0;
    nNeverInfPatients = 0;
    while(in.read_row(tmp_ptId, tmp_sampleTimes)) {
        if(tmp_sampleTimes!="NA") {
            ptLookup.insert( { tmp_ptId, i });
            ptLookupRev.insert( { i, tmp_ptId });            
            i++;
            sampleTimes.push_back(stoi(tmp_sampleTimes)-1);
        }
        else {
            nNeverInfPatients ++;
        }
    }
    nInfPatients = i;
    
    //resize infection source, types, recovery times to be correct size
    infTimes.resize(nInfPatients);
    infSources.resize(nInfPatients);
    infSourceType.resize(nInfPatients);
    recoverTimes.resize(nInfPatients);
    
    
    //now generate list of infection source, types, recovery times
    io::CSVReader<6> inMore(filePath);
    inMore.read_header(io::ignore_extra_column, "patient_id", "t_inf", "source", "source_type", "t_sample", "t_recover");
    string tmp_infTimes, tmp_infSources, tmp_infSourceType, tmp_recoverTimes;
    while(inMore.read_row(tmp_ptId, tmp_infTimes, tmp_infSources, tmp_infSourceType, tmp_sampleTimes, tmp_recoverTimes)){
        //save data - converting all times to start from 0, and patient identifiers to numeric
        int ptId;
        
        if(tmp_sampleTimes!="NA") { //for patients with an infection -
            //retrieve numeric patient id
            ptId = ptLookup.at(tmp_ptId);
            //populate vectors
            infTimes[ptId] = stoi(tmp_infTimes)-1;
            recoverTimes[ptId] = stoi(tmp_recoverTimes)-1;
            infSourceType[ptId] = stoi(tmp_infSourceType);
            //check if background source, i.e. -1
            if(tmp_infSources=="-1") {
                infSources[ptId] = -1;
            }
            else {
                int srcPtId = ptLookup.at(tmp_infSources);
                infSources[ptId] = srcPtId;
            }
        }
    }
}




//function to import ward admissions
// create two outputs -
//    wardLogInf -  // wardLogInf[time][ward] = {patients...} for infected patients
//    wardLogNeverInf - // wardLogNeverInf[time][ward] = int for count of never infected patients

void importWardLog(string filePath, unordered_map<string,int> &wardLookup, unordered_map<string,int> &ptLookup,
                   vector<vector<vector<int>>> &wardLogInf,
                   vector<vector<int>> &wardLogNeverInf,
                   int &maxTime, int &nWards, vector<int> &sampleTimes) {
    
    //temporary variables to hold text from CSV
    int tmp_admit, tmp_discharge;
    string tmp_ward, tmp_pt;
    
    //create a lookup for all the wards, and obtain the last discharge date at the same time
    io::CSVReader<2> in(filePath);
    in.read_header(io::ignore_extra_column, "ward", "t_discharge");
    int lastDischarge = 0;
    int i = 0;
    while(in.read_row(tmp_ward, tmp_discharge)){
        auto it = wardLookup.find(tmp_ward);
        if (it==wardLookup.end()) {
            wardLookup.insert( { tmp_ward, i });
            i++;
        }
        tmp_discharge -= 1; //convert tmp_discharge to number from zero
        if (tmp_discharge > lastDischarge) {
            lastDischarge = tmp_discharge;
        }
    }
    
    //get max time
    int maxTimeWard = lastDischarge;
    int maxTimeSample = *max_element(begin(sampleTimes), end(sampleTimes));
    maxTime = max({maxTimeWard, maxTimeSample});
    
    //get number of wards
    nWards = i;
    
    //set up sizes of wardLogs - wardLogInf,
    wardLogInf.resize(maxTime+1);
    wardLogNeverInf.resize(maxTime+1);
    for (int t = 0; t<=maxTime; t++) {
        wardLogInf[t].resize(nWards);
        wardLogNeverInf[t].resize(nWards);
        //set wardLogInf to be zero at all times
        for(int ward=0; ward<nWards; ward++) {
            wardLogNeverInf[t][ward] = 0;
        }
    }
    
    io::CSVReader<4> inMore(filePath);
    inMore.read_header(io::ignore_extra_column, "patient_id", "ward", "t_admit", "t_discharge");
    while(inMore.read_row(tmp_pt, tmp_ward, tmp_admit, tmp_discharge)){
        //convert admission and discharge dates
        int t_admit = tmp_admit - 1;
        int t_discharge = tmp_discharge -1;
        
        //convert ward to numeric form
        int ward = wardLookup.at(tmp_ward);
        
        //check if infected patient
        auto itPt = ptLookup.find(tmp_pt);
        if (itPt!=ptLookup.end()) {
            //infected patient
            int pt = ptLookup.at(tmp_pt);
            for (int t=t_admit; t<=t_discharge; t++) {
                wardLogInf[t][ward].push_back(pt);
            }
        }
        else {
            //uninfected patient
            for (int t=t_admit; t<=t_discharge; t++) {
                wardLogNeverInf[t][ward] += 1;
            }
        }
    }
}



//function to import genetic data
void importGeneticData(string filePathGenetic, vector<vector<double>> &geneticDist, unordered_map<string,int> &ptLookup, int nInfPatients) {
    string line;
    ifstream file(filePathGenetic);
    int lineNumber = 0;
    int fieldNumber = 0;
    string field;
    vector<int> ptIdHeader;
    
    //reszie genetic distance
    geneticDist.resize(nInfPatients);
    for (int pt =0; pt<nInfPatients; pt++) {
        geneticDist[pt].resize(nInfPatients);
    }
    
    //populate genetic distance
    if (file.is_open()) {
        while (getline(file, line)) {
            lineNumber ++;
            istringstream iss(line);
            if (lineNumber==1) {
                //read in header - which contains patient identifiers
                while(std::getline(iss, field, ' ')) {
                    string strRemove = "patient_";
                    string::size_type i = field.find(strRemove);
                    if (i != std::string::npos) {
                        field.erase(i, strRemove.length());
                    }
                    int pt = ptLookup.at(field);
                    ptIdHeader.push_back(pt);
                }
            }
            else {
                fieldNumber = 0;
                int pt1, pt2;
                while(std::getline(iss, field, ' ')) {
                    string startsWith = "patient_";
                    if(field.find(startsWith)==0) {
                        string::size_type i = field.find(startsWith);
                        if (i != std::string::npos) {
                            field.erase(i, startsWith.length());
                        }
                        pt1 = ptLookup.at(field);
                    }
                    else {
                        pt2 = ptIdHeader[fieldNumber];
                        geneticDist[pt1][pt2] = stod(field);
                        
                        fieldNumber ++;
                        
                    }
                }
            }
        }
        file.close();
    }
}
