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
            if(isPosInt(tmp_infTimes)) infTimes[ptId] = stoi(tmp_infTimes)-1;
            if(isPosInt(tmp_recoverTimes)) recoverTimes[ptId] = stoi(tmp_recoverTimes)-1;
            if(isPosInt(tmp_infSourceType)) infSourceType[ptId] = stoi(tmp_infSourceType);
            //check if background source, i.e. -1
            if(tmp_infSources=="-1") {
                infSources[ptId] = -1;
            }
            else {
                if(!tmp_infSources.empty() & tmp_infSources!="NA") {
                    int srcPtId = ptLookup.at(tmp_infSources);
                    infSources[ptId] = srcPtId;
                }

            }
        }
    }
}




//function to import ward admissions
// create two outputs -
//    wardLogInf -  // wardLogInf[time][ward] = {patients...} for infected patients
//    wardLogNeverInf - // wardLogNeverInf[time][ward] = int for count of never infected patients

void importWardLog(string filePath, unordered_map<string,int> &hospitalLookup, unordered_map<string,int> &wardLookup, unordered_map<string,int> &ptLookup,
                   vector<vector<vector<int>>> &wardLogInf,
                   vector<vector<int>> &wardLogNeverInf,
                   int &maxTime, int &nWards, vector<int> &sampleTimes, vector<vector<int>> &hospitalWards, unordered_map<int,int> &ward2Hospital, vector<vector<int>> &hospitalWardList) {
    
    //temporary variables to hold text from CSV
    int tmp_admit, tmp_discharge;
    string tmp_ward, tmp_pt, tmp_hospital;
    
    //create a lookup for all the wards, and obtain the last discharge date at the same time
    io::CSVReader<3> in(filePath);
    in.read_header(io::ignore_extra_column, "ward", "hospital", "t_discharge");
    int lastDischarge = 0;
    int i = 0;
    int h = 0;
    while(in.read_row(tmp_ward, tmp_hospital, tmp_discharge)){
        //check if have already saved hospital and if not save in lookup
        auto itHosp = hospitalLookup.find(tmp_hospital);
        if (itHosp==hospitalLookup.end()) {
            hospitalLookup.insert( {tmp_hospital, h});
            h++;
            hospitalWardList.push_back({});
        }
        
        string tmp_wardHosp = tmp_hospital+": "+tmp_ward;
        auto it = wardLookup.find(tmp_wardHosp);
        if (it==wardLookup.end()) {
            wardLookup.insert( { tmp_wardHosp, i }); //save name of ward
            int hChk = hospitalLookup.at(tmp_hospital); //save which hospital ward is in
            hospitalWardList[hChk].push_back(i);
            i++;
        }
        tmp_discharge -= 1; //convert tmp_discharge to number from zero
        if (tmp_discharge > lastDischarge) {
            lastDischarge = tmp_discharge;
        }
    }
    
    //get number of wards
    nWards = i;
    int nHospitals = h;
    
    //generate list of wards in the same hospital as each ward
    hospitalWards.resize(nWards);
    for (int hosp=0; hosp<nHospitals; hosp++) {
        for (int ptWard : hospitalWardList[hosp]) {
            for (int ward: hospitalWardList[hosp]) {
                if (ward!=ptWard) {
                    hospitalWards[ptWard].push_back(ward);
                }
            }
        }
    }
    
    //generate a list of hospitals for each ward
    for (int hospital=0; hospital<hospitalWardList.size(); hospital++) {
        for (int ward : hospitalWardList[hospital]) {
            ward2Hospital.insert( {ward, hospital} );
        }
    }
    
    
    //get max time
    int maxTimeWard = lastDischarge;
    int maxTimeSample = *max_element(begin(sampleTimes), end(sampleTimes));
    maxTime = max({maxTimeWard, maxTimeSample});
    

    
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
    
    io::CSVReader<5> inMore(filePath);
    inMore.read_header(io::ignore_extra_column, "patient_id", "ward", "hospital", "t_admit", "t_discharge");
    while(inMore.read_row(tmp_pt, tmp_ward, tmp_hospital, tmp_admit, tmp_discharge)){
        //convert admission and discharge dates
        int t_admit = tmp_admit - 1;
        int t_discharge = tmp_discharge -1;
        
        //convert ward to numeric form
        string tmp_wardHosp = tmp_hospital+": "+tmp_ward;
        int ward = wardLookup.at(tmp_wardHosp);
        
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
            istringstream iss(line);
            if (lineNumber==0) {
                //read in header - which contains patient identifiers
                while(std::getline(iss, field, ' ')) {
                    int pt = ptLookup.at(field);
                    ptIdHeader.push_back(pt);
                }
            }
            else {
                //read in non-header rows
                fieldNumber = 0;
                int pt1, pt2;
                //for each field in row
                while(std::getline(iss, field, ' ')) {
                    if(fieldNumber==0) {
                        //for first field, get patient id from row label
                        pt1 = ptLookup.at(field);
                    }
                    else {
                        //for subsequent fields save SNPs
                        pt2 = ptIdHeader[fieldNumber-1];
                        geneticDist[pt1][pt2] = stod(field);
                        
                    }
                    fieldNumber ++;
                } //end fields in row loop
            } //end non-header row if...
            lineNumber ++;
        }
        file.close();
    }
}
