
//
//  tools.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//


#include "tools.hpp"






//function to obtain the number of infectious individuals on a given ward at a specific timepoint
int getI(int t, int ward, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLogInf) {
    int nInf = 0;
    for (int pt : wardLogInf[t][ward]) {
        if (infTimes[pt] < t & recTimes[pt] > t) {
            nInf ++;
        }
    }
    return nInf;
}

//function to return a 2d vector of the number of infectious individuals on each ward at all time points
vector<vector<int>> getWardI(int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLogInf) {
    
    vector<vector<int>> wardI;
    
    //set up sizes of 2d vector
    wardI.resize(maxTime+1);
    for (int i = 0; i<=maxTime; i++) {
        wardI[i].resize(nWards);
    }
    
    //populate
    for (int t=0; t<maxTime; t++) {
        for (int ward=0; ward< nWards; ward++) {
            wardI[t][ward] = getI(t, ward, infTimes, recTimes, wardLogInf);
        }
    }
    return wardI;
}


//function to get sporePatientI - sporePatientI[ward][pt] = vector of time intervals spore first present in (allow for multiple ward discharges while infectious)
void getSporePatientI(vector<vector<vector<int>>> &sporePatientI, int nInfPatients, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes, vector<vector<int>> &ptLocation) {
    
    //reset sporePatientI for all infected patients to empty vector
    for(int ward=0; ward<nWards; ward++) {
        for(int pt=0; pt<nInfPatients; pt++) {
            sporePatientI[ward][pt] = {};
        }
    }


    //iterate through infected patients
    for(int pt=0; pt<nInfPatients; pt++) {
        // A) set spores at recovery if inpatient - ptLocation[patient][time] = wardId
        if(recTimes[pt]<=maxTime) {
            int recoveryWard = ptLocation[pt][recTimes[pt]];
            if(recoveryWard>-1) {
                sporePatientI[recoveryWard][pt].push_back(recTimes[pt]);
            }
        }
        // B) set spores for each discharge while still infectious, i.e. from time step after infected to time step before recovery
        int t0 = max({0, infTimes[pt]+1});
        int currentWard = ptLocation[pt][t0];
        for(int t = t0; t<= min({recTimes[pt], maxTime}); t++) { //infectious from day after infected, until day before recover (need to go to next day to find discharge on day of recovery)
            if(ptLocation[pt][t] != currentWard & currentWard!=-1) {
                //discharge has occured at time, t-1, set spores to start at time t on previous location
                sporePatientI[currentWard][pt].push_back(t);
            }
            currentWard = ptLocation[pt][t];
        }
    }// end of patients loop
}

//udpate sporePatientI for a single patient - sporePatientI[ward][pt] = vector of time intervals spore first present in (allow for multiple ward discharges while infectious)
void updateSporePatientI(vector<vector<vector<int>>> &sporePatientI, int updatePt, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes, vector<vector<int>> &ptLocation) {

    //set current values for sporePatientI for this patient to empty vector
    for(int ward=0; ward<nWards; ward++) {
        sporePatientI[ward][updatePt] = {};
    }
    
    //set spores at recovery if inpatient - ptLocation[patient][time] = wardId
    if(recTimes[updatePt]<=maxTime) {
        int recoveryWard = ptLocation[updatePt][recTimes[updatePt]];
        if(recoveryWard>-1) {
            sporePatientI[recoveryWard][updatePt].push_back(recTimes[updatePt]);
        }
    }
    
    //set spores for each discharge while still infectious, i.e. from time step after infected to time step before recovery
    int t0 = max({0, infTimes[updatePt]+1});
    int currentWard = ptLocation[updatePt][t0];
    for(int t = t0; t<= min({recTimes[updatePt], maxTime}); t++) { //infectious from day after infected, until day before recover (need to go to next day to find discharge on day of recovery)
        if(ptLocation[updatePt][t] != currentWard & currentWard!=-1) {
            //discharge has occured at time, t-1, set spores to start at time t, set spores on previous ward
            sporePatientI[currentWard][updatePt].push_back(t);
        }
        currentWard = ptLocation[updatePt][t];
    }
}


//function to pre-calculate the sum of spore force of infection - sporeForceSummary[t][ward]
void getSporeForceSummary(vector<vector<double>> &sporeForceSummary, vector<vector<vector<int>>> &sporePatientI, int maxTime, int minTime, set<int> wardsToUpdate, int nInfPatients, vector<int> &infTimes, vector<vector<int>> &ptLocation, Parm parm) {

    //pre-calculate all spore probabilities
    double probSpore = 1-getSporeP(parm);
    vector<double> probSporeDay;
    probSporeDay.resize(maxTime+1);
    for (int t=0; t<=(maxTime-minTime); t++) {
        probSporeDay[t] = pow(probSpore, t);
    }
    
    
    //iterate through wards
    for (int ward : wardsToUpdate) {
        //reset sporeForceSummary for all times between min and max time
        for (int t=minTime; t<=maxTime; t++) {
            sporeForceSummary[t][ward] = 0;
        }
        
        //for each ward get any patient who sets spores
        for (int sporePt=0; sporePt< nInfPatients; sporePt++) {
            if(!sporePatientI[ward][sporePt].empty()) {
                // for each time spores are set
                for (int sporeTime : sporePatientI[ward][sporePt]) {
                    for (int t=sporeTime; t<=maxTime; t++) {
                        int sporeDuration = t - sporeTime + 1;
                        sporeForceSummary[t][ward] += probSporeDay[sporeDuration];
                    }
                }
            }
        }
    }
}




//function to get vector of days an intpatient - inPtDays[patient][ward] = {times...} (whereas wardLogInf[time][ward] = {patients...})
vector<vector<vector<int>>> getInPtDays(int nInfPatients, int maxTime, int nWards, vector<vector<vector<int>>> &wardLogInf) {
    vector<vector<vector<int>>> inPtDays;
    
    //set up 3d vector
    inPtDays.resize(nInfPatients);
    for (int i=0; i<nInfPatients; i++) {
        inPtDays[i].resize(nWards);
    }
    
    //loop through wardLog
    for (int t=0; t<=maxTime; t++) {
        for (int ward=0; ward< nWards; ward++) {
            for(int pt : wardLogInf[t][ward]) {
                inPtDays[pt][ward].push_back(t);
            }
        }
    }
    return inPtDays;
    
}

//function to store the location of patients for infected patients - ptLocation[patient][time] = wardId
vector<vector<int>> getPtLocation(int nInfPatients, int maxTime, int nWards, vector<vector<vector<int>>> &inPtDays) {
    vector<vector<int>> ptLocation;
    
    //set up 3d vector
    ptLocation.resize(nInfPatients);
    for (int pt=0; pt<nInfPatients; pt++) {
        ptLocation[pt].resize(maxTime+1);
        for (int t=0; t<=maxTime; t++) {
            ptLocation[pt][t] = -1; //set default location as in community
        }
    }
    
    //loop through inPtDays[patient][ward] = {times}
    for (int pt=0; pt<nInfPatients; pt++) {
        for (int ward=0; ward< nWards; ward++) {
            for (int t: inPtDays[pt][ward]) {
                ptLocation[pt][t] = ward;
            }
        }
    }

    return ptLocation;
    
}





// Incomplete gamma function, as defined in Maple - required in llGenetic
double igamma (double a, double z) {
    return pgamma(z,a,1,0,0)*tgamma(a);
}

double logIgamma (double a, double z) {
    return pgamma(z,a,1,0,1) + lgamma(a);
}

//factorial function
double factorial (double x) {
    return tgamma(x+1);
}


//function to return log factorial
double logFactorial(double x) {
    return lgamma(x+1);
}


//remove the first element from a vector that matches
void removeFirst(int n, vector<int> &vect){
    vector<int>::iterator found = find(vect.begin(), vect.end(), n) ;
    if (found!=vect.end())
        vect.erase(found);
}




//create a normalised probability vector from vector of log likelihoods
vector<double> normaliseLL(vector<double> sourceLikelihood) {
    
    // see - http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
    
    int n = (int)sourceLikelihood.size(); //number of likelihoods
    double limit = log(pow(10, -16)) - log(n); //limit of accuracy
    
    double maxLikelihood = *max_element(begin(sourceLikelihood), end(sourceLikelihood));
    double totalLikelihood = 0;
    
    for (double likelihood : sourceLikelihood) {
        double nl = likelihood-maxLikelihood;
        if (nl > limit) {
            
            totalLikelihood += exp(nl);
        }
    }
    
    vector<double> sourceProbability;
    for (double likelihood : sourceLikelihood) {
        double nl = likelihood-maxLikelihood;
        if (nl <= limit) {
            sourceProbability.push_back(0);
        }
        else {
            sourceProbability.push_back(exp(nl)/totalLikelihood);
        }
    }
    
    return sourceProbability;
}


//choose item based on vector of probabilities, return the index of the chosen item
int sampleProbVector(vector<double> sourceProbability) {
    vector<double> sourceCumulProb;
    double cumul = 0;
    //get cumulative probability
    for (double prob : sourceProbability) {
        cumul += prob;
        sourceCumulProb.push_back(cumul);
    }
    
    double test = runif(0, 1);
    int i = 0;
    for (double cumul : sourceCumulProb) {
        if( test < cumul) {
            //this is the probability to use
            break;
        }
        else {
            i ++; //incremenet counter
        }
    }
    
    return i;
}


//check if a string is a postive integer
bool isPosInt(string& s)
{
    return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}



