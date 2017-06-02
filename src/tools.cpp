
//
//  tools.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//


#include "tools.hpp"






//function to obtain the number of infectious individuals on a given ward at a specific timepoint
int getI(int t, int ward, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLog) {
    int nInf = 0;
    for (int pt : wardLog[t][ward]) {
        if (infTimes[pt] > -1 & infTimes[pt] < t & recTimes[pt] > t) {
            nInf ++;
        }
    }
    return nInf;
}

//function to return a 2d vector of the number of infectious individuals on each ward at all time points
vector<vector<int>> getWardI(int nPatients, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes, vector<vector<vector<int>>> &wardLog) {
    vector<vector<int>> wardI;
    
    //set up sizes of 2d vector
    wardI.resize(maxTime+1);
    for (int i = 0; i<=maxTime; i++) {
        wardI[i].resize(nWards);
    }
    
    //populate
    for (int t=0; t<maxTime; t++) {
        for (int ward=0; ward< nWards; ward++) {
            wardI[t][ward] = getI(t, ward, infTimes, recTimes, wardLog);
        }
    }
    return wardI;
}


//function to get sporeI - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
void getSporeI(vector<vector<vector<int>>> &sporeI, vector<int> &infectedPatients, int nPatients, int maxTime, int nWards, vector<int> &infTimes, vector<int> &recTimes,
                                      vector<vector<int>> &ptLocation) {
    
    //reset sporeI for all infected patients
    for(int pt : infectedPatients) {
        for(int t=0; t<=maxTime; t++) {
            for(int ward=0; ward<nWards; ward++) {
                sporeI[t][ward][pt] = 0;
            }
        }
    }


    //iterate through infected patients
    for(int pt : infectedPatients) {
        //set spores at recovery if inpatient - ptLocation[patient][time] = wardId
        if(recTimes[pt]<=maxTime) {
            int recoveryWard = ptLocation[pt][recTimes[pt]];
            if(recoveryWard>-1) {
                int sporeStart = recTimes[pt];
                int sporeEnd = maxTime;
                int sporeAge = 1;
                for (int t=sporeStart; t<=sporeEnd; t++) {
                    sporeI[t][recoveryWard][pt] = sporeAge;
                    sporeAge++;
                }
            }
        }
        
        //set spores for each discharge while still infectious, i.e. from time step after infected to time step before recovery
        int currentWard = ptLocation[pt][infTimes[pt]+1];
        for(int t = infTimes[pt]+1; t<= min({recTimes[pt], maxTime}); t++) { //infectious from day after infected, until day before recover (need to go to next day to find discharge on day of recovery)
            if(ptLocation[pt][t] != currentWard & currentWard!=-1) {
                
                //discharge has occured
                int sporeStart = t;
                int sporeEnd = maxTime;
                int sporeAge = 1;
                for (int tt=sporeStart; tt<=sporeEnd; tt++) {
                    sporeI[tt][currentWard][pt] = sporeAge;
                    sporeAge ++;
                }
            }
            currentWard = ptLocation[pt][t];
        }
    }// end of patients loop
}

//udpate sporeI for a single patient, updating proposed copy passed as sporeI
void updateSporeI(vector<vector<vector<int>>> &sporeI, int updatePt, int maxTime, int nPatients, int nWards, vector<int> &infTimes, vector<int> &recTimes,
          vector<vector<int>> &ptLocation) {

    //change sporeI sporeI[t][ward][pt] = spore_duration (numbering from 1...) for patient being updated
    
    //set current values for sporeI for this patient to zero
    for(int t=0; t<=maxTime; t++) {
        for(int ward=0; ward<nWards; ward++) {
            sporeI[t][ward][updatePt] = 0;
        }
    }
    
    //set spores at recovery if inpatient - ptLocation[patient][time] = wardId
    if(recTimes[updatePt]<=maxTime) {
        int recoveryWard = ptLocation[updatePt][recTimes[updatePt]];
        if(recoveryWard>-1) {
            int sporeStart = recTimes[updatePt];
            int sporeEnd = maxTime;
            int sporeAge = 1;
            for (int t=sporeStart; t<sporeEnd; t++) {
                sporeI[t][recoveryWard][updatePt] = sporeAge;
                sporeAge++;
            }
        }
    }
    
    //set spores for each discharge while still infectious, i.e. from time step after infected to time step before recovery
    int currentWard = ptLocation[updatePt][infTimes[updatePt]+1];
    for(int t = infTimes[updatePt]+1; t<= min({recTimes[updatePt], maxTime}); t++) { //infectious from day after infected, until day before recover (need to go to next day to find discharge on day of recovery)
        if(ptLocation[updatePt][t] != currentWard & currentWard!=-1) {
            
            //discharge has occured
            int sporeStart = t;
            int sporeEnd = maxTime;
            int sporeAge = 1;
            for (int tt=sporeStart; tt<=sporeEnd; tt++) {
                sporeI[tt][currentWard][updatePt] = sporeAge;
                sporeAge ++;
            }
        }
        currentWard = ptLocation[updatePt][t];
    }
}



void getSporeForceSummary(vector<vector<double>> &sporeForceSummary, vector<int> &infectedPatients, vector<vector<vector<int>>> &sporeI, int maxTime, int nWards, int nPatients, vector<int> &infTimes, Parm parm) {
    //pre-calculate the sum of spore force of infection - [t][ward]

    //reset current values of spore force summary to be zero
    for (int t = 0; t<=maxTime; t++) { //can only add spore after infected, hence start from there
        for (int ward = 0; ward<nWards; ward++) {
            sporeForceSummary[t][ward] = 0;
        }
    }
    
    for (int sporePt: infectedPatients) {
        for (int t = max({0, infTimes[sporePt]}); t<=maxTime; t++) { //can only add spore after infected, hence start from there or t=0 if later, as spore only set after t=0
            for (int ward = 0; ward<nWards; ward++) {
                int sporeDuration = sporeI[t][ward][sporePt];
                if(sporeDuration>0) {
                    sporeForceSummary[t][ward] += pow((1-getSporeP(parm)), sporeDuration);
                }
            }
        }
    }
}







//function to get vector of days an intpatient - inPtDays[patient][ward] = {times...} (whereas wardLog[time][ward] = {patients...})
vector<vector<vector<int>>> getInPtDays(int nPatients, int maxTime, int nWards, vector<vector<vector<int>>> &wardLog) {
    vector<vector<vector<int>>> inPtDays;
    
    //set up 3d vector
    inPtDays.resize(nPatients);
    for (int i=0; i<nPatients; i++) {
        inPtDays[i].resize(nWards);
    }
    
    //loop through wardLog
    for (int t=0; t<=maxTime; t++) {
        for (int ward=0; ward< nWards; ward++) {
            for(int pt : wardLog[t][ward]) {
                inPtDays[pt][ward].push_back(t);
            }
        }
    }
    return inPtDays;
    
}

//function to store the location of patients - ptLocation[patient][time] = wardId
vector<vector<int>> getPtLocation(int nPatients, int maxTime, int nWards, vector<vector<vector<int>>> &inPtDays) {
    vector<vector<int>> ptLocation;
    
    //set up 3d vector
    ptLocation.resize(nPatients);
    for (int pt=0; pt<nPatients; pt++) {
        ptLocation[pt].resize(maxTime+1);
        for (int t=0; t<=maxTime; t++) {
            ptLocation[pt][t] = -1; //set default location as in community
        }
    }
    
    //loop through inPtDays[patient][ward] = {times}
    for (int pt=0; pt<nPatients; pt++) {
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

//factorial function
double factorial (double x) {
    return tgamma(x+1);
}


//function to return log factorial
double logFactorial(double x) {
    return lgamma(x+1);
}


//function to get vector of infected patients
vector<int> getInfectedPatients(vector<int> &sampleTimes, int nPatients) {
    vector<int> infectedPatients;
    for(int patient=0; patient<nPatients; patient++) {
        if(sampleTimes[patient]>=0) {
            infectedPatients.push_back(patient);
        }
    }
    return infectedPatients;
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






