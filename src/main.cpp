//
//  main.cpp
//  inference_test
//
//  Created by David Eyre on 02/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

#include "import.hpp"
#include "tools.hpp"
#include "proposals.hpp"
#include "likelihood.hpp"
#include "export.hpp"
#include "testing.hpp"
#include "struct.hpp"

#include "optionparser.hpp"


#include <sys/time.h>

using namespace std;


//function to run MCMC
void doMCMC(vector<Parm> &chain, vector<vector<int>> &chainInfTimes, vector<vector<int>> &chainRecTimes, vector<vector<int>> &chainInfSources,
              vector<vector<int>> &chainInfSourceTypes,
              int steps, Parm startParm, Parm startSigma,
              vector<int> &sampleTimes,
              vector<vector<vector<int>>> &wardLogInf,
              vector<vector<int>> &wardLogNeverInf,
              int nInfPatients, int nNeverInfPatients, int nWards, int maxTime,
              vector<vector<double>> &geneticDist,
              vector<int> &infTimes, vector<int> &infSources, vector<int> &infSourceTypes, vector<int> &recoveryTimes,
              unordered_map<string,int> &ptLookup, unordered_map<string,int> &wardLookup, vector<vector<int>> &hospitalWards,
              vector<int> &ward2Hospital, vector<vector<int>> &hospitalWardList) {
    
    //note the following are supplied for testing purposes only -
    // vector<int> &infTimes, vector<int> &infSources, vector<int> &infSourceTypes, vector<int> &recoveryTimes
    
    
    int updateN = round(nInfPatients*0.1); // select a subset of patients to update times for at each iteration, e.g. 10%
    bool debugPt = false;
    

    //reshape the ward log - inPtDays[patient][ward] = {times...} (whereas wardLog[time][ward] = {patients...})
    vector<vector<vector<int>>> inPtDays = getInPtDays(nInfPatients, maxTime, nWards, wardLogInf);
    //reshape again: ptLocation[patient][time] = wardId
    vector<vector<int>> ptLocation = getPtLocation(nInfPatients, maxTime, nWards, inPtDays);
    
    //get list of wards ever visited by a patient - wardEver[pt] = {ward1, ward2, ...}
    vector<vector<int>> wardEver = getWardEver(inPtDays, nInfPatients, nWards);

    //onward transmission log, indexed by source patient (contains time of onward transmission and victim)
    vector<vector<trans>> onwardTransmission;
    onwardTransmission.resize(nInfPatients);
    trans defaultTrans;
    defaultTrans.t = maxTime;
    defaultTrans.victim = -1;
    defaultTrans.srcType = -1;
    defaultTrans.ward = -1;
    for(int patient=0; patient<nInfPatients; patient++) {
        onwardTransmission[patient].push_back(defaultTrans);
    }
    
    //set starting values for infection times and infection sources
    vector<int> startInfTimes, startInfSources, startInfSourceTypes, startRecTimes;
    
    //loop over all patients to set starting values for infection times, sources, and recovery times
    for(int patient=0; patient<nInfPatients; patient++) {

        //infection time and source
        int infTime = proposeInfectionTimeInitial(patient, sampleTimes, startParm);
        startInfTimes.push_back(infTime);
        startInfSources.push_back(-1); //background source - assume all from background initially
        
        //recovery time
        double recSize = startParm.recSize;
        double recMu = startParm.recMu;
        double recProb = recSize / (recSize + recMu);
        int t_rec = sampleTimes[patient] + rnbinom(recSize, recProb); //assume initial recovery times (approx as assume for this infected at t=0)
        startRecTimes.push_back(t_rec);
        
        //infection source type
        if(infTime<0) {
            //start infected
            startInfSourceTypes.push_back(SrcType::START_POS); //start infected source type
        }
        else {
            //check location at sampling time
            if(ptLocation[patient][infTime]>-1) {
                startInfSourceTypes.push_back(SrcType::BGROUND_HOSP); //inpatient
            }
            else {
                startInfSourceTypes.push_back(SrcType::BGROUND_COMM); //community
            }
        }

    } //end patient loop
    

    //initialise chain with parameter starting values
    chain[0] = startParm;
    chainInfTimes[0] = startInfTimes;
    chainRecTimes[0] = startRecTimes;
    chainInfSources[0] = startInfSources;
    chainInfSourceTypes[0] = startInfSourceTypes;
    
    //set up acceptance counters - "beta0", "beta1", "beta2", "sampleSize", "sampleMu", "directNe", "introNe", "mu", "startInf", "betaComm", "sporeProb", "recSize", "recMu", "sporeMultiplier", "augmentation moves - infection", "aug - recovery", "disruption moves"
    vector<int> nAccepted = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    //initial sd for MH updates
    Parm sigma = startSigma;
    
    //variables used during MCMC
    Parm currentParm; currentParm = chain[0];
    Parm proposedParm;
    double proposedValue, probAccept;
    
    //current and proposed values
    vector<int> currentInfTimes, currentRecTimes, currentInfSources, currentInfSourceTypes, proposedInfTimes, proposedRecTimes, proposedInfSources, proposedInfSourceTypes;
    vector<vector<vector<SporeEvent>>> currentSporePatientI, proposedSporePatientI;
    vector<vector<double>> currentSporeForceSummary, proposedSporeForceSummary;
    vector<vector<int>> currentWardI, proposedWardI; //set up wardI log - //create array of the number of infected individuals on each ward at each time point
    
    
    //reserve memory for sporePatientI - current and proposed
    currentSporePatientI.resize(nWards);
    proposedSporePatientI.resize(nWards);
    for (int ward = 0; ward<nWards; ward++) {
        currentSporePatientI[ward].resize(nInfPatients);
        proposedSporePatientI[ward].resize(nInfPatients);
    }
    
    //reserve memory for sporeForceSummary
    //current
    currentSporeForceSummary.resize(maxTime+1);
    for (int t = 0; t<=maxTime; t++) {
        currentSporeForceSummary[t].resize(nWards);
        for (int ward = 0; ward<nWards; ward++) {
            currentSporeForceSummary[t][ward] = 0.00;
        }
    }
    //proposed
    proposedSporeForceSummary.resize(maxTime+1);
    for (int t = 0; t<=maxTime; t++) {
        proposedSporeForceSummary[t].resize(nWards);
        for (int ward = 0; ward<nWards; ward++) {
            proposedSporeForceSummary[t][ward] = 0.00;
        }
    }
    
    currentInfTimes = startInfTimes;
    currentInfSources = startInfSources;
    currentInfSourceTypes = startInfSourceTypes;
    currentRecTimes = startRecTimes;
    
    //TEMP OVER_RIDE FOR TESTING - use for parameter checking only as doesn't set up onward transmission log
    //currentInfTimes = infTimes;
    //currentInfSources = infSources;
    //currentInfSourceTypes = infSourceTypes;
    //currentRecTimes = recoveryTimes;
    
    
    currentWardI = getWardI(maxTime, nWards, currentInfTimes, currentRecTimes, wardLogInf); //get initial values
    proposedWardI = currentWardI;
    getSporePatientI(currentSporePatientI, nInfPatients, maxTime, nWards, currentInfTimes, currentRecTimes, ptLocation);
    int minTime = 0;
    set<int> allWards;
    for (int ward=0; ward<nWards; ward++) {
        allWards.insert(ward);
    }
    
    
    getSporeForceSummary(currentSporeForceSummary, currentSporePatientI, maxTime, minTime, allWards, nInfPatients, currentInfTimes, ptLocation, currentParm);
    //also set proposedSporePatientI to start as current SporePatientI as only update this at the accept / reject step
    proposedSporePatientI = currentSporePatientI;
    
   
    //calculate initial targetDist value
    double currentLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes,
                                  currentInfSources, currentInfSourceTypes,
                                  currentSporePatientI, currentSporeForceSummary,
                                  wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
    printf("Starting target distribution value: %0.4f\n\n", currentLL);
    double proposedLL;
    
    
    //get time at start of MCMC
    struct timeval previousSystemTime;
    gettimeofday(&previousSystemTime, NULL);

    struct timeval previousSystemTimeIter;
    gettimeofday(&previousSystemTimeIter, NULL);
    
    //generate random walk
    for(int i=1; i<steps; i++) {
        
        /*
        //timer for every iteration for testing
        printf("Starting iteration %d\n", i);
        struct timeval currentSystemTimeIter;
        gettimeofday(&currentSystemTimeIter, NULL);
        if(i>1) {            
            double deltaIter = ((currentSystemTimeIter.tv_sec  - previousSystemTimeIter.tv_sec) * 1000000u +
                            currentSystemTimeIter.tv_usec - previousSystemTimeIter.tv_usec) / 1.e6;
            printf("iteration processing time: %0.1f seconds\n\n",deltaIter);
            
        }
        previousSystemTimeIter = currentSystemTimeIter;
        */
        
        //PARAMETER AND TUNING REPORTING
        //every 100 steps
        if(i % 100 == 0) {
            printf("Starting iteration %d\n", i);
            //tune proposal distributions - update sigma based on nAccepted
            printf("Accepted Moves\n beta0: %d\tbeta1: %d\tbeta2: %d\tsampleSize: %d\tsampleMu: %d\tdirectNe: %d\tintroNe: %d\tmu: %d\t betaComm: %d\nsporeProb: %d\tsporeMultiplier: %d\tpStartInf: %d\trecSize %d\trecMu %d\n",
                   nAccepted[0], nAccepted[1], nAccepted[2], nAccepted[3], nAccepted[4], nAccepted[5], nAccepted[6], nAccepted[7], nAccepted[8], nAccepted[9], nAccepted[13], nAccepted[10], nAccepted[11], nAccepted[12]);
            for (int j=0; j<14; j++) {
                if(nAccepted[j]<20) {
                    sigma[j] = sigma[j] * 0.8;
                } else if (nAccepted[j]>40) {
                    sigma[j] = sigma[j] * 1.25;
                }
                nAccepted[j] = 0;
            }
            
            //print augmentation moves
            printf("Augmentation (infection) moves accepted %d\n", nAccepted[14]);
            nAccepted[14] = 0;
            printf("Augmentation (recovery) moves accepted %d\n", nAccepted[15]);
            nAccepted[15] = 0;
            
            //print disruption moves - cumulative total
            printf("Disruption moves accepted %d\n", nAccepted[16]);
            nAccepted[16] = 0;
            
            
            //display current parameter values
            chain[i-1].displayLog();
            
            //display current LL for different components of likelihood
            double currentLLSample = llSample(nInfPatients, currentInfTimes, sampleTimes, currentParm);
            printf("LL sample, current: %0.3f\n", currentLLSample);
            
            double currentLLGenetic = llGenetic(currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, nInfPatients, currentParm);
            printf("LL genetic, current: %0.3f\n", currentLLGenetic);
            
            double currentLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporePatientI, currentSporeForceSummary, wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
            printf("LL trans, current: %0.3f\n", currentLLTrans);
            
            double currentLLRecovery = llRecover(nInfPatients, sampleTimes, currentRecTimes, currentParm);
            printf("LL recovery, current: %0.3f\n\n", currentLLRecovery);
            
            //estimate time remaining
            struct timeval currentSystemTime;
            gettimeofday(&currentSystemTime, NULL);
            double delta = ((currentSystemTime.tv_sec  - previousSystemTime.tv_sec) * 1000000u +
                     currentSystemTime.tv_usec - previousSystemTime.tv_usec) / 1.e6;
            
            double estTotalRunTime = delta*steps/100/60;
            double estRemainingRunTime = delta*(steps-i)/100/60;
            printf("estimated MCMC total processing time: %0.1f minutes, remaining time: %0.1f minutes\n\n",estTotalRunTime,estRemainingRunTime);
            previousSystemTime = currentSystemTime;
            
            
        }

        
        //PARAMETER UPDATES
        for (int parmIndex=0; parmIndex<14; parmIndex++) {
                
            //iterate through parmameters - see Parm class for list in struct.hpp
            //1. proposal distribution, i.e. proposed change in value of beta
            proposedValue = currentParm[parmIndex] + rnorm(0 , sigma[parmIndex]);
            
            /*
            if(parmIndex==13) {
                proposedValue=0.4; //temporary over-ride to fix at simualted value
            }
            */
            
            proposedParm = currentParm; proposedParm[parmIndex] = proposedValue;
            if(proposedValue <= 0 & parmIndex !=9 & parmIndex !=10 & parmIndex !=13) {
                //if invalid proposed value skip this step - probabilities on logit scale so can be negative
                chain[i][parmIndex] = chain[i-1][parmIndex];
            }
            else {
                
                if(parmIndex==9) {
                    getSporeForceSummary(proposedSporeForceSummary, currentSporePatientI, maxTime, minTime, allWards, nInfPatients, currentInfTimes, ptLocation, proposedParm);
                    proposedLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                            currentSporePatientI, proposedSporeForceSummary,
                                            wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, proposedParm);
                }
                else {
                    //2. Compute the probability of accepting the proposed jump.
                    proposedLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                            currentSporePatientI, currentSporeForceSummary,
                                            wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, proposedParm);
                }
                
                
                
                probAccept = min( log(1), proposedLL -  currentLL);
                //3. Accept proposed jump with probability pAccept
                if ( log(runif(0,1)) < probAccept ) {
                    
                    chain[i][parmIndex] = proposedValue; // accept the proposed jump
                    nAccepted[parmIndex] += 1; // increment the accepted counter
                    //printf("Move accepted for parameter %d at iteration %d (current: %0.4f, proposed %0.4f\n", parmIndex, i, currentParm[parmIndex], proposedValue);
                    currentParm = proposedParm; //update for next parameter
                    currentLL = proposedLL;
                    
                    //if updating spore probabity update the sporeforce and summary
                    if(parmIndex ==9) {
                        currentSporeForceSummary = proposedSporeForceSummary;
                    }
                    
                }
                else {
                    // reject the proposed jump, stay at current position
                    chain[i][parmIndex] = currentParm[parmIndex];
                    
                    //reset proposedSporeForceSummary
                    if(parmIndex ==9) {
                        proposedSporeForceSummary = currentSporeForceSummary;
                    }
                    
                    //printf("Move rejected for parameter %d at iteration %d (current: %0.4f, proposed %0.4f\n", parmIndex, i, currentParm[parmIndex], proposedValue);
                }
            }
            
        } //end 14 parameter updates
        
        
        
        //over-ride - fix paramters
        
        //TRUE values
        
//        chain[i].betaBgroundHosp = 0.002;
//        chain[i].betaWard = 0.005;
//        chain[i].betaHosp = 0.0000001;
//        chain[i].sampleSize = 5;
//        chain[i].sampleMu = 10;
//        chain[i].directNe = 1;
//        chain[i].introNe = 10000;
//        chain[i].mu = 0.005475702;
//        chain[i].betaComm = 0.000001;
//        chain[i].sporeProbLogit = logit(0.2);
//        chain[i].probStartInfLogit = logit(0.00001);
//        chain[i].recSize = 3;
//        chain[i].recMu = 30;
        
        
        /*
        //over-ride for data augmentation
        
        currentInfTimes = infTimes;
        currentRecTimes = recoveryTimes;
        currentInfSources = infSources;
        currentInfSourceTypes = infSourceTypes;
        currentWardI = getWardI(maxTime, nWards, currentInfTimes, currentRecTimes, wardLogInf);
        getSporePatientI(currentSporePatientI, nInfPatients, maxTime, nWards, currentInfTimes, currentRecTimes, ptLocation);
        getSporeForceSummary(currentSporeForceSummary, currentSporePatientI, maxTime, minTime, allWards, nInfPatients, currentInfTimes, ptLocation, currentParm);
     
        //generic over-ride
    
        currentParm = chain[i];
        currentLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                               currentSporePatientI, currentSporeForceSummary,
                               wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
        
        
        currentLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes,
                               currentInfSources, currentInfSourceTypes,
                               currentSporePatientI, currentSporeForceSummary,
                               wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
        
         getSporeForceSummary(currentSporeForceSummary, currentSporePatientI, maxTime, minTime, allWards, nInfPatients, currentInfTimes, ptLocation, currentParm);
        // end over-ride
*/
        
        //DATA AUGMENTATION
        
        // get list of patients to update at this iteration
        vector<int> patientUpdateList(nInfPatients);
        iota(begin(patientUpdateList), end(patientUpdateList), 0); //create a list of all infected patients
        random_shuffle(begin(patientUpdateList), end(patientUpdateList)); //shuffle the patient list - UPDATE TO TAKE SEED
        patientUpdateList.resize(updateN); //keep only the first 10%

        //loop through patients to update infection times and sources
        for (int proposedPatient : patientUpdateList) {
            
            //debugging code
            debugPt = false;
            if (proposedPatient==-1) {
                printf("i=%d, patient=%d, current LL: %0.3f\n", i, proposedPatient, currentLL);
                double currentLLCheck = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                                   currentSporePatientI, currentSporeForceSummary, wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
                printf("LL check: %0.3f\n", currentLLCheck);
                debugPt = true;
            }
            
            
            //set proposed values to current values
            proposedInfTimes = currentInfTimes;
            proposedInfSources = currentInfSources;
            proposedInfSourceTypes = currentInfSourceTypes;
            
            //proposedSporePatientI = currentSporePatientI; //costly to repeat here - handle this at each accept / reject step
            
            vector<vector<trans>> proposedOnwardTransmission = onwardTransmission;
            
            double hastingsRatio = 0; //store log hastings ratio
            
            //propose infection infection time
            double sdInfTime = 5;
            if(i % 10 == 0) {
                sdInfTime = 25; // for 1 in 10 updates propose a bigger SD to improve chances of moving between adjacent admissions
            }
             //sd for normal distribution for infection time updates
            int proposedInfTime = proposeInfectionTime(proposedPatient, currentInfTimes[proposedPatient], onwardTransmission,
                                                      sampleTimes, maxTime, sdInfTime, ptLocation);

            //proposedInfTime = infTimes[proposedPatient]; //over-ride for testing
            
            //check if proposed time invalid - if so, reject move, and move to next patient
            if(proposedInfTime == std::numeric_limits<int>::min()) {
                //reject move and exit
                chainInfTimes[i] = currentInfTimes;
                chainInfSources[i] = currentInfSources;
                chainInfSourceTypes[i] = currentInfSourceTypes;
                if (debugPt) {
                    printf("Patient %d: new infection time proposed invalid (previous %d), move rejected\n\n",
                           proposedPatient, currentInfTimes[proposedPatient]);
                }
                continue;
            }
            
            
            proposedInfTimes[proposedPatient] = proposedInfTime;
            
            
            if (debugPt) {
                printf("Patient %d: new infection time proposed: %d (previous %d) [sample time %d]\n",
                        proposedPatient, proposedInfTime, currentInfTimes[proposedPatient], sampleTimes[proposedPatient]);
            }
            
            
            //propose a new source based on new infection time
            // the only changes by this time point are to the infection time of the proposed patient
            // and so can use currentInfTimes, currentRecTimes, and currentSporePatientI to determine the new sources, as the proposed patient specifically excluded in the function
            Src proposedSource = proposeConditionalSource(proposedPatient, proposedInfTime, currentInfTimes, sampleTimes, currentRecTimes,
                                                                  wardLogInf, currentInfSourceTypes, currentSporePatientI, ptLocation,
                                                                  geneticDist, nWards, nInfPatients, hospitalWards, currentParm);
            
            proposedInfSources[proposedPatient] = proposedSource.srcIndex;
            proposedInfSourceTypes[proposedPatient] = proposedSource.srcRoute;
            
            
            //Hastings ratio calculation - helper - probability of proposing new source (current --> proposed)
            double pChooseSourceProposed = proposedSource.srcP;
            
            //calculate the probability of choosing the current source with the current sampling time (proposed --> current)
            // i.e. assuming that next move is to reverse the one above
            // as only thing that changes is proposed patient can still use currentInfTimes, currentRecTimes, currentSporePatientI
            double pChooseSourceCurrent;
            SrcList srcProbsCurrent = getSourceProb(proposedPatient, currentInfTimes[proposedPatient], currentInfTimes, sampleTimes, currentRecTimes,
                                  wardLogInf, currentInfSourceTypes, currentSporePatientI, ptLocation,
                                  geneticDist,
                                  nWards, nInfPatients, hospitalWards, currentParm);
            
            for (int ii=0; ii<srcProbsCurrent.sourceList.size(); ii++) {
                //iterate over current sources
                if(srcProbsCurrent.sourceList[ii] == currentInfSources[proposedPatient] & srcProbsCurrent.sourceTypeList[ii] == currentInfSourceTypes[proposedPatient]) {
                    pChooseSourceCurrent = srcProbsCurrent.sourceProbabilities[ii];
                    break;
                }
            }
            
            
            //log hastings ratio for conditional algorithm, prob(-->current) / prob(-->proposed)
            // p(a-->b) = prob(choose Tinf) * p(choose src | Tinf)
            // prob(choose Tinf) cancels top and bottom
            hastingsRatio += log(pChooseSourceCurrent) - log(pChooseSourceProposed);
            
            
            if (debugPt) {
                printf("Patient %d: new source %d (route %d) [previous %d (%d)] - Hastings ratio: %0.4f\n", proposedPatient,
                       proposedSource.srcIndex, proposedSource.srcRoute,
                       currentInfSources[proposedPatient], currentInfSourceTypes[proposedPatient],
                       hastingsRatio);
            }
            
            //update list of onward transmissions
            //remove old onward transmission
            int oldSource = currentInfSources[proposedPatient];
            if (oldSource!=-1) { //if old source is a patient
                int oldInfTime = currentInfTimes[proposedPatient];
                trans oldTrans;
                oldTrans.t = oldInfTime;
                oldTrans.victim = proposedPatient;
                
                proposedOnwardTransmission[oldSource].erase(remove(proposedOnwardTransmission[oldSource].begin(),
                                                           proposedOnwardTransmission[oldSource].end(),
                                                           oldTrans),
                                                    proposedOnwardTransmission[oldSource].end());
                
            }
            
            //add new onward transmission
            if (proposedInfSources[proposedPatient]!=-1) { //if new source is a patient
                trans newTrans;
                newTrans.t = proposedInfTime;
                newTrans.victim = proposedPatient;
                newTrans.srcType = proposedInfSourceTypes[proposedPatient];
                newTrans.ward = ptLocation[proposedPatient][proposedInfTime];
                
                proposedOnwardTransmission[proposedInfSources[proposedPatient]].push_back(newTrans);
            }
            
            
            //get limits for ward I and sporeForceSummary update
            //first time could have set spores under current or proposed scenario
            int minPtTime = min({currentInfTimes[proposedPatient], proposedInfTime});
            if(minPtTime<0) {
                minPtTime = 0;
            }
            
            //get vector of wards visited
            set<int> wardsToUpdate;
            for (int t=minPtTime; t<=maxTime; t++) {
                wardsToUpdate.insert(ptLocation[proposedPatient][t]);
            }
            wardsToUpdate.erase(-1);
            
            
            //determine proposed infectious individuals on wards under new infection time
            updateWardI(proposedWardI, maxTime, minPtTime, wardsToUpdate, proposedInfTimes, currentRecTimes, wardLogInf);
            
            //update sporePatientI to reflect new infection times (as ward discharges while infectious may have changed)
            updateSporePatientI(proposedSporePatientI, proposedPatient, maxTime, nWards, proposedInfTimes, currentRecTimes, ptLocation);
            
            //update sporeForceSummary
            getSporeForceSummary(proposedSporeForceSummary, proposedSporePatientI, maxTime, minPtTime, wardsToUpdate, nInfPatients, proposedInfTimes, ptLocation, currentParm);
            
            //debugging code
            if(debugPt) {
                double currentLLSample = llSample(nInfPatients, currentInfTimes, sampleTimes, currentParm);
                double proposedLLSample = llSample(nInfPatients, proposedInfTimes, sampleTimes, currentParm);
                printf("LL sample, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLSample, proposedLLSample, proposedLLSample-currentLLSample);
                
                double currentLLGenetic = llGenetic(currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, nInfPatients, currentParm);
                double proposedLLGenetic = llGenetic(proposedInfTimes, sampleTimes, proposedInfSources, proposedInfSourceTypes, geneticDist, nInfPatients, currentParm);
                printf("LL genetic, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLGenetic, proposedLLGenetic, proposedLLGenetic-currentLLGenetic);
                
                double currentLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporePatientI, currentSporeForceSummary,
                                                wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
                double proposedLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, proposedInfTimes, proposedInfSourceTypes, proposedInfSources, proposedSporePatientI, proposedSporeForceSummary, wardLogInf, wardLogNeverInf, inPtDays, ptLocation, proposedWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
                printf("LL transmission, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLTrans, proposedLLTrans, proposedLLTrans-currentLLTrans);
                
                double currentLLRecovery = llRecover(nInfPatients, sampleTimes, currentRecTimes, currentParm);
                double proposedLLRecovery = llRecover(nInfPatients, sampleTimes, currentRecTimes, currentParm);
                printf("LL recovery, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLRecovery, proposedLLRecovery, proposedLLRecovery-currentLLRecovery);
                
                double prior = getPrior(currentParm);
                
                double checkCurrentLL = currentLLSample + currentLLGenetic + currentLLTrans + currentLLRecovery + prior;
                double checkProposedLL = proposedLLSample + proposedLLGenetic + proposedLLTrans + proposedLLRecovery + prior;
                double diff = checkProposedLL - checkCurrentLL;
                printf("check LL current: %0.4f\ncheck LL proposed %0.4f\ndiff: %0.4f\nhastings ratio: %0.4f\n", checkCurrentLL, checkProposedLL, diff, hastingsRatio);
                
                
            }
            

            //2. Compute the probability of accepting the proposed jump.
            proposedLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, proposedInfTimes, sampleTimes, currentRecTimes, proposedInfSources, proposedInfSourceTypes,
                                    proposedSporePatientI, proposedSporeForceSummary,
                                    wardLogInf, wardLogNeverInf, inPtDays, ptLocation, proposedWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
            
            probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
            
            
            //3. Accept proposed jump with probability pAccept
            if ( log(runif(0,1)) < probAccept ) {
                if (debugPt) printf("Move accepted (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                //accept the jump and log the change
                chainInfTimes[i] = proposedInfTimes;
                chainInfSources[i] = proposedInfSources;
                chainInfSourceTypes[i] = proposedInfSourceTypes;
                
                //increment accepted counter
                nAccepted[14] += 1;
                
                //update onward transmission log
                onwardTransmission = proposedOnwardTransmission;
                
                //set current values to the proposed values
                currentInfTimes = proposedInfTimes;
                currentInfSources = proposedInfSources;
                currentInfSourceTypes = proposedInfSourceTypes;;
                currentWardI = proposedWardI;
                currentLL = proposedLL;
                
                //currentSporePatientI = proposedSporePatientI; //this is costly - just update the bit we need to
                updateSporePatientI(currentSporePatientI, proposedPatient, maxTime, nWards, proposedInfTimes, currentRecTimes, ptLocation);
                //update currentSporeI
                updateWardI(currentWardI, maxTime, minPtTime, wardsToUpdate, proposedInfTimes, currentRecTimes, wardLogInf);
                
                
                currentSporeForceSummary = proposedSporeForceSummary;
                
                //printf("Current LL check %f\n", currentLLCheck);
                //printf("Current LL updated to %f\n", currentLL);

            }
            else {
                if (debugPt) printf("Move rejected (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                // reject the proposed jump, stay at current position
                chainInfTimes[i] = currentInfTimes;
                chainInfSources[i] = currentInfSources;
                chainInfSourceTypes[i] = currentInfSourceTypes;
                
                //reset proposedSporePatientI back to currentSporePatientI
                updateSporePatientI(proposedSporePatientI, proposedPatient, maxTime, nWards, currentInfTimes, currentRecTimes, ptLocation);
                //reset proposedWardI
                updateWardI(proposedWardI, maxTime, minPtTime, wardsToUpdate, currentInfTimes, currentRecTimes, wardLogInf);
                
                
            }
            
            
            
        } //end of loop over proposed patients to update infection times and sources
        
        
        // UPDATE RECOVERY TIMES //
        
        // Hastings ratio here = 0, as updates are proposed symmetrically
        double hastingsRatio = 0;
        
        //loop over all patients to be updated
        for (int proposedPatient : patientUpdateList) {
            
            debugPt = false;
            if (proposedPatient==-1) debugPt= true;
            
            //set proposed recovery times to current times, and spore similarlly
            proposedRecTimes = currentRecTimes;
            proposedSporePatientI = currentSporePatientI;
            
            //propose new recovery time
            int sdRecTime = 5;
            int proposedRecTime = proposeRecoveryTime(proposedPatient, currentRecTimes[proposedPatient],
                                                      onwardTransmission, sampleTimes, maxTime, minTime, sdRecTime, ptLocation, currentInfTimes);
            
            //proposedRecTime = recoveryTimes[proposedPatient];//temp OVER-RIDE for TESTING
            
            if(proposedRecTime == std::numeric_limits<int>::min()) {
                //reject move and exit
                chainRecTimes[i] = currentRecTimes;
                if (debugPt) {
                    printf("Patient %d: new recovery time proposed invalid (previous %d), move rejected\n\n",
                           proposedPatient, currentRecTimes[proposedPatient]);
                }
                continue;
            }
            
            if (debugPt) {
                printf("Patient %d: new recovery time proposed: %d (previous %d) [sample time %d]\n",
                       proposedPatient, proposedRecTime, currentRecTimes[proposedPatient], sampleTimes[proposedPatient]);
            }
            
            proposedRecTimes[proposedPatient] = proposedRecTime;
            

            
            
            //get limits for wardI and sporeForceSummary update
            //first time could have set spores under current or proposed scenario
            int minPtTime = currentInfTimes[proposedPatient];
            if(minPtTime<0) {
                minPtTime = 0;
            }
            
            //get vector of wards visited
            set<int> wardsToUpdate;
            for (int t=minPtTime; t<=maxTime; t++) {
                wardsToUpdate.insert(ptLocation[proposedPatient][t]);
            }
            wardsToUpdate.erase(-1);
            
            //update wardI
            updateWardI(proposedWardI, maxTime, minPtTime, wardsToUpdate, currentInfTimes, proposedRecTimes, wardLogInf);
            
            //update sporePatientI to reflect new recovery time
            updateSporePatientI(proposedSporePatientI, proposedPatient, maxTime, nWards, currentInfTimes, proposedRecTimes, ptLocation);
            
            //update sporeForceSummary
            getSporeForceSummary(proposedSporeForceSummary, proposedSporePatientI, maxTime, minPtTime, wardsToUpdate, nInfPatients, currentInfTimes, ptLocation, currentParm);
            
            
            if(debugPt) {
                double currentLLRecovery = llRecover(nInfPatients, sampleTimes, currentRecTimes, currentParm);
                double proposedLLRecovery = llRecover(nInfPatients, sampleTimes, proposedRecTimes, currentParm);
                printf("LL recovery, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLRecovery, proposedLLRecovery, proposedLLRecovery-currentLLRecovery);
            }
            
            
            //2. Compute the probability of accepting the proposed jump.
            proposedLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, sampleTimes, proposedRecTimes, currentInfSources, currentInfSourceTypes,
                                    proposedSporePatientI, proposedSporeForceSummary,
                                    wardLogInf, wardLogNeverInf, inPtDays, ptLocation, proposedWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
            
            probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
            
            
            //3. Accept proposed jump with probability pAccept
            if ( log(runif(0,1)) < probAccept ) {
                if (debugPt) printf("Move accepted (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                //accept the jump and log the change
                chainRecTimes[i] = proposedRecTimes;
                
                //increment accepted counter
                nAccepted[15] += 1;
                
                //set current values to the proposed values
                currentRecTimes = proposedRecTimes;
                currentWardI = proposedWardI;
                currentLL = proposedLL;
                
                //currentSporePatientI = proposedSporePatientI; //again this is costly - just update the bit we need to
                updateSporePatientI(currentSporePatientI, proposedPatient, maxTime, nWards, currentInfTimes, proposedRecTimes, ptLocation);
                //update wardI
                updateWardI(currentWardI, maxTime, minPtTime, wardsToUpdate, currentInfTimes, proposedRecTimes, wardLogInf);
                
                currentSporeForceSummary = proposedSporeForceSummary;
                
                //printf("Current LL check %f\n", currentLLCheck);
                //printf("Current LL updated to %f\n", currentLL);
                
            }
            else {
                if (debugPt) printf("Move rejected (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                // reject the proposed jump, stay at current position
                chainRecTimes[i] = currentRecTimes;;
                
                //reset proposedSporePatientI back to currentSporePatientI
                updateSporePatientI(proposedSporePatientI, proposedPatient, maxTime, nWards, currentInfTimes, currentRecTimes, ptLocation);
                //reset proposed wardI
                updateWardI(proposedWardI, maxTime, minPtTime, wardsToUpdate, currentInfTimes, currentRecTimes, wardLogInf);
            }
            
            
        } //end of loop over proposed patietns to update


        
        
        // DISRUPTION STEP TO DATA AUGMENTATION - UPDATES INFECTION TIMES AS A BLOCK FOR A SUBSET OF CASES TO ALLOW REVERSAL OF TRANSMISSION DIRECTIONS //
        int useDisruption = 0;
        if(i>20) useDisruption = 1;
        if (useDisruption == 1) {
            int debugDisrupt = 0;
            bool saveDebugPt = debugPt;
            if (debugDisrupt) {
                debugPt = true;
                printf("\n\nStarting disruption step, at iteration %d\n", i);
            }
            

            // disruption step - enables reversal of transmission directions otherwise highly unlikely using the main data augmentation steps
            // 1. select case a random, Z
            // 2. determine all cases that have Z has their direct, or indirect transmission source, together with Z,
            //      work down the transmission chain, until setSizeMax nodes (including Z) are included, let this be the set C
            int setSizeMax = 10;
            // 3. break all transmission links into nodes from the set C
            // 4. sample new infection times conditional on epsilon for each case within set C, record the component of the hastings ratio for this change
            // 5. working forwards in time using the new infection times determine the source of each case
            // 6. each change in source generates a component of the hastings ratio as before
            // 7. calculate a new likelhood, an acceptance ratio, and accept or reject all updates as a block
            
            //rationale -
            // include the source of the parent node in the changes
            //    if current state A --> B
            //    (true B --> A)
            //
            //    Choose A
            //    To support reversal need to resample source of both
            //
            //
            //    X --> Z --> A --> B
            //
            //    Choose A
            //    Keep
            //    X --> Z
            //    Resample infection sources and times for A and B
            
            // only need to remove transmission into nodes
            // --> Z --> B --> C --> D
            
            // if do this can limit, e.g. to two
            // --> Z --> B - drop links otherwise keep: Z  B --> C --> D
            
            proposedInfTimes = currentInfTimes;
            proposedInfSources = currentInfSources;
            proposedInfSourceTypes = currentInfSourceTypes;
            vector<vector<trans>> proposedOnwardTransmission = onwardTransmission;
            
            
            
            
            //0. find a list of nodes which are sources of transmission, restrict to only include these nodes
            vector<int> transmissionPatients;
            for (int patient =0; patient<nInfPatients; patient++) {
                for (trans transmission : proposedOnwardTransmission[patient]) {
                    if (transmission.t < maxTime) {
                        transmissionPatients.push_back(patient);
                        break; //only need to continue loop as far as first onward transmission
                    }
                }
            }

            if(transmissionPatients.size()>0) {
                

                //1. select a case at random, Z
                int parentNodeIndex = floor(runif(0, transmissionPatients.size()));
                int parentNode = transmissionPatients[parentNodeIndex];
                
                if (debugDisrupt) printf("Parent node for update %d\n", parentNode);
                
                //2. determine all cases that have Z has their direct, or indirect transmission source, together with Z, the set C
                //restrict size of set C to max of setSizeMax
                vector<int> nodeSet, nodeSetPrelim;
                nodeSetPrelim.push_back(parentNode); //add the parent node to preliminary set from which sources will be sought
                while(nodeSetPrelim.size()>0 & nodeSet.size()<=setSizeMax) {
                    //while there remain nodes in the preliminary set seek onward transmissions
                    //loop through the infection sources looking for parent node in the vector
                    for(int victim=0; victim<nInfPatients; victim++) {
                        if(proposedInfSources[victim] == nodeSetPrelim[0]) {
                            nodeSetPrelim.push_back(victim);
                        }
                    }
                    //add source for this iteration to the final node set
                    nodeSet.push_back(nodeSetPrelim[0]);
                    
                    //TEMP ERROR CHECKING
                    //if(nodeSetPrelim[0]==23) debugDisrupt=1;
                    
                    
                    
                    //and remove it from the preliminary set, i.e. remove the first item
                    nodeSetPrelim.erase(nodeSetPrelim.begin());
                }
                if (debugDisrupt) {
                    printf("Disruption node set: ");
                    for (int node :nodeSet) {
                        printf("%d, ", node);
                    }
                    printf("\n");
                }
                

                // 3. break all transmission links into nodes from the set C
                for (int node : nodeSet) {
                    //for all nodes set the source of transmission to be background for now, will return to later
                    proposedInfSources[node] = -1;
                    //remove all record of onward transmission into these nodes
                    
                    for(int patient=0;patient<nInfPatients;patient++) { //iterate over all possible sources of infection for node
                        trans removeTrans;
                        bool sourceFound = false;
                        for (trans transmission : proposedOnwardTransmission[patient]) { //iterate over each transmission event for each source
                            if (transmission.victim == node) {
                                removeTrans = transmission; //if find the source of transmission for this node, then exit
                                sourceFound = true;
                                break;
                            }
                        }
                        if (sourceFound) {
                            proposedOnwardTransmission[patient].erase(remove(proposedOnwardTransmission[patient].begin(),
                                                                                proposedOnwardTransmission[patient].end(),
                                                                                removeTrans),
                                                                      proposedOnwardTransmission[patient].end());
                            break;
                        }
                    }
                }

                
                // 4. sample new infection times - this is done symmetrically so Hastings ratio for time update is 0
                for (int node : nodeSet) {
                    
                    int sdDisruptInfTime = 4;
                    int proposedInfTime = proposeInfectionTime(node, currentInfTimes[node], proposedOnwardTransmission,
                                                           sampleTimes, maxTime, sdDisruptInfTime, ptLocation);
                    if(proposedInfTime == std::numeric_limits<int>::min()) {
                        proposedInfTime = currentInfTimes[node]; //keep current infection time
                    }
                    
                    proposedInfTimes[node] = proposedInfTime;
                }
                
                if (debugDisrupt) {
                    printf("New infection times - \n");
                    for (int node :nodeSet) {
                        printf("\tNode %d, infection time: %d (previous %d)\n", node, proposedInfTimes[node], currentInfTimes[node]);
                    }
                }
                
                
                // 4B) sample new recovery times as well, to avoid recovery times being fixed by on-ward transmission requirements
                proposedRecTimes = currentRecTimes;
                for(int node : nodeSet) {
                    int sdDisruptRecTime = 5;
                    int proposedRecoveryTime = proposeRecoveryTime(node, currentRecTimes[node], proposedOnwardTransmission,
                                                                   sampleTimes, maxTime, minTime, sdDisruptRecTime, ptLocation, proposedInfTimes);
                    if(proposedRecoveryTime == std::numeric_limits<int>::min()) {
                        //keep current infection time
                        proposedRecoveryTime = currentRecTimes[node];
                    }
                    proposedRecTimes[node] = proposedRecoveryTime;
                }
                
                
                // 5. working forwards in time using the new infection times determine the source of each case
                // 6. each change in source generates a component of the hastings ratio as before
                
                // update all as a block
                // for proposedInfTimes, also need to update sporePatientI based on these - unlike above where only updating one case
                proposedWardI = getWardI(maxTime, nWards, proposedInfTimes, proposedRecTimes, wardLogInf);
                getSporePatientI(proposedSporePatientI, nInfPatients, maxTime, nWards, proposedInfTimes, proposedRecTimes, ptLocation);
                getSporeForceSummary(proposedSporeForceSummary, proposedSporePatientI, maxTime, minTime, allWards, nInfPatients, proposedInfTimes, ptLocation, currentParm);
                
                //currentInfSourceTypes - is used to determine genetic likelihood - use this for all updates current --> proposed updates
                // and converse (proposedInfSourceTypes) for proposed --> current updates
                
                double hastingsRatio = 0;
                
                for (int node : nodeSet) {
                    
                    //get new source for each case
                    Src proposedSource = proposeConditionalSource(node, proposedInfTimes[node], proposedInfTimes,
                                                sampleTimes, proposedRecTimes, wardLogInf, currentInfSourceTypes, proposedSporePatientI, ptLocation,
                                                geneticDist, nWards, nInfPatients, hospitalWards, currentParm);
                    proposedInfSources[node] = proposedSource.srcIndex;
                    proposedInfSourceTypes[node] = proposedSource.srcRoute;
                    
                    //record onward transmission if from non-background
                    if (proposedInfSources[node]!=-1) {
                        trans newTrans;
                        newTrans.t = proposedInfTimes[node];
                        newTrans.victim = node;
                        newTrans.srcType = proposedInfSourceTypes[node];
                        newTrans.ward = ptLocation[node][proposedInfTimes[node]];
                        proposedOnwardTransmission[proposedInfSources[node]].push_back(newTrans);
                    }
                    
                    
                    //log hastings ratio for conditional algorithm, prob(-->current) / prob(-->proposed)
                    // p(a-->b) = prob(choose Tinf) * p(choose src | Tinf)
                    // prob(choose Tinf) cancels top and bottom
                    
                    //Hastings ratio denominator calculation - helper - probability of proposing new source
                    double pChooseSourceProposed = proposedSource.srcP;
                    hastingsRatio -= log(pChooseSourceProposed);
                    
                    if (debugDisrupt) {
                        printf("Patient %d: new source %d (route %d) [previous %d (%d)] - cumulative Hastings ratio: %0.4f\n", node,
                               proposedSource.srcIndex, proposedSource.srcRoute,
                               currentInfSources[node], currentInfSourceTypes[node],
                               hastingsRatio);
                    }
                    
                }

                //now find hastings ratio numerator, as a block, using proposedInfSourceTypes
                
                for (int node : nodeSet) {
                    // determine hastings ratio numerator: for (proposed --> current) for all updated cases
                    // using currentInfTimes, currentRecTimes, currentSporePatientI
                    double pChooseSourceCurrent;
                    SrcList srcProbsCurrent = getSourceProb(node, currentInfTimes[node], currentInfTimes,
                                sampleTimes, currentRecTimes, wardLogInf, proposedInfSourceTypes, currentSporePatientI, ptLocation,
                                geneticDist, nWards, nInfPatients, hospitalWards, currentParm);
                    
                    for (int ii=0; ii<srcProbsCurrent.sourceList.size(); ii++) {
                        //iterate over possible sources until find current source
                        if(srcProbsCurrent.sourceList[ii] == currentInfSources[node] & srcProbsCurrent.sourceTypeList[ii] == currentInfSourceTypes[node]) {
                            pChooseSourceCurrent = srcProbsCurrent.sourceProbabilities[ii];
                            break;
                        }
                    }
                    hastingsRatio += log(pChooseSourceCurrent);
                }

                
                
                // 7. calculate a new likelhood, an acceptance ratio, and accept or reject all updates as a block
                
                //determine proposed infectious individuals on wards under new infection time
                
                
                proposedLL = targetDist(wardEver, hospitalWards, ward2Hospital, hospitalWardList, proposedInfTimes, sampleTimes, proposedRecTimes, proposedInfSources, proposedInfSourceTypes,
                                        proposedSporePatientI, proposedSporeForceSummary,
                                        wardLogInf, wardLogNeverInf, inPtDays, ptLocation, proposedWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist, currentParm);
                
                probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
                
                if(debugDisrupt) {
                    double currentLLSample = llSample(nInfPatients, currentInfTimes, sampleTimes, currentParm);
                    double proposedLLSample = llSample(nInfPatients, proposedInfTimes, sampleTimes, currentParm);
                    printf("LL sample, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLSample, proposedLLSample, proposedLLSample-currentLLSample);
                    
                    double currentLLGenetic = llGenetic(currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, nInfPatients, currentParm);
                    double proposedLLGenetic = llGenetic(proposedInfTimes, sampleTimes, proposedInfSources, proposedInfSourceTypes, geneticDist, nInfPatients, currentParm);
                    printf("LL genetic, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLGenetic, proposedLLGenetic, proposedLLGenetic-currentLLGenetic);
                    
                    double currentLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporePatientI, currentSporeForceSummary,
                                                    wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
                    double proposedLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, proposedInfTimes, proposedInfSourceTypes, proposedInfSources,
                                                     proposedSporePatientI, proposedSporeForceSummary,
                                                     wardLogInf, wardLogNeverInf, inPtDays, ptLocation, proposedWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
                    printf("LL transmission, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLTrans, proposedLLTrans, proposedLLTrans-currentLLTrans);
                    
                    printf("Acceptance probability: %0.4f\n", probAccept);
                    
                    
                }
                
                
                if ( log(runif(0,1)) < probAccept ) {
                    if (debugDisrupt) printf("Move accepted (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                    //accept the jump and log the change
                    chainInfTimes[i] = proposedInfTimes;
                    chainInfSources[i] = proposedInfSources;
                    chainInfSourceTypes[i] = proposedInfSourceTypes;
                    chainRecTimes[i] = proposedRecTimes;
                    
                    //increment distruption moves accepted counter
                    nAccepted[16] += 1;
                    
                    //update onward transmission log
                    onwardTransmission = proposedOnwardTransmission;
                    
                    //set current values to the proposed values
                    currentInfTimes = proposedInfTimes;
                    currentInfSources = proposedInfSources;
                    currentInfSourceTypes = proposedInfSourceTypes;
                    currentRecTimes = proposedRecTimes;
                    currentWardI = proposedWardI;
                    currentLL = proposedLL;
                    currentSporePatientI = proposedSporePatientI;
                    currentSporeForceSummary = proposedSporeForceSummary;
                    
                }
                else {
                    if (debugDisrupt) printf("Move rejected (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                    // reject the proposed jump, stay at current position, already saved in main data augmentation step
                    
                    //reset proposed sporePatientI back to currentSporePatientI, this is done for other variables at the start of the data augmentation step
                    proposedSporePatientI = currentSporePatientI;
                }
                
                if (debugDisrupt) {
                    debugPt = saveDebugPt;
                    printf("\n\n");
                }
            }
            
        } //end of disruption section
        
        
        double currentLLSample = llSample(nInfPatients, currentInfTimes, sampleTimes, currentParm);
        double currentLLGenetic = llGenetic(currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, nInfPatients, currentParm);
        double currentLLTrans = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporePatientI,
                                        currentSporeForceSummary, wardLogInf, wardLogNeverInf, inPtDays, ptLocation, currentWardI, nInfPatients, nNeverInfPatients, nWards, maxTime, currentParm);
        double currentLLRecovery = llRecover(nInfPatients, sampleTimes, currentRecTimes, currentParm);

        //log posterior
        chain[i][14] = currentLL;
        //log trans likelihood
        chain[i][15] = currentLLTrans;
        //log genetic likelihood
        chain[i][16] = currentLLGenetic;
        // log samling likelihood
        chain[i][17] = currentLLSample;
        //log recovery likelihood
        chain[i][18] = currentLLRecovery;
        
        
    } //end of MCMC steps
} //end MCMC function






//options for arugement parser
struct Arg: public option::Arg
{
    static void printError(const char* msg1, const option::Option& opt, const char* msg2)
    {
        fprintf(stderr, "ERROR: %s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
    }
    
    static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
    {
        if (option.arg != 0 && option.arg[0] != 0)
            return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires a non-empty argument\n");
        return option::ARG_ILLEGAL;
    }
    
    static option::ArgStatus Optional(const option::Option& option, bool msg)
    {
        return option::ARG_OK;
    }
    
    static option::ArgStatus Numeric(const option::Option& option, bool msg)
    {
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
        if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires a numeric argument\n");
        return option::ARG_ILLEGAL;
    }
};

enum  optionIndex { UNKNOWN, HELP, PATH, ITER, SEED};
const option::Descriptor usage[] =
{
    {UNKNOWN, 0,"" , ""    , option::Arg::None, "USAGE: test_program [options]\n\n"
        "Options:" },
    {HELP,    0, "" , "help", option::Arg::None, "  --help  \tPrint usage and exit." },
    {PATH,    0, "p", "path", Arg::NonEmpty, "  --path, -p  \tPath to the input files." },
    {ITER,    0, "i", "iter", Arg::Numeric, "  --iter, -i  \tNumber of MCMC iterations." },
    {SEED,    0, "s", "seed", Arg::Optional, "  --seed, -s  \tRandom seed for MCMC [default==Random], if set to 0 then random seed is used." },
    {0,0,0,0,0,0}
};




int main(int argc, const char * argv[]) {
    
    //option parsing and error checking
    argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
    option::Stats stats(usage, argc, argv);
    option::Option* options = (option::Option*)calloc(stats.options_max, sizeof(option::Option));
    option::Option* buffer  = (option::Option*)calloc(stats.buffer_max,  sizeof(option::Option));
    option::Parser parse(usage, argc, argv, options, buffer);
    
    if (parse.error())
        return 1;
    
    if (options[HELP] || argc == 0) {
        option::printUsage(std::cout, usage);
        return 0;
    }
    
    bool err = false;
    for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next()) {
        cout << "Unknown option: " << string(opt->name,opt->namelen) << "\n";
        err = true;
    }
    
    for (int i = 0; i < parse.nonOptionsCount(); ++i) {
        std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";
        err = true;
    }
    
    if(err) {
        cout << "\n";
        option::printUsage(std::cout, usage);
        return 1;
    }
    
    
    //read in options
    string path;
    int steps;
    if(options[PATH]) path = options[PATH].arg;
    if(options[ITER]) steps = stoi(options[ITER].arg);
    
    //read in random seed
    int seed, rndSeed;
    if(options[SEED]) seed = stoi(options[SEED].arg);
    if(seed==0) {
        rndSeed = (int)time(NULL);
    }
    else {
        rndSeed = seed;
    }
    
    set_seed(rndSeed, rndSeed*2-1);
    srand(rndSeed);
    
    
    
    
    // **** IMPORT DATA **** //
    //read in epi data CSV files
    
    //vectors for sample times for infected cases
    vector <int> sampleTimes;
    
    //vectors for infection times and recovery times, infection sources, and source types for infected cases
    vector <int> infTimes, recoverTimes;
    vector <int> infSources; //recorded using number source indentifiers
    vector <int> infSourceTypes; //hospital background = 0, ward =1, hospital = 2, comm background = 3, started infection = 4
    
    
    int nInfPatients; //number of infected patients
    int nNeverInfPatients; //number of patients never infected, but contributing to admission logs
    
    unordered_map<string,int> ptLookup; //lookup that converts patient identifers as strings into integers numbered from zero
    unordered_map<int,string> ptLookupRev; //reverse for export
    
    //import details for infected patients
    string filePath = path + "/input/patientLog.csv";
    importPatientLog(filePath, ptLookup, ptLookupRev, infTimes, infSources, infSourceTypes, sampleTimes, recoverTimes, nInfPatients, nNeverInfPatients);
    
    
    //set up ward logs
    vector<vector<vector<int>>> wardLogInf;
    vector<vector<int>> wardLogNeverInf;
    //    wardLogInf -  // wardLogInf[time][ward] = {patients...} for infected patients
    //    wardLogNeverInf - // wardLogNeverInf[time][ward] = int for count of never infected patients
    
    //ward lookup that converts text ward names into integers starting at zero, hospitalLookup likewise
    unordered_map<string,int> wardLookup;
    unordered_map<string,int> hospitalLookup;
    
    //last time to track
    int maxTime;
    
    //number of wards
    int nWards;
    
    //vector of wards in the same hospital, by ward, excluding the ward itself - used for finding cases on all other words in the same hospital
    vector<vector<int>> hospitalWards;
    vector<int> ward2Hospital; //lookup of hopsitals indexed by ward id
    vector<vector<int>> hospitalWardList; //list of wards in each hospital
    
    //import ward admission details
    filePath = path + "/input/wardLog.csv";
    importWardLog(filePath, hospitalLookup, wardLookup, ptLookup, wardLogInf, wardLogNeverInf, maxTime, nWards, sampleTimes, hospitalWards, ward2Hospital, hospitalWardList);
    
    //read in genetic data
    string filePathGenetic = path + "/input/geneticDistances_snps.txt";
    vector<vector<double>> geneticDist; //pairwise genetic distances, indexed using integers from ptLookup
    importGeneticData(filePathGenetic, geneticDist, ptLookup, nInfPatients);

 
    
    // **** RUN MCMC **** //
    // MCMC setup
    //clock_t startMCMC,finishMCMC;
    //startMCMC = clock();
    struct timeval startMCMC, endMCMC;
    
    gettimeofday(&startMCMC, NULL);
    double deltaMCMC;

    
    
    //starting values for parameters
    Parm startParm;
    startParm.betaBgroundHosp = 0.0005;
    startParm.betaWard = 0.0005;
    startParm.betaHosp = 0.0005;
    startParm.sampleSize = 3;
    startParm.sampleMu = 30;
    startParm.directNe =  1;
    startParm.introNe = 500;
    startParm.mu = 2/365.25;
    startParm.probStartInfLogit = logit(0.01);
    startParm.betaComm = 0.0005;
    startParm.sporeProbLogit = 0.5;
    startParm.recSize = 3;
    startParm.recMu = 90;
    startParm.sporeMultiplier = 0.4;
    
    //starting values for standard deviation for metropolis updates to parameters
    Parm startSigma;
    startSigma.betaBgroundHosp = 0.0001;
    startSigma.betaWard = 0.0001;
    startSigma.betaHosp = 0.0001;
    startSigma.sampleSize = 0.1;
    startSigma.sampleMu = 0.1;
    startSigma.directNe =  0.1;
    startSigma.introNe = 100;
    startSigma.mu = 2/365.25/10;
    startSigma.probStartInfLogit = 0.05;
    startSigma.betaComm = 0.0001;
    startSigma.sporeProbLogit = 0.1;
    startSigma.recSize = 0.1;
    startSigma.recMu = 1;
    startSigma.sporeMultiplier = 0.05;

    
    //2d vector for results = columns: "beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm"; rows = each iteration
    //2d vectors for infection times, infection sources, and source types, columns for each patient, row = each iteration
    vector<Parm> chain;
    vector<vector<int>> chainInfTimes, chainInfSources, chainInfSourceTypes, chainRecTimes;
    
    //set up sizes of chain for parameters and various times as 2d vectors
    chain.resize(steps);
    chainInfTimes.resize(steps); chainInfSources.resize(steps);
    chainInfSourceTypes.resize(steps); chainRecTimes.resize(steps);
    
    for (int i = 0; i<steps; i++) {
        chainInfTimes[i].resize(nInfPatients);
        chainInfSources[i].resize(nInfPatients);
        chainInfSourceTypes[i].resize(nInfPatients);
        chainRecTimes[i].resize(nInfPatients);

    }
    
    vector<vector<vector<int>>> inPtDays = getInPtDays(nInfPatients, maxTime, nWards, wardLogInf);
    vector<vector<int>> ptLocation = getPtLocation(nInfPatients, maxTime, nWards, inPtDays);
    vector<vector<int>> wardI = getWardI(maxTime, nWards, infTimes, recoverTimes, wardLogInf);
    
    
    
    //runTest(infTimes, sampleTimes, recoverTimes, startParm,
    //        infSources, infSourceTypes,
    //        sporePatientI, sporeForce, sporeForceSummary, geneticDist, geneticMap, wardLog,
    //        inPtDays, ptLocation, wardI, nPatients, nWards, maxTime);
    
    
    doMCMC(chain, chainInfTimes, chainRecTimes, chainInfSources, chainInfSourceTypes, steps, startParm, startSigma,
           sampleTimes, wardLogInf, wardLogNeverInf, nInfPatients, nNeverInfPatients, nWards, maxTime, geneticDist,
           infTimes, infSources, infSourceTypes, recoverTimes, ptLookup, wardLookup, hospitalWards, ward2Hospital, hospitalWardList);
    
    
    //export the chain to a file
    const string mcmcLog = path + "/inference"; //path for log files
    exportChain(chain, chainInfTimes, chainInfSources, chainInfSourceTypes, chainRecTimes, steps, mcmcLog, ptLookupRev);
    
    
    //report MCMC run time
    //finishMCMC = clock();			//Stop the clock.
    //double timerMCMC = (double(finishMCMC)-double(startMCMC))/CLOCKS_PER_SEC;
    //printf("MCMC Processing time: %0.4f seconds\n",timerMCMC);
    gettimeofday(&endMCMC, NULL);
    deltaMCMC = ((endMCMC.tv_sec  - startMCMC.tv_sec) * 1000000u +
             endMCMC.tv_usec - startMCMC.tv_usec) / 1.e6;
    printf("MCMC Processing time: %0.4f seconds\n",deltaMCMC);
    
    
    return 0;
}
