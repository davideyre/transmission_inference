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
              int steps, Parm startParm, vector<double> startSigma,
              vector<int> &sampleTimes,
              vector<vector<vector<int>>> &wardLog,
              int nPatients, int nWards, int maxTime,
              vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap,
              vector<int> &infTimes, vector<int> &infSources, vector<int> &infSourceTypes, vector<int> &recoveryTimes) {
    
    //note are supplied for testing purposes only -
    // vector<int> &infTimes, vector<int> &infSources, vector<int> &infSourceTypes, vector<int> &recoveryTimes
    
    
    //get a list of infected and uninfected patients - this could be dynamic if augmenting which patients are infected, but fix for now
    vector<int> infectedPatients, uninfectedPatients;
    for (int pt=0; pt<nPatients; pt++) {
        if(sampleTimes[pt]==-1) {
            uninfectedPatients.push_back(pt);
        } else {
            infectedPatients.push_back(pt);
        }
    }
    
    //number infected
    int nInfected = (int)infectedPatients.size();
    int updateN = round(nInfected*0.1); // select a subset of patients to update times for at each iteration, e.g. 10%
    bool debugPt = false;
    
    

    
    
    //reshape the ward log - inPtDays[patient][ward] = {times...} (whereas wardLog[time][ward] = {patients...})
    vector<vector<vector<int>>> inPtDays = getInPtDays(nPatients, maxTime, nWards, wardLog);
    //reshape again: ptLocation[patient][time] = wardId
    vector<vector<int>> ptLocation = getPtLocation(nPatients, maxTime, nWards, inPtDays);
    
    
    //onward transmission log, indexed by source patient (contains time of onward transmission and victim)
    vector<vector<trans>> onwardTransmission;
    onwardTransmission.resize(nPatients);
    trans defaultTrans;
    defaultTrans.t = maxTime;
    defaultTrans.victim = -1;
    defaultTrans.srcType = -1;
    defaultTrans.ward = -1;
    for(int patient=0;patient<nPatients;patient++) {
        onwardTransmission[patient].push_back(defaultTrans);
    }
    
    //set starting values for infection times and infection sources
    vector<int> startInfTimes, startInfSources, startInfSourceTypes, startRecTimes;
    
    //vector<double> startPChooseSource;
    for(int patient=0; patient<nPatients; patient++) {
        if(sampleTimes[patient]==-1) {
            //uninfected patients remain uninfected
            startInfTimes.push_back(-1);
            startInfSources.push_back(-1);
            startInfSourceTypes.push_back(-1);
            startRecTimes.push_back(-1);
        }
        else {
            //infected patients
            //sample an infection time using the starting value of epsilon
            //parameters - "beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm"
            double e  = startParm[3];
            double recSize = startParm[11];
            double recMu = startParm[12];
            double recProb = recSize / (recSize + recMu);
            int infTime = sampleTimes[patient] - rpois(e);
            
            if(infTime<1) {
                //start infected
                startInfTimes.push_back(0); //store as zero, but actually may be earlier
                startInfSources.push_back(-1); //background source
                startInfSourceTypes.push_back(4); //start infected source type
                int t_rec = sampleTimes[patient] + rnbinom(recSize, recProb); //assume initial recovery times (approx as assume for this infected at t=0)
                startRecTimes.push_back(t_rec);
            }
            else {
                startInfTimes.push_back(infTime);
                //set source of all infections as background initially, check location at sampling time
                startInfSources.push_back(-1);
                if(ptLocation[patient][infTime]>-1) {
                    startInfSourceTypes.push_back(0); //inpatient
                }
                else {
                    startInfSourceTypes.push_back(3); //community
                }
                int t_rec = sampleTimes[patient] + rnbinom(recSize, recProb);
                startRecTimes.push_back(t_rec);
            }

        } // end infected patient
    } //end patient loop
    
    //startSporeDurations = sporeDurations; // TEMP OVER_RIDE TO FIX WHEN IMPLEMENT AUGMENTATION
    

    
    //initialise chain with parameter starting values
    chain[0] = startParm;
    chainInfTimes[0] = startInfTimes;
    chainRecTimes[0] = startRecTimes;
    chainInfSources[0] = startInfSources;
    chainInfSourceTypes[0] = startInfSourceTypes;
    
    //set up acceptance counters - "beta0", "beta1", "beta2", "sampleSize", "sampleMu", "directNe", "introNe", "mu", "startInf", "betaComm", "sporeProb", "recSize", "recMu", , "augmentation moves - infection", "aug - recovery", "disruption moves"
    vector<int> nAccepted = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    //initial sd for MH updates
    vector<double> sigma = startSigma;
    
    //variables used during MCMC
    Parm currentParm; currentParm = chain[0];
    Parm proposedParm;
    double proposedValue, probAccept;
    
    //current and proposed values
    vector<int> currentInfTimes, currentRecTimes, currentInfSources, currentInfSourceTypes, proposedInfTimes, proposedRecTimes, proposedInfSources, proposedInfSourceTypes;
    vector<vector<vector<int>>> currentSporeI, proposedSporeI;
    vector<vector<double>> currentSporeForceSummary, proposedSporeForceSummary;
    vector<vector<int>> currentWardI, proposedWardI; //set up wardI log - //create array of the number of infected individuals on each ward at each time point
    
    
    //reserve memory for sporeI
    //set up size of currentSporeI
    currentSporeI.resize(maxTime+1);
    for (int i = 0; i<=maxTime; i++) {
        currentSporeI[i].resize(nWards);
        for (int ward = 0; ward<nWards; ward++) {
            currentSporeI[i][ward].resize(nPatients);
        }
    }
    //set up size of proposedSporeI
    proposedSporeI.resize(maxTime+1);
    for (int i = 0; i<=maxTime; i++) {
        proposedSporeI[i].resize(nWards);
        for (int ward = 0; ward<nWards; ward++) {
            proposedSporeI[i][ward].resize(nPatients);
        }
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
    
    
    currentWardI = getWardI(nPatients, maxTime, nWards, currentInfTimes, currentRecTimes, wardLog); //get initial values
    getSporeI(currentSporeI, infectedPatients, nPatients, maxTime, nWards, currentInfTimes, currentRecTimes, ptLocation);
    getSporeForceSummary(currentSporeForceSummary, infectedPatients, currentSporeI, maxTime, nWards, nPatients, currentInfTimes, currentParm);
    
    
   
    //calculate initial targetDist value
    double currentLL = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, currentRecTimes,
                                  currentInfSources, currentInfSourceTypes,
                                  currentSporeI, currentSporeForceSummary,
                                  wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, currentParm);
    printf("Starting target distribution value: %0.4f\n\n", currentLL);
    double proposedLL;
    
    
    //clock_t previousSystemTime = clock();
    
    struct timeval previousSystemTime;
    gettimeofday(&previousSystemTime, NULL);


    
    //generate random walk
    for(int i=1; i<steps; i++) {
        //printf("Starting iteration %d\n", i);
        
        //PARAMETER AND TUNING REPORTING
        //every 100 steps
        if(i % 100 == 0) {
            printf("Starting iteration %d\n", i);
            //tune proposal distributions - update sigma based on nAccepted
            printf("Accepted Moves\n beta0: %d\tbeta1: %d\tbeta2: %d\tsampleSize: %d\tsampleMu: %d\tdirectNe: %d\tintroNe: %d\tmu: %d\t betaComm: %d\nsporeProb: %d\tpStartInf: %d\trecSize %d\trecMu %d\n",
                   nAccepted[0], nAccepted[1], nAccepted[2], nAccepted[3], nAccepted[4], nAccepted[5], nAccepted[6], nAccepted[7], nAccepted[8], nAccepted[9], nAccepted[10], nAccepted[11], nAccepted[12]);
            for (int j=0; j<13; j++) {
                if(nAccepted[j]<20) {
                    sigma[j] = sigma[j] * 0.8;
                } else if (nAccepted[j]>40) {
                    sigma[j] = sigma[j] * 1.25;
                }
                nAccepted[j] = 0;
            }
            
            //print augmentation moves
            printf("Augmentation (infection) moves accepted %d\n", nAccepted[13]);
            nAccepted[13] = 0;
            printf("Augmentation (recovery) moves accepted %d\n", nAccepted[14]);
            nAccepted[14] = 0;
            
            //print disruption moves - cumulative total
            printf("Disruption moves accepted %d\n", nAccepted[15]);
            nAccepted[15] = 0;
            
            
            //display current parameter values
            chain[i-1].displayLog();
            
            //display current LL for different components of likelihood
            double currentLLSample = llSample(infectedPatients, currentInfTimes, sampleTimes, currentParm);
            printf("LL sample, current: %0.3f\n", currentLLSample);
            
            double currentLLGenetic = llGenetic(infectedPatients, currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
            printf("LL genetic, current: %0.3f\n", currentLLGenetic);
            
            double currentLLTrans = llTrans(infectedPatients, uninfectedPatients, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporeI, currentSporeForceSummary, wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, currentParm);
            printf("LL trans, current: %0.3f\n", currentLLTrans);
            
            double currentLLTransAvoid = llTransAvoid(infectedPatients, uninfectedPatients, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporeI, currentSporeForceSummary, wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, currentParm);
            printf("LL trans avoid, current: %0.3f\n", currentLLTransAvoid);
            
            
            double currentLLRecovery = llRecover(infectedPatients, sampleTimes, currentRecTimes, currentParm);
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
        for (int parmIndex=0; parmIndex<13; parmIndex++) {
                
            //iterate through parmameters - see Parm class for list in struct.hpp
            //1. proposal distribution, i.e. proposed change in value of beta
            proposedValue = currentParm[parmIndex] + rnorm(0 , sigma[parmIndex]);
            proposedParm = currentParm; proposedParm[parmIndex] = proposedValue;
            if(proposedValue <= 0 & parmIndex !=9 & parmIndex !=10) {
                //if invalid proposed value skip this step - probabilities on logit scale so can be negative
                chain[i][parmIndex] = chain[i-1][parmIndex];
            }
            else {
                
                if(parmIndex==9) {
                    getSporeForceSummary(proposedSporeForceSummary, infectedPatients, currentSporeI, maxTime, nWards, nPatients, currentInfTimes, proposedParm);
                    proposedLL = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                            currentSporeI, proposedSporeForceSummary,
                                            wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, proposedParm);
                }
                else {
                    //2. Compute the probability of accepting the proposed jump.
                    proposedLL = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                            currentSporeI, currentSporeForceSummary,
                                            wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, proposedParm);
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
                    //printf("Move rejected for parameter %d at iteration %d (current: %0.4f, proposed %0.4f\n", parmIndex, i, currentParm[parmIndex], proposedValue);
                }
            }
            
        } //end 13 parameter updates
        
        
        
        //over-ride - fix 9 paramters
        
        //TRUE values  - "beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm"
        
        //chain[i][0] = 0.005;
        //chain[i][1] = 0.005;
        //chain[i][2] = 0.002;
        //chain[i][3] = 10;
        //chain[i][4] = 1;
        //chain[i][5] = 10000;
        //chain[i][6] = 0.005475702;
        //chain[i][7] = 0.000001;
        //chain[i][8] = 0.000001;
        //chain[i][9] = 0.2;
        //chain[i][10] = 0.002;
        //chain[i][11] = 3;
        //chain[i][12] = 30;
        
        
        
        //over-ride for data augmentation
        /*
        currentInfTimes = infTimes;
        currentRecTimes = recoveryTimes;
        currentInfSources = infSources;
        currentInfSourceTypes = infSourceTypes;
        currentWardI = getWardI(nPatients, maxTime, nWards, currentInfTimes, currentRecTimes, wardLog);
        getSporeI(currentSporeI, infectedPatients, nPatients, maxTime, nWards, infTimes, currentRecTimes, ptLocation);
        getSporeForceSummary(currentSporeForceSummary, infectedPatients, currentSporeI, maxTime, nWards, nPatients, currentInfTimes, currentParm);
       
        //generic over-ride
    
        currentParm = chain[i];
        currentLL = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                               currentSporeI, currentSporeForceSummary,
                             wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, geneticDist, geneticMap,currentParm);
         getSporeForceSummary(currentSporeForceSummary, infectedPatients, currentSporeI, maxTime, nWards, nPatients, currentInfTimes, currentParm);
         */
        // end over-ride

        
        //DATA AUGMENTATION
        
        // get list of patients to update at this iteration
        vector<int> patientUpdateList(nInfected);
        iota(begin(patientUpdateList), end(patientUpdateList), 0); //create a list of all infected patients
        random_shuffle(begin(patientUpdateList), end(patientUpdateList)); //shuffle the patient list - UPDATE TO TAKE SEED
        patientUpdateList.resize(updateN); //keep only the first 10%

        
        for (int proposedPatientIndex : patientUpdateList) {
            
            double hastingsRatio; //store log hastings ratio
            int proposedInfTime; //store the new proposed infection time
            
            //loop through patient infection times to update
            //printf("Updating patient: %d\n", proposedPatientIndex);
            int proposedPatient = infectedPatients[proposedPatientIndex];
            
            debugPt = false;
            if (proposedPatient==-1) {
                printf("i=%d, patient=%d, current LL: %0.3f\n", i, proposedPatient, currentLL);
                double currentLLCheck = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, currentRecTimes, currentInfSources, currentInfSourceTypes,
                                                   currentSporeI, currentSporeForceSummary, wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, currentParm);
                printf("LL check: %0.3f\n", currentLLCheck);
                debugPt = true;
            }
            
            
            proposedInfTimes = currentInfTimes;
            proposedInfSources = currentInfSources;
            proposedInfSourceTypes = currentInfSourceTypes;
            
            vector<vector<trans>> proposedOnwardTransmission = onwardTransmission;
            
            
            //propose the patient was positive at admission with probability q - IMPLEMENT THIS AND HR TO GO WITH IT
            double q = 0.01;
            
            if (runif(0,1) < q) {
                //propose patient was positive at time zero, infTime=0, infSource=-1, infSourceType=4
                proposedInfTime = 0;
                proposedInfTimes[proposedPatient] = proposedInfTime;
                proposedInfSources[proposedPatient] = -1;
                proposedInfSourceTypes[proposedPatient] = 4;
                
                if (currentInfTimes[proposedPatient] == 0) {
                    hastingsRatio = 0;
                }
                else {
                    //hastings ratio - prob(-->current) / prob(-->proposed)
                    hastingsRatio = log((1-q) / q);
                }
                
            }
            
            else {
                //patient acquired infection after t=0
                //calculate hastings ratio for proposing start time >0 and then propose infection time
                if (currentInfTimes[proposedPatient] > 0) {
                    //if already have infection time after t=0
                    hastingsRatio = 0;
                    //new time - adjust using random walk
                    double delta = 2.5;
                    proposedInfTime = proposeInfectionTime(proposedPatient, currentInfTimes[proposedPatient], onwardTransmission,
                                                           sampleTimes, maxTime, delta, ptLocation);
                    
                }
                else {
                    //if currently have infection time, t=0
                    hastingsRatio = log(q / (1-q));
                    //propose infection time conditional on epsilon - this avoids poor mixing where get stuck at t=0 as infection time
                    //proposedInfTime =  proposeInfectionTimeConditional(proposedPatient, onwardTransmission, sampleTimes, maxTime, currentParm[3]);
                    
                    proposedInfTime = proposeInfectionTimeInitial(proposedPatient, onwardTransmission, sampleTimes, maxTime);
                    
                    //adjust hastings ratio - log PDF at current - log PDF at proposed
                    int currentInterval = sampleTimes[proposedPatient] - currentInfTimes[proposedPatient];
                    int proposedInterval = sampleTimes[proposedPatient] - proposedInfTime;
                    hastingsRatio += dpois(currentInterval, currentParm[3], 1) - dpois(proposedInterval, currentParm[3], 1);
                    
                }
                
                //proposedInfTime = infTimes[proposedPatient]; //over-ride for testing
                
                if(proposedInfTime==-1) { 
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
                
                
                //propose a new source
                // the only changes here are to the infection time of the proposed patient
                // and so can use currentInfTimes, currentRecTimes, and currentSporeI to determine the new sources, as the proposed patient specifically excluded in the function
                Src proposedSource = proposeConditionalSource(infectedPatients, proposedPatient, proposedInfTime, currentInfTimes, sampleTimes, currentRecTimes,
                                                                      wardLog, currentInfSourceTypes, currentSporeI, ptLocation,
                                                                      geneticDist, geneticMap, nWards, nPatients, currentParm);
                
                
                proposedInfSources[proposedPatient] = proposedSource.srcIndex;
                proposedInfSourceTypes[proposedPatient] = proposedSource.srcRoute;
                
                
                //Hastings ratio calculation - helper - probability of proposing new source (current --> proposed)
                double pChooseSourceProposed = proposedSource.srcP;
                
                //calculate the probability of choosing the current source with the current sampling time (proposed --> current)
                // i.e. assuming that next move is to reverse the one above
                // as only thing that changes is proposed patient can still use currentInfTimes, currentRecTimes, currentSporeI
                double pChooseSourceCurrent;
                SrcList srcProbsCurrent = getSourceProb(infectedPatients, proposedPatient, currentInfTimes[proposedPatient], currentInfTimes, sampleTimes, currentRecTimes,
                                      wardLog, currentInfSourceTypes, currentSporeI, ptLocation,
                                      geneticDist, geneticMap,
                                      nWards, nPatients, currentParm);
                
                for (int ii=0; ii<srcProbsCurrent.sourceList.size(); ii++) {
                    //iterate over current sources
                    if(srcProbsCurrent.sourceList[ii] ==  currentInfSources[proposedPatient]) {
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
                
                
            } //end of else vs propose positive at t=0
            
            
            
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
            
            
            //determine proposed infectious individuals on wards under new infection time
            proposedWardI = getWardI(nPatients, maxTime, nWards, proposedInfTimes, currentRecTimes, wardLog);
            //update sporeI to reflect new infection times (as ward discharges while infectious may have changed)
            proposedSporeI = currentSporeI;
            updateSporeI(proposedSporeI, proposedPatient, maxTime, nPatients, nWards, proposedInfTimes, currentRecTimes, ptLocation);

            getSporeForceSummary(proposedSporeForceSummary, infectedPatients, proposedSporeI, maxTime, nWards, nPatients, proposedInfTimes, currentParm);
            

            if(debugPt) {
                double currentLLSample = llSample(infectedPatients, currentInfTimes, sampleTimes, currentParm);
                double proposedLLSample = llSample(infectedPatients, proposedInfTimes, sampleTimes, currentParm);
                printf("LL sample, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLSample, proposedLLSample, proposedLLSample-currentLLSample);
                
                double currentLLGenetic = llGenetic(infectedPatients, currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
                double proposedLLGenetic = llGenetic(infectedPatients, proposedInfTimes, sampleTimes, proposedInfSources, proposedInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
                printf("LL genetic, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLGenetic, proposedLLGenetic, proposedLLGenetic-currentLLGenetic);
                
                double currentLLTrans = llTrans(infectedPatients, uninfectedPatients, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporeI, currentSporeForceSummary,
                                                wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, currentParm);
                double proposedLLTrans = llTrans(infectedPatients, uninfectedPatients, proposedInfTimes, proposedInfSourceTypes, proposedInfSources, proposedSporeI, proposedSporeForceSummary, wardLog, inPtDays, ptLocation, proposedWardI, nPatients, nWards, maxTime, currentParm);
                printf("LL transmission, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLTrans, proposedLLTrans, proposedLLTrans-currentLLTrans);
                
                double currentLLRecovery = llRecover(infectedPatients, sampleTimes, currentRecTimes, currentParm);
                double proposedLLRecovery = llRecover(infectedPatients, sampleTimes, currentRecTimes, currentParm);
                printf("LL recovery, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLRecovery, proposedLLRecovery, proposedLLRecovery-currentLLRecovery);
                
                double prior = getPrior(currentParm);
                
                double checkCurrentLL = currentLLSample + currentLLGenetic + currentLLTrans + currentLLRecovery + prior;
                double checkProposedLL = proposedLLSample + proposedLLGenetic + proposedLLTrans + proposedLLRecovery + prior;
                double diff = checkProposedLL - checkCurrentLL;
                printf("check LL current: %0.4f\ncheck LL proposed %0.4f\ndiff: %0.4f\nhastings ratio: %0.4f\n", checkCurrentLL, checkProposedLL, diff, hastingsRatio);
                
                
            }
            
            
            
            //2. Compute the probability of accepting the proposed jump.
            proposedLL = targetDist(infectedPatients, uninfectedPatients, proposedInfTimes, sampleTimes, currentRecTimes, proposedInfSources, proposedInfSourceTypes,
                                    proposedSporeI, proposedSporeForceSummary,
                                    wardLog, inPtDays, ptLocation, proposedWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, currentParm);
            
            probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
            
            
            //3. Accept proposed jump with probability pAccept
            if ( log(runif(0,1)) < probAccept ) {
                if (debugPt) printf("Move accepted (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                //accept the jump and log the change
                chainInfTimes[i] = proposedInfTimes;
                chainInfSources[i] = proposedInfSources;
                chainInfSourceTypes[i] = proposedInfSourceTypes;
                
                //increment accepted counter
                nAccepted[13] += 1;
                
                //update onward transmission log
                onwardTransmission = proposedOnwardTransmission;
                
                //set current values to the proposed values
                currentInfTimes = proposedInfTimes;
                currentInfSources = proposedInfSources;
                currentInfSourceTypes = proposedInfSourceTypes;;
                currentWardI = proposedWardI;
                currentLL = proposedLL;
                currentSporeI = proposedSporeI;
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
            }
            
            
            
        } //end of loop over proposed patietns to update
        
        

        
        
        //UPDATE RECOVERY TIMES
            // Hastings ratio here = 0, as updates are proposed symmetrically
        for (int proposedPatientIndex : patientUpdateList) {
            int proposedPatient = infectedPatients[proposedPatientIndex];
            
            debugPt = false;
            if (proposedPatient==-1) debugPt= true;
            
            
            int proposedRecTime; //store the new proposed recovery time
            
            double hastingsRatio = 0;
            proposedRecTimes = currentRecTimes;
            
            //propose new recovery time - RETURN TO THIS TO REVIEW WHAT HAPPENS WHEN SET TO INFECTED AT t=0
            int deltaRec = 10;
            proposedRecTime = proposeRecoveryTime(proposedPatient, currentRecTimes[proposedPatient],
                                                  onwardTransmission, sampleTimes, maxTime, deltaRec, ptLocation, currentInfTimes);
            
            //proposedRecTime = recoveryTimes[proposedPatient];//temp OVER-RIDE for TESTING
            
            if(proposedRecTime==-1) {
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
            
            //determine proposed infectious individuals on wards under new recovery time
            proposedWardI = getWardI(nPatients, maxTime, nWards, currentInfTimes, proposedRecTimes, wardLog);
            //update sporeI to reflect new recovery time
            proposedSporeI = currentSporeI;
            updateSporeI(proposedSporeI, proposedPatient, maxTime, nPatients, nWards, currentInfTimes, proposedRecTimes, ptLocation);
            getSporeForceSummary(proposedSporeForceSummary, infectedPatients, proposedSporeI, maxTime, nWards, nPatients, currentInfTimes, currentParm);
            
            
            if(debugPt) {
                double currentLLRecovery = llRecover(infectedPatients, sampleTimes, currentRecTimes, currentParm);
                double proposedLLRecovery = llRecover(infectedPatients, sampleTimes, proposedRecTimes, currentParm);
                printf("LL recovery, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLRecovery, proposedLLRecovery, proposedLLRecovery-currentLLRecovery);
            }
            
            
            //2. Compute the probability of accepting the proposed jump.
            proposedLL = targetDist(infectedPatients, uninfectedPatients, currentInfTimes, sampleTimes, proposedRecTimes, currentInfSources, currentInfSourceTypes,
                                    proposedSporeI, proposedSporeForceSummary,
                                    wardLog, inPtDays, ptLocation, proposedWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, currentParm);
            
            probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
            
            
            //3. Accept proposed jump with probability pAccept
            if ( log(runif(0,1)) < probAccept ) {
                if (debugPt) printf("Move accepted (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                //accept the jump and log the change
                chainRecTimes[i] = proposedRecTimes;
                
                //increment accepted counter
                nAccepted[14] += 1;
                
                //set current values to the proposed values
                currentRecTimes = proposedRecTimes;
                currentWardI = proposedWardI;
                currentLL = proposedLL;
                currentSporeI = proposedSporeI;
                currentSporeForceSummary = proposedSporeForceSummary;
                
                //printf("Current LL check %f\n", currentLLCheck);
                //printf("Current LL updated to %f\n", currentLL);
                
            }
            else {
                if (debugPt) printf("Move rejected (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                // reject the proposed jump, stay at current position
                chainRecTimes[i] = currentRecTimes;;
            }
            
            
        } //end of loop over proposed patietns to update


        
        
        // DISRUPTION STEP TO DATA AUGMENTATION - UPDATES INFECTION TIMES TO ALLOW REVERSAL OF TRANSMISSION DIRECTIONS
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
            for (int patient : infectedPatients) {
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
                    for(int victim=0; victim<nPatients; victim++) {
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
                    
                    for(int patient=0;patient<nPatients;patient++) { //iterate over all possible sources of infection for node
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
                    
                    int delta = 4;
                    int proposedInfTime = proposeInfectionTime(node, currentInfTimes[node], proposedOnwardTransmission,
                                                           sampleTimes, maxTime, delta, ptLocation);
                    if(proposedInfTime ==-1) {
                        //keep current infection time
                        proposedInfTime = currentInfTimes[node];
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
                    int delta = 5;
                    int proposedRecoveryTime = proposeRecoveryTime(node, currentRecTimes[node], proposedOnwardTransmission, sampleTimes, maxTime, delta, ptLocation, proposedInfTimes);
                    if(proposedRecoveryTime ==-1) {
                        //keep current infection time
                        proposedRecoveryTime = currentRecTimes[node];
                    }
                    proposedRecTimes[node] = proposedRecoveryTime;
                }
                
                
                // 5. working forwards in time using the new infection times determine the source of each case
                // 6. each change in source generates a component of the hastings ratio as before
                
                // update all as a block
                // for proposedInfTimes, also need to update sporeI based on these - unlike above where only updating one case
                proposedWardI = getWardI(nPatients, maxTime, nWards, proposedInfTimes, proposedRecTimes, wardLog);
                getSporeI(proposedSporeI, infectedPatients, nPatients, maxTime, nWards, proposedInfTimes, proposedRecTimes, ptLocation);
                getSporeForceSummary(proposedSporeForceSummary, infectedPatients, proposedSporeI, maxTime, nWards, nPatients, proposedInfTimes, currentParm);
                
                //currentInfSourceTypes - is used to determine genetic likelihood - use this for all updates current --> proposed updates
                // and converse (proposedInfSourceTypes) for proposed --> current updates
                
                double hastingsRatio = 0.00;
                
                for (int node : nodeSet) {
                    
                    //get new source for each case
                    Src proposedSource = proposeConditionalSource(infectedPatients, node, proposedInfTimes[node], proposedInfTimes,
                                                sampleTimes, proposedRecTimes, wardLog, currentInfSourceTypes, proposedSporeI, ptLocation,
                                                geneticDist, geneticMap, nWards, nPatients, currentParm);
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
                    // using currentInfTimes, currentRecTimes, currentSporeI
                    double pChooseSourceCurrent;
                    SrcList srcProbsCurrent = getSourceProb(infectedPatients, node, currentInfTimes[node], currentInfTimes,
                                sampleTimes, currentRecTimes, wardLog, proposedInfSourceTypes, currentSporeI, ptLocation,
                                geneticDist, geneticMap, nWards, nPatients, currentParm);
                    
                    for (int ii=0; ii<srcProbsCurrent.sourceList.size(); ii++) {
                        //iterate over possible sources until find current source
                        if(srcProbsCurrent.sourceList[ii] ==  currentInfSources[node]) {
                            pChooseSourceCurrent = srcProbsCurrent.sourceProbabilities[ii];
                            break;
                        }
                    }
                    hastingsRatio += log(pChooseSourceCurrent);
                }

                
                
                // 7. calculate a new likelhood, an acceptance ratio, and accept or reject all updates as a block
                
                //determine proposed infectious individuals on wards under new infection time
                
                
                proposedLL = targetDist(infectedPatients, uninfectedPatients, proposedInfTimes, sampleTimes, proposedRecTimes, proposedInfSources, proposedInfSourceTypes,
                                        proposedSporeI, proposedSporeForceSummary,
                                        wardLog, inPtDays, ptLocation, proposedWardI, nPatients, nWards, maxTime, geneticDist, geneticMap, currentParm);
                
                probAccept = min( log(1), proposedLL -  currentLL + hastingsRatio);
                
                if(debugDisrupt) {
                    double currentLLSample = llSample(infectedPatients, currentInfTimes, sampleTimes, currentParm);
                    double proposedLLSample = llSample(infectedPatients, proposedInfTimes, sampleTimes, currentParm);
                    printf("LL sample, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLSample, proposedLLSample, proposedLLSample-currentLLSample);
                    
                    double currentLLGenetic = llGenetic(infectedPatients, currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
                    double proposedLLGenetic = llGenetic(infectedPatients, proposedInfTimes, sampleTimes, proposedInfSources, proposedInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
                    printf("LL genetic, current: %0.3f, proposed: %0.3f, difference = %0.3f\n", currentLLGenetic, proposedLLGenetic, proposedLLGenetic-currentLLGenetic);
                    
                    double currentLLTrans = llTrans(infectedPatients, uninfectedPatients, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporeI, currentSporeForceSummary,
                                                    wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, currentParm);
                    double proposedLLTrans = llTrans(infectedPatients, uninfectedPatients, proposedInfTimes, proposedInfSourceTypes, proposedInfSources,
                                                     proposedSporeI, proposedSporeForceSummary,
                                                     wardLog, inPtDays, ptLocation, proposedWardI, nPatients, nWards, maxTime, currentParm);
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
                    nAccepted[15] += 1;
                    
                    //update onward transmission log
                    onwardTransmission = proposedOnwardTransmission;
                    
                    //set current values to the proposed values
                    currentInfTimes = proposedInfTimes;
                    currentInfSources = proposedInfSources;
                    currentInfSourceTypes = proposedInfSourceTypes;
                    currentRecTimes = proposedRecTimes;
                    currentWardI = proposedWardI;
                    currentLL = proposedLL;
                    currentSporeI = proposedSporeI;
                    currentSporeForceSummary = proposedSporeForceSummary;
                    
                }
                else {
                    if (debugDisrupt) printf("Move rejected (acceptance prob: %f, currentLL: %f, proposedLL: %f)\n\n", probAccept, currentLL, proposedLL);
                    // reject the proposed jump, stay at current position, already saved in main data augmentation step
                }
                
                if (debugDisrupt) {
                    debugPt = saveDebugPt;
                    printf("\n\n");
                }
            }
            
        } //end of disruption section
        
        double currentLLSample = llSample(infectedPatients, currentInfTimes, sampleTimes, currentParm);
        double currentLLGenetic = llGenetic(infectedPatients, currentInfTimes, sampleTimes, currentInfSources, currentInfSourceTypes, geneticDist, geneticMap, nPatients, currentParm);
        double currentLLTrans = llTrans(infectedPatients, uninfectedPatients, currentInfTimes, currentInfSourceTypes, currentInfSources, currentSporeI,
                                        currentSporeForceSummary, wardLog, inPtDays, ptLocation, currentWardI, nPatients, nWards, maxTime, currentParm);
        double currentLLRecovery = llRecover(infectedPatients, sampleTimes, currentRecTimes, currentParm);

        //log posterior
        chain[i][13] = currentLL;
        //log trans likelihood
        chain[i][14] = currentLLTrans;
        //log genetic likelihood
        chain[i][15] = currentLLGenetic;
        // log sampling likelihood
        chain[i][16] = currentLLSample;
        //log recovery likelihood
        chain[i][17] = currentLLRecovery;
        
        
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
    vector <int> infTimes; //if no infection then set to -1
    vector <int> infSources; //if no infection or background/starts infected set to -1
    vector <int> infSourceTypes; //hospital background = 0, ward =1, hospital = 2, comm background = 3, started infection = 4, no infection = -1
    vector <int> sampleTimes; // if not infection set to -1
    vector <int> recoverTimes; // if not infection set to -1
    
    int nPatients;
    string filePath = path + "/input/patientLog.csv";
    importPatientLog(filePath, infTimes, infSources, infSourceTypes, sampleTimes, recoverTimes, nPatients);
    
    filePath = path + "/input/wardLog.csv";
    vector<vector<vector<int>>> wardLog;
    int nWards, maxTime;
    importWardLog(filePath, wardLog, maxTime, nWards, nPatients, sampleTimes);  // wardLog[time][ward] = {patients...}
    vector<int> infectedPatients = getInfectedPatients(sampleTimes, nPatients);
    
    //read in genetic data
    string filePathGenetic = path + "/input/simDistances_snps.txt";
    vector<vector<double>> geneticDist; //pairwise genetic distances
    unordered_map<int,int> geneticMap; //maps the row and column of the 2d vector back to the patient ids
    importGeneticData(filePathGenetic, geneticDist, geneticMap);

    // **** RUN MCMC **** //
    // MCMC setup
    //clock_t startMCMC,finishMCMC;
    //startMCMC = clock();
    struct timeval startMCMC, endMCMC;
    
    gettimeofday(&startMCMC, NULL);
    double deltaMCMC;

    
    
    //parameters: "beta0", "beta1", "beta2", "sampleSize", "sampleMu", "directNe", "introNe", "mu", "startInf", "betaComm", "sporeProb", "recoverSize", "recMu",
    Parm startParm;
    startParm.betaBgroundHosp = 0.005;
    startParm.betaWard = 0.005;
    startParm.betaHosp = 0.005;
    startParm.sampleSize = 3;
    startParm.sampleMu = 5;
    startParm.directNe =  1;
    startParm.introNe = 500;
    startParm.mu = 2/365.25;
    startParm.probStartInfLogit = logit(0.01);
    startParm.betaComm = 0.005;
    startParm.sporeProbLogit = 0.5;
    startParm.recSize = 3;
    startParm.recMu = 30;
    
    vector<double> startSigma = {0.001, 0.001, 0.001, 0.2, 0.2, 0.2, 100, 0.01, 0.0005, 0.1, 0.1, 0.1, 0.1};
    
    
    //2d vector for results = columns: "beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu", "startInf", "betaComm"; rows = each iteration
    //2d vectors for infection times, infection sources, and source types, columns for each patient, row = each iteration
    vector<Parm> chain;
    vector<vector<int>> chainInfTimes, chainInfSources, chainInfSourceTypes, chainRecTimes;
    
    //set up sizes of chain for parameters and various times as 2d vectors
    chain.resize(steps);
    chainInfTimes.resize(steps); chainInfSources.resize(steps);
    chainInfSourceTypes.resize(steps); chainRecTimes.resize(steps);
    
    for (int i = 0; i<steps; i++) {
        chainInfTimes[i].resize(nPatients);
        chainInfSources[i].resize(nPatients);
        chainInfSourceTypes[i].resize(nPatients);
        chainRecTimes[i].resize(nPatients);

    }
    
    vector<vector<vector<int>>> inPtDays = getInPtDays(nPatients, maxTime, nWards, wardLog);
    vector<vector<int>> ptLocation = getPtLocation(nPatients, maxTime, nWards, inPtDays);
    vector<vector<int>> wardI = getWardI(nPatients, maxTime, nWards, infTimes, recoverTimes, wardLog);
    
    
    //vector<vector<vector<int>>> sporeI = getSporeI(nPatients, maxTime, nWards, infTimes, recoverTimes, ptLocation);
    //vector<vector<vector<double>>> sporeForce = getSporeForce(sporeI, maxTime, nWards, nPatients, startParm);
    //vector<vector<double>> sporeForceSummary = getSporeForceSummary(sporeForce, maxTime, nWards, nPatients);
    
    
    //runTest(infTimes, sampleTimes, recoverTimes, startParm,
    //        infSources, infSourceTypes,
    //        sporeI, sporeForce, sporeForceSummary, geneticDist, geneticMap, wardLog,
    //        inPtDays, ptLocation, wardI, nPatients, nWards, maxTime);
    
    
    
    doMCMC(chain, chainInfTimes, chainRecTimes, chainInfSources, chainInfSourceTypes, steps, startParm, startSigma,
           sampleTimes, wardLog, nPatients, nWards, maxTime, geneticDist, geneticMap,
           infTimes, infSources, infSourceTypes, recoverTimes);
    
    
    //export the chain to a file
    const string mcmcLog = path + "/inference"; //path for log files
    exportChain(chain, chainInfTimes, chainInfSources, chainInfSourceTypes, chainRecTimes, steps, mcmcLog);
    
    
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
