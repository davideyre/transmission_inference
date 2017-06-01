//
//  proposals.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#include "proposals.hpp"



//function to propose a new infection time - based on uniform dis  (note will update recovery time at same iteration, and so not constrained by recovery)
int proposeInfectionTime(int proposedPatient, int currentInfTime,
                         vector<vector<trans>> &onwardTransmission,
                         vector<int> &sampleTimes,
                         int maxTime,
                         double delta,
                         vector<vector<int>> &ptLocation) {
    
    //over-ride for patients already positive at start
    if(sampleTimes[proposedPatient]==0) {
        return 0;
    }
    
    //acceptable interval for infection
    // from t=1
    // to max Time or time step before first onward transmission
    
    //nb updates to infection time can only alter onward transmission from spore that is set by ward discharge, as recovery time is not changed

    
    // initially only restrict for onward transmission from routes other than spore
    int onwardTransmissionMin = maxTime;
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.t<onwardTransmissionMin & transmission.srcType!=5) {
            //set new minimum onward transmission, excluding transmissions from spore
            onwardTransmissionMin = transmission.t;
        }
    }
    int endTime = min({sampleTimes[proposedPatient], onwardTransmissionMin-1});
    
    //sample new delta
    //int proposedInfTime = currentInfTime + floor(runif(-delta, delta+1));
    int proposedInfTime = currentInfTime + round(rnorm(0, delta)); //sample from normal distribution centred on zero, with sd delta
    
    if (proposedInfTime<1 | proposedInfTime>endTime) {
        //impossible move - reject
        proposedInfTime = std::numeric_limits<int>::min();;
        
        return proposedInfTime;
    }
    
    
    //now check for onward transmission from spore arising from ward discharge whlie infectious - moving earlier will not result in loss of opportunity to leave spore, therefore...
    //only applies if moving the infection time later
    if(proposedInfTime>currentInfTime) {
        //check for spore transmission from this patient
        for (trans transmission : onwardTransmission[proposedPatient]) {
            if(transmission.srcType==5) {
                //check if any of these spore transmissions depend on discharge between the current and proposed infection times
                int dependDischarge = 0;
                int currentLocation = ptLocation[proposedPatient][currentInfTime+1]; //location at time first infectious
                
                
                for(int t=currentInfTime+2; t<=proposedInfTime+1; t++) {
                    if(ptLocation[proposedPatient][t] != currentLocation &  currentLocation!=-1) {
                        //a discharge has occured
                        int dischargeDate = t-1;
                        int dischargeLocation = ptLocation[proposedPatient][dischargeDate];
                        
                        if(transmission.ward==dischargeLocation & dischargeDate<transmission.t) {
                            //discharge from the ward of transmission, before the spore transmission between current and proposed times becomes infectious
                            dependDischarge=1;
                            
                            //check no subsequent discharge from this ward before transmission.t
                            //first check no re-admission
                            for(int tt=t; tt<=proposedInfTime+1; tt++) {
                                if(ptLocation[proposedPatient][tt] == transmission.ward) {
                                    //readmission has occured, now check for discharge
                                    for (int ttt=tt; ttt<=proposedInfTime+1; ttt++) {
                                        if(ptLocation[proposedPatient][ttt]!= transmission.ward) {
                                            ///discharge has occured
                                            dependDischarge = 0;
                                            break;
                                        }
                                    }
                                }
                            }
                            
                            if(dependDischarge==1) {
                                break;
                            }

                        }
                    }
                    currentLocation = ptLocation[proposedPatient][t];
                }
                
                if(dependDischarge==1) {
                    
                    //onward transmission depends on being discharged while infectious
                    //printf("Proposed infection time does not allow for a transmission requiring a discharge while infectious, move rejected\n");
                    proposedInfTime = std::numeric_limits<int>::min(); //impossible move, reject and exit
                    return(proposedInfTime);
                }

            } //end of spore transmissions
        } //end of transmissions loop
    }
    
    //if not problem with infection time found -
    return proposedInfTime;
}



//function to propose a new recovery time
  // this is constrained by {last onward transmission and sampling} to {first spore transmission}
int proposeRecoveryTime(int proposedPatient, int currentRecTime,
                        vector<vector<trans>> &onwardTransmission, vector<int> &sampleTimes, int maxTime, double delta,
                        vector<vector<int>> &ptLocation, vector<int> &infTimes) {
    
    //propose a new recovery time
    //int proposedRecTime = currentRecTime + floor(runif(-delta, delta+1));
    int proposedRecTime = currentRecTime + round(rnorm(0, delta)); //sample from normal distribution centred on zero, with sd delta
    
    //define minimum step step can recover in,
    //  i.e. the time step after the latest of last onward transmission from direct infection (rather than spores) or sampling
    vector<int> transmissionTimes; //vector of onward transmission times from direct infection
    int minRecoveryTime; //minimum time step can recover in
    
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.srcType != SrcType::SPORE){
            transmissionTimes.push_back(transmission.t); //add to vector of transmission times provided not spore transmission
        }
    }
    
    if(transmissionTimes.size()>0)  {
        int lastTransmission = *max_element(begin(transmissionTimes), end(transmissionTimes));
        minRecoveryTime = max({lastTransmission+1,sampleTimes[proposedPatient]});  //i.e. not infectious in time step recover in
    }
    else {
        minRecoveryTime = sampleTimes[proposedPatient];
    }
    
    if(proposedRecTime< minRecoveryTime) {
        proposedRecTime = std::numeric_limits<int>::min();;; //impossible move, reject and exit
        return(proposedRecTime);
    }
    
    
    //determine if new recovery time interfers with any existing spore transmission
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.srcType == SrcType::SPORE){
            //iterate over spore transmissions from this patient
            int sporeError = 9; //default spore error
            
            //check if the spore could be created by patient recovery under the new recovery time
            if(proposedRecTime<=maxTime) {
                if(transmission.ward==ptLocation[proposedPatient][proposedRecTime]) {
                    //patient recovered on the same ward as the spore transmission, check timings
                    //patient must recover in or before the time step of spore transmission
                    if(proposedRecTime>transmission.t) {
                        //check the proposed recovery time does not occur after a spore transmission
                        sporeError = 1;
                    }
                    else {
                        //no error
                        sporeError = 0;
                    }
                }
            }
            
            if(sporeError!=0) {
                //need to check if the spore could have been created by a ward discharge between infection time +1 and the new recovery time
                int dischargeFound = 0;
                
                //search for ward discharges
                int currentLocation = ptLocation[proposedPatient][infTimes[proposedPatient]+1];
                for(int t=infTimes[proposedPatient]+2; t<=proposedRecTime+1; t++) {
                    if(ptLocation[proposedPatient][t] != currentLocation &  currentLocation!=-1) {
                        //a discharge has occured
                        int dischargeDate = t-1;
                        int dischargeLocation = ptLocation[proposedPatient][dischargeDate];
                        
                        if(transmission.ward==dischargeLocation & dischargeDate<transmission.t) {
                            //discharge from the ward of transmission, before the spore transmission
                            //no error
                            sporeError = 0;
                            dischargeFound = 1;
                            break;
                        }
                    }
                    currentLocation = ptLocation[proposedPatient][t];
                }
                
                if(dischargeFound==0) {
                    //no discharge to account for spore transmission found
                    sporeError = 2;
                }
            } //end of spore transmissions dependent on discharges
            
            if(sporeError!=0) {
                /*
                if(sporeError==1) {
                    printf("Proposed recovery time is after a spore transmisssion, move rejected\n");
                }
                else if(sporeError==2) {
                    printf("Proposed recovery time does not allow for a transmission requiring a discharge, move rejected\n");
                }
                else {
                    printf("ERROR - HALT\n\n");
                }
                */
                proposedRecTime = std::numeric_limits<int>::min(); //impossible move, reject and exit
                return(proposedRecTime);
            }
            
            
            
        } //end of spore transmissions
    } //end of transmissions loop
    
    
    //if no problems with proposed recovery time identified, allow this to be returned
    return proposedRecTime;
}







/*
//function for independence sampler for infection time, based on current value of episilon
int proposeInfectionTimeConditional(int proposedPatient,
                                    vector<vector<trans>> &onwardTransmission,
                                    vector<int> &sampleTimes,
                                    int maxTime,
                                    double epsilon) {
    
    if(sampleTimes[proposedPatient]==0) {
        return 0;
    }
    
    //define interval to move in, can be from t=1 onwards until time interval before first onward transmission or sampling time
    int onwardTransmissionMin = maxTime;
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.t<onwardTransmissionMin & transmission.srcType!=5) {
            onwardTransmissionMin = transmission.t; //update min onward transmisson, excluding spore transmissions
        }
    }
    int endTime = min({sampleTimes[proposedPatient], onwardTransmissionMin-1});
    
    //sample new infection time conditional on value of epsilon
    int proposedInfTime = sampleTimes[proposedPatient] - rpois(epsilon);
    
    if (proposedInfTime<1 | proposedInfTime>endTime) {
        //impossible move - reject
        proposedInfTime = -1;
    }
    return proposedInfTime;
    
}
*/

//function for initiating infection times (at start or after move from being infected at t=0
int proposeInfectionTimeInitial(int proposedPatient,
                                vector<vector<trans>> &onwardTransmission,
                                vector<int> &sampleTimes,
                                int maxTime)  {
    
    //define interval to move in, can be from t=1 onwards until time interval before first onward transmission or sampling time
    int onwardTransmissionMin = maxTime;
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.t<onwardTransmissionMin & transmission.srcType != SrcType::SPORE) {
            onwardTransmissionMin = transmission.t; // exclude spore transmission
        }
    }
    int endTime = min({sampleTimes[proposedPatient], onwardTransmissionMin-1});
    
    //propose infection time between t=0 and point must have been infectious from (sampling or onward transmission)
    int proposedInfTime = floor(runif(0, endTime+1));
    
    return proposedInfTime;
    
}


//get list of sources and their probabilities for a given patient and infection time, and parameter set
SrcList getSourceProb(vector<int> &infectedPatients, int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                             vector<vector<vector<int>>> &wardLog, vector<int> infSourceType,
                             vector<vector<vector<int>>> &sporeI,
                             vector<vector<int>> &ptLocation,
                             vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap,
                             int nWards,
                             int nPatients, Parm &parm) {
    SrcList output; //struct of type SrcList to store output
    
    vector<int> sourceList; //vector of sources
    vector<int> sourceTypeList; //vector of source types
    vector<double> sourceLikelihood; //vector of likelihood values for each source
    
   
    //potential ward sources - wardLog[time][ward], ptLocation[patient][time]
    int ward = ptLocation[proposedPatient][proposedInfTime];
    
    if (ward==-1) {
        //patient is an outpatient, only one possible source - community background
        output.sourceList = {-1};
        output.sourceTypeList = {3};
        output.sourceProbabilities = {1};
    }
    else {
        //add background source to lists
        sourceList.push_back(-1);
        sourceTypeList.push_back(SrcType::BGROUND_HOSP);
        
        //potential ward sources
        for (int srcPt : wardLog[proposedInfTime][ward]) {
            if (infTimes[srcPt]!=-1 & infTimes[srcPt]<proposedInfTime & recoveryTimes[srcPt]>proposedInfTime & srcPt!=proposedPatient) {
                sourceList.push_back(srcPt);
                sourceTypeList.push_back(SrcType::WARD);
                
            }
        }
        
        //potential hospital sources
        for (int nonWard=0; nonWard<nWards; nonWard++) {
            if (nonWard!=ward) {
                for (int srcPt : wardLog[proposedInfTime][nonWard]) {
                    if (infTimes[srcPt]!=-1 & infTimes[srcPt]<proposedInfTime & recoveryTimes[srcPt]>proposedInfTime & srcPt!=proposedPatient) {
                        sourceList.push_back(srcPt);
                        sourceTypeList.push_back(SrcType::HOSP);
                        
                    }
                }
            }
        }
        
        //potential spore sources
        for (int srcPt = 0; srcPt<nPatients; srcPt ++) {
            if(sporeI[proposedInfTime][ward][srcPt]>0 & srcPt!=proposedPatient) {
                sourceList.push_back(srcPt);
                sourceTypeList.push_back(SrcType::SPORE);
                
            }
        }
        
        
        if(sourceList.size()==1) {
            //if no ward or hospital sources assign to background
            output.sourceList = {-1};
            output.sourceTypeList = {SrcType::BGROUND_COMM};
            output.sourceProbabilities = {1.0};
        }
        else {
            
            //determine the probability of each source
            vector<double> sourceLogLikelihood; //vector of log likelihood values for each source
            
            //parm: "beta0", "beta1", "beta2", "epsilon", "directNe", "introNe", "mu"
            
            
            //get the total spore force at this time-point
            /*
            double sporeIt = 0.00;
            for (int sporePt =0; sporePt<nPatients; sporePt++) {
                int sporeDuration = sporeI[proposedInfTime][ward][sporePt];
                if(sporeDuration>0) {
                    sporeIt += pow((1-parm.sporeProbLogit), sporeDuration);
                }
            }
            */
            
            //calculate the background likelihood
            double backgroundLogLikelihood = log(parm.betaBgroundHosp) + llGeneticSingle(infectedPatients, sampleTimes, proposedPatient, -1, infSourceType, geneticDist, geneticMap, nPatients, parm);
            
            sourceLogLikelihood.resize(sourceList.size());
            sourceLogLikelihood[0] = backgroundLogLikelihood;
            
            //calculate the other likelihoods
            for (int srcIndex = 1; srcIndex<sourceList.size(); srcIndex++) {
                if( sourceTypeList[srcIndex] == SrcType::WARD) {
                    //infection from same ward
                    sourceLogLikelihood[srcIndex] = log(parm.betaWard) + llGeneticSingle(infectedPatients, sampleTimes, proposedPatient, sourceList[srcIndex], infSourceType, geneticDist, geneticMap, nPatients, parm);
                }
                if( sourceTypeList[srcIndex] == SrcType::HOSP) {
                    //infection from hospital
                    sourceLogLikelihood[srcIndex] = log(parm.betaHosp) + llGeneticSingle(infectedPatients, sampleTimes, proposedPatient, sourceList[srcIndex], infSourceType, geneticDist, geneticMap, nPatients, parm);
                }
                if( sourceTypeList[srcIndex] == SrcType::SPORE) {
                    //infection from spore
                    
                    //get the relative infectious level of this spore
                    int specificSporeDuration = sporeI[proposedInfTime][ward][sourceList[srcIndex]];
                    double specificSporeLevel = pow((1-getSporeP(parm)), specificSporeDuration);
                    
                    sourceLogLikelihood[srcIndex] = log(parm.betaWard * specificSporeLevel) + llGeneticSingle(infectedPatients, sampleTimes, proposedPatient, sourceList[srcIndex], infSourceType, geneticDist, geneticMap, nPatients, parm);
                    
                }
                
            }
            
            //calculate the normalised likelihood of each source
            vector<double> sourceProbability = normaliseLL(sourceLogLikelihood);
            
            output.sourceList = sourceList;
            output.sourceTypeList = sourceTypeList;
            output.sourceProbabilities = sourceProbability;
            }
    }
    

    
    return output;
}


//function to determine proposed source of infection conditionally
Src proposeConditionalSource(vector<int> &infectedPatients, int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                             vector<vector<vector<int>>> &wardLog, vector<int> infSourceType, vector<vector<vector<int>>> &sporeI,
                             vector<vector<int>> &ptLocation,
                             vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap,
                             int nWards, int nPatients, Parm &parm) {
    Src output; //struct of type Src to store output
    
    if(sampleTimes[proposedPatient]==0) {
        output = {-1, SrcType::START_POS, 1};
    }
    
    else {
        
        SrcList sourceList = getSourceProb(infectedPatients, proposedPatient, proposedInfTime, infTimes, sampleTimes, recoveryTimes,
                              wardLog, infSourceType, sporeI, ptLocation, geneticDist, geneticMap, nWards, nPatients, parm);
        
        if(sourceList.sourceList.size()==1) {
            if(sourceList.sourceTypeList[0] == 0) {
                //if no ward or hospital sources assign to background - hospital
                output = {-1, SrcType::BGROUND_HOSP, 1};
            }
            else {
                //in community, i.e. only one source, assign to community
                output = {-1, SrcType::BGROUND_COMM, 1};
            }
            
        }
        else {
            //chose a source based on probability
            int i = sampleProbVector(sourceList.sourceProbabilities);
            
            //set ouput
            output = {sourceList.sourceList[i], sourceList.sourceTypeList[i], sourceList.sourceProbabilities[i]};
        }
    }
    
    return output;
}

