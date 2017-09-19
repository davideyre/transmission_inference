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
    
    //acceptable interval for infection
    // from t = -inf, i.e. can be infected before t=0 when have data
    // to max Time or time step before first onward transmission
    
    //NB updates to infection time can only alter onward transmission from spore
    // that is set by ward discharge, as recovery time is not changed

    
    // initially only restrict for onward transmission from routes other than spore
    int onwardTransmissionMin = maxTime;
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.t<onwardTransmissionMin & transmission.srcType!= SrcType::SPORE) {
            //set new minimum onward transmission, excluding transmissions from spore
            onwardTransmissionMin = transmission.t;
        }
    }
    int endTime = min({sampleTimes[proposedPatient], onwardTransmissionMin-1});
        //i.e. must have been infected in interval before onward transmission
    
    //sample new delta
    //int proposedInfTime = currentInfTime + floor(runif(-delta, delta+1));
    int proposedInfTime = currentInfTime + round(rnorm(0, delta)); //sample from normal distribution centred on zero, with sd delta
    
    if (proposedInfTime>endTime) {
        //impossible move - reject
        proposedInfTime = std::numeric_limits<int>::min();;
        return proposedInfTime;
    }
    
    
    //now check for onward transmission from spore arising from ward discharge whlie infectious - moving earlier will not result in loss of opportunity to leave spore, therefore...
    //only applies if moving the infection time later
    if(proposedInfTime>currentInfTime) {
        //check for spore transmission from this patient
        for (trans transmission : onwardTransmission[proposedPatient]) {
            if(transmission.srcType == SrcType::SPORE) {
                
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
    
    //if no problem with infection time found then...
    return proposedInfTime;
}




//function for initiating infection times (at start of inference) - all initialised to be at or after start of data, i.e. t=0
int proposeInfectionTimeInitial(int proposedPatient, vector<int> &sampleTimes, Parm parm)  {
    
    //as all initial infections are set to be from background, don't need to check for onward transmission
    
    //propose infection time sampling from negative binomial distribution
    double sampleProb = parm.sampleSize / (parm.sampleSize + parm.sampleMu);
    int proposedInfTime = sampleTimes[proposedPatient] - rnbinom(parm.sampleSize, sampleProb);
    return proposedInfTime;
}




//function to propose a new recovery time
  // this is constrained by {last onward transmission and sampling} to {first spore transmission}
int proposeRecoveryTime(int proposedPatient, int currentRecTime,
                        vector<vector<trans>> &onwardTransmission, vector<int> &sampleTimes, int maxTime, int minTime, double delta,
                        vector<vector<int>> &ptLocation, vector<int> &infTimes) {
    
    //propose a new recovery time
    //int proposedRecTime = currentRecTime + floor(runif(-delta, delta+1));
    int proposedRecTime = currentRecTime + round(rnorm(0, delta)); //sample from normal distribution centred on zero, with sd delta
    
    //define minimum step step can recover in,
    //  i.e. the time step after the latest of last onward transmission from direct infection (rather than spores) or sampling
    vector<int> transmissionTimes; //vector of onward transmission times from direct infection
    int minRecoveryTime; //minimum time step can recover in
    
    for (trans transmission : onwardTransmission[proposedPatient]) {
        if(transmission.srcType != SrcType::SPORE & transmission.victim != -1){
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
                int startTime = max({minTime, infTimes[proposedPatient]+1}); //find search window - set to start at minTime for cases infected before minTime
                int currentLocation = ptLocation[proposedPatient][startTime];
                for(int t=startTime+1; t<=proposedRecTime+1; t++) {
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



//get list of sources and their probabilities for a given patient and infection time, and parameter set
SrcList getSourceProb(int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                             vector<vector<vector<int>>> &wardLogInf, vector<int> infSourceType,
                             vector<vector<vector<SporeEvent>>> &sporePatientI,
                             vector<vector<int>> &ptLocation,
                             vector<vector<double>> &geneticDist,
                             int nWards,
                             int nInfPatients, vector<vector<int>> &hospitalWards, Parm &parm) {
    SrcList output; //struct of type SrcList to store output
    
    vector<int> sourceList; //vector of sources
    vector<int> sourceTypeList; //vector of source types
    vector<double> sourceLikelihood; //vector of likelihood values for each source
    
    if (proposedInfTime<0) {
        //patient is infected before start of data - record this
        output.sourceList = {-1};
        output.sourceTypeList = {SrcType::START_POS};
        output.sourceProbabilities = {1};

    }
    else {
        //patient infected on or after t=0
        
        //potential ward sources - wardLog[time][ward], ptLocation[patient][time]
        int ward = ptLocation[proposedPatient][proposedInfTime];
        
        if (ward==-1) {
            //patient is an outpatient, only one possible source - community background
            output.sourceList = {-1};
            output.sourceTypeList = {SrcType::BGROUND_COMM};
            output.sourceProbabilities = {1};
        }
        else {
            //add hospital background source to lists
            sourceList.push_back(-1);
            sourceTypeList.push_back(SrcType::BGROUND_HOSP);
            
            //potential ward sources
            for (int srcPt : wardLogInf[proposedInfTime][ward]) {
                if (infTimes[srcPt]!=-1 & infTimes[srcPt]<proposedInfTime & recoveryTimes[srcPt]>proposedInfTime & srcPt!=proposedPatient) {
                    sourceList.push_back(srcPt);
                    sourceTypeList.push_back(SrcType::WARD);
                    
                }
            }
            
            //potential hospital sources
            for (int nonWard : hospitalWards[ward]) {
                for (int srcPt : wardLogInf[proposedInfTime][nonWard]) {
                    if (infTimes[srcPt]!=-1 & infTimes[srcPt]<proposedInfTime & recoveryTimes[srcPt]>proposedInfTime & srcPt!=proposedPatient) {
                        sourceList.push_back(srcPt);
                        sourceTypeList.push_back(SrcType::HOSP);
                        
                    }
                }
            }
            
            //potential spore sources
            for (int srcPt = 0; srcPt<nInfPatients; srcPt ++) {
                if(!sporePatientI[ward][srcPt].empty()) {
                    for (SporeEvent sporeEvent: sporePatientI[ward][srcPt]) {
                        if(sporeEvent.start<=proposedInfTime & sporeEvent.end>=proposedInfTime & srcPt!=proposedPatient) {
                            sourceList.push_back(srcPt);
                            sourceTypeList.push_back(SrcType::SPORE);
                            break;
                        }
                    }
                }
            }
            
            
            if(sourceList.size()==1) {
                //if no ward or hospital sources assign to background
                output.sourceList = {-1};
                output.sourceTypeList = {SrcType::BGROUND_HOSP};
                output.sourceProbabilities = {1.0};
            }
            else {
                
                //determine the probability of each source
                vector<double> sourceLogLikelihood; //vector of log likelihood values for each source
                
                //calculate the background likelihood
                double backgroundLogLikelihood = log(parm.betaBgroundHosp) + llGeneticSingle(sampleTimes, proposedPatient, -1,
                                                                                             infSourceType, geneticDist, nInfPatients, parm);
                
                sourceLogLikelihood.resize(sourceList.size());
                sourceLogLikelihood[0] = backgroundLogLikelihood;
                
                //calculate the other likelihoods
                for (int srcIndex = 1; srcIndex<sourceList.size(); srcIndex++) {
                    if( sourceTypeList[srcIndex] == SrcType::WARD) {
                        //infection from same ward
                        sourceLogLikelihood[srcIndex] = log(parm.betaWard) + llGeneticSingle(sampleTimes, proposedPatient, sourceList[srcIndex],
                                                                                             infSourceType, geneticDist, nInfPatients, parm);
                        
                    }
                    if( sourceTypeList[srcIndex] == SrcType::HOSP) {
                        //infection from hospital
                        sourceLogLikelihood[srcIndex] = log(parm.betaHosp) + llGeneticSingle(sampleTimes, proposedPatient, sourceList[srcIndex],
                                                                                             infSourceType, geneticDist, nInfPatients, parm);
                    }
                    if( sourceTypeList[srcIndex] == SrcType::SPORE) {
                        //infection from spore
                        
                        //get the relative infectious level of this spore
                        double specificSporeLevel = 0;
                        for(SporeEvent sporeEvent: sporePatientI[ward][sourceList[srcIndex]]) {
                            if (sporeEvent.start <= proposedInfTime & sporeEvent.end >= proposedInfTime) {
                                int specificSporeDuration = proposedInfTime - sporeEvent.start + 1;
                                specificSporeLevel += pow((1-getSporeP(parm)), specificSporeDuration);
                            }
                        }

                        sourceLogLikelihood[srcIndex] = log(parm.betaWard * getSporeMultiplier(parm) * specificSporeLevel) + llGeneticSingle(sampleTimes, proposedPatient, sourceList[srcIndex],
                                                                                                                  infSourceType, geneticDist, nInfPatients, parm);
                    }
                    
                }
                
                //calculate the normalised likelihood of each source
                vector<double> sourceProbability = normaliseLL(sourceLogLikelihood);
                
                output.sourceList = sourceList;
                output.sourceTypeList = sourceTypeList;
                output.sourceProbabilities = sourceProbability;
            }
            
        } //end of ward==-1 condition
    } //end of t<0 condition
    
    return output;
}


//function to determine proposed source of infection conditionally
Src proposeConditionalSource(int proposedPatient, int proposedInfTime, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoveryTimes,
                             vector<vector<vector<int>>> &wardLog, vector<int> infSourceType, vector<vector<vector<SporeEvent>>> &sporePatientI,
                             vector<vector<int>> &ptLocation,
                             vector<vector<double>> &geneticDist, int nWards, int nInfPatients, vector<vector<int>> &hospitalWards, Parm &parm) {
    Src output; //struct of type Src to store output
    
    if(proposedInfTime<0) {
        output = {-1, SrcType::START_POS, 1}; //{source, source type, probability}
    }
    
    else {
        
        SrcList sourceList = getSourceProb(proposedPatient, proposedInfTime, infTimes, sampleTimes, recoveryTimes,
                              wardLog, infSourceType, sporePatientI, ptLocation, geneticDist, nWards, nInfPatients, hospitalWards, parm);
        
        if(sourceList.sourceList.size()==1) {
            if(sourceList.sourceTypeList[0] == SrcType::BGROUND_HOSP) {
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

