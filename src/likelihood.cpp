//
//  likelihood.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#include "likelihood.hpp"


double llTransAvoid(vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &infSourceType, vector<int> &infSources,
               vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
               vector<vector<vector<int>>> &wardLog,
               vector<vector<vector<int>>> &inPtDays,
               vector<vector<int>> &ptLocation,
               vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, Parm &parm) {
    
    double ll = 0.00;
    
    //ll for uninfectedPatients
    int nUninfected = (int)uninfectedPatients.size();
    int wardInfDays = 0;
    int hospInfDays = 0;
    double sporeInfDays = 0.00;
    int admissionDuration = 0;
    
    // Pr(not infected at start) //
    double probNotStartInf = 1 - getStartInfP(parm);
    ll += nUninfected*log(probNotStartInf);
    
    //get counts of number of each type of day uninfected patients exposed to
    
    for(int patient : uninfectedPatients) {
        //loop through all days an inpatient - inPtDays[patient][ward] = {times...}
        for(int ward=0; ward < nWards; ward++) {
            if (inPtDays[patient][ward].size()>0) { // if inpatient on that ward
                for (int t: inPtDays[patient][ward]) {
                    //log of admission duration for background exposure
                    admissionDuration ++;
                    
                    //ward infected person-days exposed to - wardI[t][ward] = nInf
                    wardInfDays += wardI[t][ward];
                    
                    //get infectious numbers from other wards
                    for(int nonWard=0; nonWard < nWards; nonWard++) {
                        if(nonWard!=ward) {
                            hospInfDays += wardI[t][nonWard];
                        }
                    }
                    
                    //get number of spore day equivalents - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
                    //spores decay according to geometric distribution
                    sporeInfDays += sporeForceSummary[t][ward];
                }
            }
        }
    }
    
    int nonAdmissionDuration = (nUninfected*maxTime) - admissionDuration;
    
    double logProbAvoidInf = -((parm.betaBgroundHosp * admissionDuration) //hospital background
                               + (parm.betaWard * wardInfDays) //ward infection force
                               + (parm.betaHosp * hospInfDays) //across hospital infection force
                               + (parm.betaComm * nonAdmissionDuration) //community background
                               + (parm.betaWard * sporeInfDays ) //ward-based spore force
                               );
    ll += logProbAvoidInf;
    return (ll);
}



//log likelihood contribution from transmission model - p(I | parm)

double llTrans(vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &infSourceType, vector<int> &infSources,
               vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
               vector<vector<vector<int>>> &wardLog,
               vector<vector<vector<int>>> &inPtDays,
               vector<vector<int>> &ptLocation,
               vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, Parm &parm) {
     
    
    double ll = 0.00;
    
    //ll for uninfectedPatients
    int nUninfected = (int)uninfectedPatients.size();
    int wardInfDays = 0;
    int hospInfDays = 0;
    double sporeInfDays = 0.00;
    int admissionDuration = 0;
    
    // Pr(not infected at start) //
    double probNotStartInf = 1 - getStartInfP(parm);
    ll += nUninfected*log(probNotStartInf);
    
    //get counts of number of each type of day uninfected patients exposed to
    
    for(int patient : uninfectedPatients) {
        //loop through all days an inpatient - inPtDays[patient][ward] = {times...}
        for(int ward=0; ward < nWards; ward++) {
            if (inPtDays[patient][ward].size()>0) { // if inpatient on that ward
                for (int t: inPtDays[patient][ward]) {
                    //log of admission duration for background exposure
                    admissionDuration ++;
                    
                    //ward infected person-days exposed to - wardI[t][ward] = nInf
                    wardInfDays += wardI[t][ward];
                    
                    //get infectious numbers from other wards
                    for(int nonWard=0; nonWard < nWards; nonWard++) {
                        if(nonWard!=ward) {
                            hospInfDays += wardI[t][nonWard];
                        }
                    }
                    
                    //get number of spore day equivalents - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
                    //spores decay according to geometric distribution
                    sporeInfDays += sporeForceSummary[t][ward];
                }
            }
        }
    }

    int nonAdmissionDuration = (nUninfected*maxTime) - admissionDuration;

    double logProbAvoidInf = -((parm.betaBgroundHosp * admissionDuration) //hospital background
                                         + (parm.betaWard * wardInfDays) //ward infection force
                                         + (parm.betaHosp * hospInfDays) //across hospital infection force
                                         + (parm.betaComm * nonAdmissionDuration) //community background
                                         + (parm.betaWard * sporeInfDays ) //ward-based spore force
                                         );
    ll += logProbAvoidInf;

    
    
    
    //loop over infected patients
    
    for(int patient : infectedPatients) {
        //infected pateints
        if (infTimes[patient] < 0) {
            //patients infected before start
            ll += log(getStartInfP(parm)); //probabilty that start infected
        }
        
        else {
            //patients with infection after t=0
            
            // Pr(not infected at start) * Pr(avoid infection up until t-1) * Pr(infecetd at t) * Pr(infected by specific source)
            
            // Pr(not infected at start) //
            double probNotStartInf = 1 - getStartInfP(parm);
            
            // Pr(avoid infection up until t-1)
            // as written can account for patients infected on day of admission, i.e. probAvoidInf==1
            
            //get counts of number of infectious patient days exposed to before infected
            int wardInfDays = 0;
            int hospInfDays = 0;
            double sporeInfDays = 0.00;
            int admissionDuration = 0;
            
            //loop through all days an inpatient - inPtDays[patient][ward] = {times...}
            for(int ward=0; ward < nWards; ward++) {
                if (inPtDays[patient][ward].size()>0) { // if inpatient on that ward
                    for (int t: inPtDays[patient][ward]) {
                        if(t < infTimes[patient]) { //if not yet infected
                            //use wardI to look this up - wardI[t][ward] = nInf
                            wardInfDays += wardI[t][ward];
                            admissionDuration ++;
                            
                            //get infectious numbers from other wards
                            for(int nonWard=0; nonWard < nWards; nonWard++) {
                                if(nonWard!=ward) {
                                    hospInfDays += wardI[t][nonWard];
                                }
                            }
                            
                            //get number of spore day equivalents - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
                            //spores decay according to geometric distribution
                            sporeInfDays += sporeForceSummary[t][ward];
                            
                        }
                    }
                }
            }
            
            int nonAdmissionDuration = infTimes[patient] - 1 - admissionDuration + 1;
            
            double logProbAvoidInf = -((parm.betaBgroundHosp * admissionDuration) //hospital background
                                         + (parm.betaWard * wardInfDays) //ward infection force
                                         + (parm.betaHosp * hospInfDays) //across hospital infection force
                                         + (parm.betaComm * nonAdmissionDuration) //community background
                                         + (parm.betaWard * sporeInfDays ) //ward-based spore force
                                         );
            
            //Pr(infecetd at t)
            int t = infTimes[patient]; //infection time
            
            //if in community
            double betaI;
            if(ptLocation[patient][t]==-1) {
                betaI = parm.betaComm;
            } else {
            //if in hospital
                int ward = ptLocation[patient][t]; //location on when infected
                int wardIt = wardI[t][ward]; //number infectious at this location
                
                int hospIt = 0; //number infectious on other wards
                for(int nonWard=0; nonWard < nWards; nonWard++) {
                    if(nonWard!=ward) {
                        hospIt += wardI[t][nonWard];
                    }
                }
                
                //get number of spore day equivalents - sporeI[t][ward][pt] = spore_duation (numbering from 1...)
                //spores decay according to geometric distribution
                double sporeIt = sporeForceSummary[t][ward];
                
                betaI = parm.betaBgroundHosp + (parm.betaWard * wardIt) + (parm.betaHosp * hospIt) + (parm.betaWard * sporeIt);
            }
            
            double probInf = 1 - exp(-betaI);
            //if(patient==338) printf("%0.4f\t%0.4f\n", probInf, betaI);
            
            //Pr(infected by specific source)
            double probInfSource;
            if (infSourceType[patient] == SrcType::BGROUND_COMM) {
                //background - community
                probInfSource = parm.betaComm/betaI * probInf;
            } else if (infSourceType[patient] == SrcType::BGROUND_HOSP) {
                //background - hospital
                probInfSource = parm.betaBgroundHosp/betaI * probInf;
            } else if (infSourceType[patient] == SrcType::WARD) {
                //ward  - all sources equally likelly so number of infectious inidivudals cancels top and bottom
                probInfSource = parm.betaWard/betaI * probInf;
            } else if (infSourceType[patient] == SrcType::HOSP) {
                //hospital - as above
                probInfSource = parm.betaHosp/betaI * probInf;
            } else if (infSourceType[patient] == SrcType::SPORE) {
                //spore
                int ward = ptLocation[patient][t];
                double specificSporeDuration = sporeI[t][ward][infSources[patient]];
                double specificSporeLevel = pow((1-getSporeP(parm)), specificSporeDuration);
                // prob(infected) * prob(infected by spore) * prob(infected by specific spore), total spore force cancels top and bottom
                probInfSource = probInf * (parm.betaWard  / betaI ) * (specificSporeLevel) ;
                //if(patient==338) printf("%0.4f\t%0.4f\t%0.4f\n", probInf, betaI, specificSporeLevel);
                
                if(specificSporeLevel==0) {
                    printf("patient: %d\n", patient);
                    printf("source %d on ward %d at t=%d (spore duration check %d)\n", infSources[patient], ward, t, sporeI[t][ward][infSources[patient]]);
                    printf("source diagnosis: %d\n source ward stays on ward %d: ", infTimes[infSources[patient]], ward);
                    for(int tt=0; tt<=maxTime; tt++) {
                        if(ptLocation[infSources[patient]][tt]==ward) {
                            printf("%d, ", tt);
                        }
                    }
                    printf("\n");
                }
                
            }
            
            ll += logProbAvoidInf + log(probInfSource) + log(probNotStartInf);

        }
    }
    
    return ll;
}



//log likelihood contribution from sampling times - p(S | I, parm)
double llSample(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, Parm &parm) {
   
    // neg. binomial distributed time between infection and sampling
    double sampleSize = parm.sampleSize;
    double sampleMu = parm.sampleMu;
    double sampleProb = sampleSize / (sampleSize + sampleMu);
    
    double ll = 0;
    //#pragma omp parallel for reduction(+:ll) num_threads(4)
    for(int pt : infectedPatients) {
        //for all infected patients
        int dd = sampleTimes[pt] - infTimes[pt];
        ll += dnbinom(dd, sampleSize, sampleProb, 1);
    }
    
    return ll;
}



//log likelihood contribution for recovery, p(R | I, parm)
double llRecover(vector<int> &infectedPatients, vector<int> &sampleTimes, vector<int> &recTimes, Parm &parm) {
    //neg binomial distributed time between sampling and recovery
    //"recSize", "recMu" = parm.recSize and parm.recMu
    double recSize = parm.recSize;
    double recMu = parm.recMu;
    double recProb = recSize / (recSize + recMu);

    double ll = 0;
    //#pragma omp parallel for reduction(+:ll) num_threads(4)
    for(int pt : infectedPatients) {
        //for those patients with infections
        int dd = recTimes[pt] - sampleTimes[pt];
        ll += dnbinom(dd, recSize, recProb, 1);
        }
    
    return ll;
}



//genetic log likelihood for single patient pair
double llGeneticSingle(vector<int> &infectedPatients, vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, Parm &parm) {
    double ll;
    
    if (transmissionSource==-1) {
        //likelihood for genetic distance based on import:
        
        //for a pair of samples taken from the overall population at the same time -
        //snp ~ Pois(2 * mu * v) and v ~Exp(1/Ne)
        //hence snp ~ Geom(p), where p = p = 1/(1+(1/Ne))
        double p = 1/(1+(1/(2*parm.mu*parm.introNe)));
        
        //have one member of pair of samples, multiple possible ways to choose other member of pair
            // 1. pick one at random - estimate of llGeneticSingle then too unstable
            // 2. determine the mean SNP distance and use this to calculate ll - CURRENTLY implemented
            // 3. determine the mean probability and use the log of this as as the ll - worse performance
            // 4. determine the mean(log(prob)), numerically equivalent to 2.
            // 5. as for 4, and multiply by sum(1/(1:(n-1))) to reflect information available - considerably worse performance

        
        //find all possible comparison samples, the 'nearest neighbours'
            //nearest neighbours defined for now as the first case of every cluster
        vector<int> potentialNN;
        for (int nn : infectedPatients) {
            if (nn!=patient & (infSourceType[nn]== SrcType::BGROUND_HOSP
                               | infSourceType[nn] == SrcType::BGROUND_COMM
                               | infSourceType[nn] == SrcType::START_POS)) {
                potentialNN.push_back(nn);
            }
        }
        int nPotentialNN = (int)potentialNN.size();
        int ptIndex = geneticMap.at(patient);
        
        //get average SNP difference to all potential NN
        double snpSum = 0;
        for(int nn: potentialNN) {
            int srcIndex = geneticMap.at(nn);
            double snp = geneticDist[ptIndex][srcIndex];
            snpSum += snp;
        }
        //use mean SNPs to calculate prob of NN
        double meanSnp = snpSum / nPotentialNN;
        double logProbNN = log(1-p) + meanSnp*log(p);
        ll = logProbNN;
    }
    else {
        //likelihood for direct transmission
        int ptIndex = geneticMap.at(patient);
        int srcIndex = geneticMap.at(transmissionSource);
        double snp = geneticDist[ptIndex][srcIndex];
        double time = fabs(sampleTimes[patient] - sampleTimes[transmissionSource]);
        double Ne = parm.directNe;
        double mu = parm.mu;
        ll = time/2/Ne + logIgamma(1+snp,mu*time+time/2/Ne) + snp*log(2*mu) - logFactorial(snp) - log(Ne) - (snp+1)*log(2*mu+1/Ne);
    }
    return ll;
}






//log likelihood contribution from the genetic distance matrix - p(G | I, S, parm)
double llGenetic(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, Parm &parm) {
    double ll = 0.00;
    //#pragma omp parallel for reduction(+:ll) num_threads(4) schedule(static)
    for (int patient : infectedPatients) {
        //log likelihood for patient not infected = 0
        //therefore, determine log likelihood for infected patients only
        
        if (infSources[patient]==-1) { //infected by hospital or community background or at start
            //likelihood for genetic distance based on import, i.e. source type =0
            double Ne = parm.introNe;
            double mu = parm.mu;
            //snp ~ Pois(2 * mu * v) and v ~Exp(1/Ne)
            //hence snp ~ Geom(p), where p = p = 1/(1+(1/Ne)), and PDF (1-p)p^snp
            double p = 1/(1+(1/(2*mu*Ne)));
            //double llNN = 0;
            vector<int> potentialNN;
            for (int nn : infectedPatients) {
                if (infSources[nn]==-1 & nn!=patient) {
                    //nearest neighbours could be defined as the first case of every cluster - problem is this changes numbers
                    // as move source, therefore link to all except self
                    potentialNN.push_back(nn);
                }
            }
            int nPotentialNN = (int)potentialNN.size();
            int ptIndex = geneticMap.at(patient);

            //get average SNP difference to all potential NN
            double snpSum = 0.0;
            for(int nn: potentialNN) {
                int srcIndex = geneticMap.at(nn);
                double snp = geneticDist[ptIndex][srcIndex];
                snpSum += snp;
            }
            
            //use mean SNPs to calculate prob of NN
            double meanSnp = snpSum / nPotentialNN;
            double logProbNN = log(1-p) + meanSnp*log(p);
            ll += logProbNN;
        }
        
        else {
            //likelihood for direct transmission - ward, hospital or spore
            
            //tight bootleneck forces early coalescence with a host and therefore forces Ne lower
            // - i.e. bottleneck puts upper limit on Ne (most marked effect when generation time small)
            
            //sampling of source long after transmission, provides diversity within 2 hosts and tends to inflate Ne
            
            int sourcePatient = infSources[patient];
            int ptIndex = geneticMap.at(patient);
            int srcIndex = geneticMap.at(sourcePatient);
            double snp = geneticDist[ptIndex][srcIndex];
            double time = fabs(sampleTimes[patient] - sampleTimes[sourcePatient]);
            double Ne = parm.directNe;
            double mu = parm.mu;
            ll += time/2/Ne + logIgamma(1+snp,mu*time+time/2/Ne) + snp*log(2*mu) - logFactorial(snp) - log(Ne) - (snp+1)*log(2*mu+1/Ne);
            
        }
    }
    
    return ll;
}





//prior
double getPrior(Parm &parm) {
    double priorBeta0 = dgamma(parm.betaBgroundHosp, 2, 0.002, 1); //dexp(parm.betaBgroundHosp, 100, 1);
    double priorBeta1 = dgamma(parm.betaWard, 2, 0.002, 1); //dexp(parm.betaWard, 100, 1);
    double priorBeta2 = dgamma(parm.betaHosp, 2, 0.002, 1); //dexp(parm.betaHosp, 100, 1);
    double priorSampleSize = dgamma(parm.sampleSize, 3, 1/0.5, 1); //favours lower values of size parameter
    double priorSampleMu = dgamma(parm.sampleMu, 3, 1/0.1, 1); //up to 100, but most mass around 10-20
    double priorDirectNe = dgamma(parm.directNe, 2, 2, 1); //dexp(parm.directNe, 1, 1);
    double priorIntroNe = dgamma(parm.introNe, 2, 10000, 1); //dexp(parm.introNe, 100, 1);
    double priorMu = dnorm(parm.mu, 2/365.25, 0.05/365.25, 1); //relatively tight prior around 2 SNPs per year
    double priorStartInfLogit = dnorm(parm.probStartInfLogit, 0, 1, 1); //relatively uniform over 0 to 1
    double priorBetaComm = dexp(parm.betaComm, 1, 1);
    double priorSporeProbLogit = dnorm(parm.sporeProbLogit, 0, 2, 1); //relatively uniform over 0 to 1
    double priorRecSize = dnorm(parm.recSize, 3, 0.5, 1); //relatively tight prior around 3, i.e. likely between 2 and 4
    double priorrecMu = dnorm(parm.recMu, 30, 3, 1);   //dgamma(parm.recMu, 5, 1/0.15, 1); //in the 5 to 60 range
    double prior = priorBeta0 + priorBeta1 + priorBeta2 + priorSampleSize + priorSampleMu +
    priorDirectNe + priorIntroNe + priorMu + priorStartInfLogit + priorBetaComm + priorSporeProbLogit + priorRecSize + priorrecMu;
    return(prior);
}

//target distribution, i.e. non-normalised posterior
double targetDist (vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes,
                   vector<int> &infSources, vector<int> &infSourceType,
                   vector<vector<vector<int>>> &sporeI, vector<vector<double>> &sporeForceSummary,
                   vector<vector<vector<int>>> &wardLog,
                   vector<vector<vector<int>>> &inPtDays,
                   vector<vector<int>> &ptLocation,
                   vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime, vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap, Parm &parm) {

    
    double td = llTrans(infectedPatients, uninfectedPatients, infTimes, infSourceType, infSources, sporeI, sporeForceSummary, wardLog, inPtDays, ptLocation, wardI, nPatients, nWards, maxTime, parm) +
                    llSample(infectedPatients, infTimes, sampleTimes, parm) +
                    llRecover(infectedPatients, sampleTimes, recoverTimes, parm) +
                    llGenetic(infectedPatients, infTimes, sampleTimes, infSources, infSourceType, geneticDist, geneticMap, nPatients, parm) +
                    getPrior(parm);
    
    return td;
}
