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
        if (infTimes[patient] == 0) {
            //patients infected at start
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
    // epsilon = parm.sampleEpsilon
    double ll = 0.00;
    //#pragma omp parallel for reduction(+:ll) num_threads(4)
    for(int i : infectedPatients) {
        //for all infected patients
        int dd = sampleTimes[i] - infTimes[i];
        if(infTimes[i]>0) {
            //infected after t=0
            ll += dpois(dd, parm.sampleEpsilon, 1);
        }
        else {
            int ddNotInfected = dd - 1; //were not infected until before this interval, e.g. sample t=2, infected t<=0
            ll += log(1 - ppois(ddNotInfected, parm.sampleEpsilon, 1, 0));
        }
    }
    return ll;
}



//log likelihood contribution for recovery, p(R | I, parm)
double llRecover(vector<int> &infectedPatients, vector<int> &sampleTimes, vector<int> &recTimes, Parm &parm) {
    //"recSize", "recMu" = parm.recSize and parm.recMu
    double recSize = parm.recSize;
    double recMu = parm.recMu;
    double recProb = recSize / (recSize + recMu);

    double ll = 0;
    //#pragma omp parallel for reduction(+:ll) num_threads(4)
    for(int i : infectedPatients) {
        //for those patients with infections
        int dd = recTimes[i] - sampleTimes[i];
        ll += dnbinom(dd, recSize, recProb, 1);
        }
    return ll;
}



//genetic log likelihood for single patient pair
double llGeneticSingle(vector<int> &infectedPatients, vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, Parm &parm) {
    double ll;
    if (transmissionSource==-1) {
        //likelihood for genetic distance based on import, i.e. source type =0
        double Ne = parm.introNe;
        double mu = parm.mu;
        //snp ~ Pois(2 * mu * v) and v ~Exp(1/Ne)
        //hence snp ~ Geom(p), where p = p = 1/(1+(1/Ne))
        double p = 1/(1+(1/(2*mu*Ne)));
        vector<int> potentialNN;
        for (int nn : infectedPatients) {
            if (nn!=patient & (infSourceType[nn]== SrcType::BGROUND_HOSP | infSourceType[nn] == SrcType::BGROUND_COMM)) {
                //nearest neighbours defined for now as any case (may wish to consider only comparing to the first case of every cluster
                potentialNN.push_back(nn);
            }
        }
        
        
        
        int nPotentialNN = (int)potentialNN.size();
        int ptIndex = geneticMap.at(patient);
        /*
        double probNN = 0;
        for(int nn: potentialNN) {
            int srcIndex = geneticMap.at(nn);
            double snp = geneticDist[ptIndex][srcIndex];
            probNN += (1-p) * pow(p,snp);

        }
        ll = log(probNN/nPotentialNN);
         */
        
        
        
        //get average SNP difference to all potential NN
        double snpSum = 0.0;
        for(int nn: potentialNN) {
            int srcIndex = geneticMap.at(nn);
            double snp = geneticDist[ptIndex][srcIndex];
            snpSum += snp;
        }
        
        //use mean SNPs to calculate prob of NN
        double meanSnp = snpSum / nPotentialNN;
        double probNN = (1-p) * pow(p,meanSnp);
        
        
        /*
        
        //pick a NN at random and get snp distance
        int nnRandomIndex = floor(runif(0, nPotentialNN));
        int nnRandom = potentialNN[nnRandomIndex];
        int srcIndex = geneticMap.at(nnRandom);
        double snp = geneticDist[ptIndex][srcIndex];
        double probNN = (1-p) * pow(p,snp);
        */
        
        ll = log(probNN);
        
    }
    else {
        //likelihood for direct transmission
        int ptIndex = geneticMap.at(patient);
        int srcIndex = geneticMap.at(transmissionSource);
        double snp = geneticDist[ptIndex][srcIndex];
        double time = fabs(sampleTimes[patient] - sampleTimes[transmissionSource]);
        double Ne = parm.directNe;
        double mu = parm.mu;
        
        
        if(snp>160) {
            //over-ride to deal with factorial having upper limit, this will be an issue for distances >170, can scale though if needed in inputs
            ll = -numeric_limits<double>::infinity();
        }
        else {
            //ll = log(exp(time/2/Ne+log(igamma(1+snp,mu*time+time/2/Ne)))
            //          * pow(2*mu, snp)
            //          / factorial(snp)/Ne/pow((2*mu+1/Ne),(snp+1)));
            ll = time/2/Ne +
                    log(igamma(1+snp,mu*time+time/2/Ne)) +
                    snp*log(2*mu) -
                    log(factorial(snp)) -
                    log(Ne) -
                    (snp+1)*log(2*mu+1/Ne);

        }
        /*
        ll = time/2/Ne +
        log(igamma(1+snp,mu*time+time/2/Ne)) +
        snp*log(2*mu) -
        logFactorial(snp) -
        log(Ne) -
        (snp+1)*log(2*mu+1/Ne);
        */
        
    }

    
    return ll;
}


//log likelihood over all pairs calling the llGeneticSingle function (slower x4, but used for testing only)
double llGeneticAlt(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, Parm &parm) {
    
    double ll = 0.00;
    for (int patient=0; patient < nPatients; patient++) {
        //log likelihood for patient not infected = 0
        //therefore, determine log likelihood for infected patients
        if(infTimes[patient]>-1) { //infected patients only
            if (infSourceType[patient]== SrcType::BGROUND_HOSP |
                infSourceType[patient]== SrcType::BGROUND_COMM |
                infSourceType[patient]== SrcType::START_POS) { //infected by hospital or community background or at start
                ll += llGeneticSingle(infectedPatients, sampleTimes, patient, -1, infSourceType, geneticDist, geneticMap, nPatients, parm);
            }
            else {
                //likelihood for direct transmission
                int sourcePatient = infSources[patient];
                ll += llGeneticSingle(infectedPatients, sampleTimes, patient, sourcePatient, infSourceType, geneticDist, geneticMap, nPatients, parm);
                
            }
        }
    }
    return ll;
}



//log likelihood contribution from the genetic distance matrix - p(G | I, S, parm)
double llGenetic(vector<int> &infectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType, vector<vector<double>> &geneticDist, unordered_map<int,int> &geneticMap, int nPatients, Parm &parm) {
    double ll = 0.00;
    #pragma omp parallel for reduction(+:ll) num_threads(4) schedule(static)
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
            /*
            double probNN = 0;
            for(int nn: potentialNN) {
                int srcIndex = geneticMap.at(nn);
                double snp = geneticDist[ptIndex][srcIndex];
                probNN += (1-p) * pow(p,snp);
                //llNN += log(1-p) + snp*log(p) - log(nPotentialNN);
                //llNN += log((1-p) * pow(p,snp) * 1/nPotentialNN); //change form to be numerically more robust

            }
            //printf("%0.6f\n", probNN);
             ll += log(probNN/nPotentialNN);
             */
            
            
            //get average SNP difference to all potential NN
            double snpSum = 0.0;
            for(int nn: potentialNN) {
                int srcIndex = geneticMap.at(nn);
                double snp = geneticDist[ptIndex][srcIndex];
                snpSum += snp;
            }
            
            //use mean SNPs to calculate prob of NN
            double meanSnp = snpSum / nPotentialNN;
            double probNN = (1-p) * pow(p,meanSnp);
             
            /*
            
            int nnRandomIndex = floor(runif(0, nPotentialNN));
            int nnRandom = potentialNN[nnRandomIndex];
            int srcIndex = geneticMap.at(nnRandom);
            double snp = geneticDist[ptIndex][srcIndex];
            double probNN = (1-p) * pow(p,snp);
             */
            
            ll += log(probNN);

            //ll += llNN;
        }
        
        else {
            //likelihood for direct transmission - ward, hospital or spore
            
            //tight bootleneck forces early coalescnce and therefore forces Ne lower - i.e. bottleneck puts upper limit on Ne (most marked effect when generation time small)
            
            //sampling of source long after transmission, provides diversity within 2 hosts and tends to inflate Ne
            
            int sourcePatient = infSources[patient];
            int ptIndex = geneticMap.at(patient);
            int srcIndex = geneticMap.at(sourcePatient);
            double snp = geneticDist[ptIndex][srcIndex];
            double time = fabs(sampleTimes[patient] - sampleTimes[sourcePatient]);
            double Ne = parm.directNe;
            double mu = parm.mu;
            
            
            if(snp>170) {
                //over-ride to deal with factorial having upper limit, this will be an issue for distances >170, can scale though if needed in inputs
                ll += -numeric_limits<double>::infinity();
            }
            else {
                //ll += log(exp(time/2/Ne+log(igamma(1+snp,mu*time+time/2/Ne)))
                //          * pow(2*mu, snp)
                //          / factorial(snp)/Ne/pow((2*mu+1/Ne),(snp+1)));
                ll += time/2/Ne +
                        log(igamma(1+snp,mu*time+time/2/Ne)) +
                        snp*log(2*mu) -
                        log(factorial(snp)) -
                        log(Ne) -
                        (snp+1)*log(2*mu+1/Ne);
            }
            
            /*
            ll += time/2/Ne +
            log(igamma(1+snp,mu*time+time/2/Ne)) +
            snp*log(2*mu) -
            logFactorial(snp) -
            log(Ne) -
            (snp+1)*log(2*mu+1/Ne);
            */
        }
    }
    
    return ll;
}





//prior
double getPrior(Parm &parm) {
    double priorBeta0 = dgamma(parm.betaBgroundHosp, 2, 0.002, 1); //dexp(parm.betaBgroundHosp, 100, 1);
    double priorBeta1 = dgamma(parm.betaWard, 2, 0.002, 1); //dexp(parm.betaWard, 100, 1);
    double priorBeta2 = dgamma(parm.betaHosp, 2, 0.002, 1); //dexp(parm.betaHosp, 100, 1);
    double priorEpsilon = dgamma(parm.sampleEpsilon, 3, 1/0.5, 1);; //rgamma(1000,3,0.5)
    double priorDirectNe = dgamma(parm.directNe, 2, 2, 1); //dexp(parm.directNe, 1, 1);
    double priorIntroNe = dgamma(parm.introNe, 2, 10000, 1); //dexp(parm.introNe, 100, 1);
    double priorMu = dnorm(parm.mu, 2/365.25, 0.05/365.25, 1); //relatively tight prior around 2 SNPs per year
    double priorStartInfLogit = dnorm(parm.probStartInfLogit, 0, 1, 1);
    double priorBetaComm = dexp(parm.betaComm, 1, 1);
    double priorSporeProbLogit = dnorm(parm.sporeProbLogit, 0, 2, 1); //dunif(parm.sporeProbLogit, 0, 1, 1); //dbeta(parm.sporeProbLogit, 2.5, 5, 1); ///dgamma(parm.sporeProbLogit, 5, 0.06, 1);
    //dunif(parm.sporeProbLogit, 0, 1, 1); //dgamma(parm.sporeProbLogit, 10, 0.03, 1);//for gamma specify by shape and scale (note default in R is shape and rate), shape*scale=mean
    //dnorm(parm.sporeProbLogit, 0.15, 0.05, 1); //spore prob between 0.05 and 0.3, spore durations tail up to 15-100 days
    //double priorBetaSpore = dgamma(parm.betaSpore, 2, 0.004, 1); //dexp(parm.betaSpore, 100, 1);
    double priorRecSize = dnorm(parm.recSize, 3, 0.5, 1); //relatively tight prior around 3, i.e. likely between 2 and 4
    double priorrecMu = dnorm(parm.recMu, 30, 3, 1);   //dgamma(parm.recMu, 5, 1/0.15, 1); //in the 5 to 60 range
    double prior = priorBeta0 + priorBeta1 + priorBeta2 + priorEpsilon +
    priorDirectNe + priorIntroNe + priorMu + priorStartInfLogit + priorBetaComm + priorSporeProbLogit + /*priorBetaSpore + */ priorRecSize + priorrecMu;
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

