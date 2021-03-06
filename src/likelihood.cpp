//
//  likelihood.cpp
//  inference_test
//
//  Created by David Eyre on 07/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#include "likelihood.hpp"

//log likelihood contribution from transmission model - p(I | parm)

double llTrans(vector<vector<int>> &wardEver, vector<vector<int>> &hospitalWards, vector<int> &ward2Hospital, vector<vector<int>> &hospitalWardList,
               vector<int> &infTimes, vector<int> &recTimes, vector<int> &infSourceType, vector<int> &infSources,
               vector<vector<vector<SporeEvent>>> &sporePatientI, vector<vector<double>> &sporeForceSummary,
               vector<vector<vector<int>>> &wardLogInf, vector<vector<int>> &wardLogNeverInf,
               vector<vector<vector<int>>> &inPtDays,
               vector<vector<int>> &ptLocation,
               vector<vector<int>> &wardI, int nInfPatients, int nNeverInfPatients, int nWards, int maxTime, Parm &parm) {

    double ll = 0.00;

    int wardInfDays = 0; // patient - days of ward pressure exposed to
    int hospInfDays = 0; // hospital pressure
    double sporeInfDays = 0.00; // spore pressure
    int admissionDuration = 0; // total days of inpatient stay
    
    //get total number of patients infectious in each hospital at each time point - hospitalI[t][hospital]
    vector<vector<int>> hospitalI;
    int nHospitals = (int)hospitalWardList.size();
    hospitalI.resize(maxTime+1);
    for(int t=0; t<=maxTime; t++) {
        hospitalI[t].resize(nHospitals);
        for( int hosp = 0; hosp<nHospitals; hosp++) {
            hospitalI[t][hosp] = 0;
            for (int ward : hospitalWardList[hosp]) {
                hospitalI[t][hosp] += wardI[t][ward];
            }
        }
    }
    
    // ** LL for never infected patients ** //
    
    // Pr(not infected at start) //
    double probNotStartInf = 1 - getStartInfP(parm);
    ll += nNeverInfPatients*log(probNotStartInf);
    
    //loop over each time step - get numbers infectious and at risk who are never infected
    for(int t=0; t<=maxTime; t++) {

      
        //loop over wards
         for(int ward=0; ward < nWards; ward++) {
             
             //log patients admitted that day admission duration for background exposure
             admissionDuration += wardLogNeverInf[t][ward];
             
             //ward infected person-days exposed to - wardI[t][ward] = nInf
             wardInfDays += wardI[t][ward] * wardLogNeverInf[t][ward]; // wardLogNeverInf[time][ward] = int for count of never infected patients
             
             //get hospital for current ward
             int hosp = ward2Hospital[ward];
             
             //get infectious numbers from other wards
             int nonWardInfs = hospitalI[t][hosp] - wardI[t][ward];
             hospInfDays += nonWardInfs * wardLogNeverInf[t][ward];
             
             //get number of spore day equivalents - from sporeForceSummary
             sporeInfDays += sporeForceSummary[t][ward] * wardLogNeverInf[t][ward];
             
         }
        
         
    }

    int nonAdmissionDuration = (nNeverInfPatients*maxTime) - admissionDuration;

    double logProbAvoidInf = -((parm.betaBgroundHosp * admissionDuration) //hospital background
                                         + (parm.betaWard * wardInfDays) //ward infection force
                                         + (parm.betaHosp * hospInfDays) //across hospital infection force
                                         + (parm.betaComm * nonAdmissionDuration) //community background
                                         + (parm.betaWard * getSporeMultiplier(parm) * sporeInfDays ) //ward-based spore force
                                         );
    ll += logProbAvoidInf;

    
    // ** loop over infected patients ** //
    
    for(int patient=0; patient<nInfPatients; patient++) {
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
            
            //loop through all wards have ever been on using wardEver[pt] = {ward1, ward2, ...}
            for (int ward : wardEver[patient]) {
                int hosp = ward2Hospital[ward]; //get hospital ward on
                for (int t: inPtDays[patient][ward]) { //loop through all days an inpatient - inPtDays[patient][ward] = {times...}
                    if(t < infTimes[patient]) { //if not yet infected
                        //use wardI to look this up - wardI[t][ward] = nInf
                        wardInfDays += wardI[t][ward];
                        admissionDuration ++;
                        
                        //get infectious numbers from other wards
                        hospInfDays += hospitalI[t][hosp] - wardI[t][ward];
                        
                        //get number of spore day equivalents - from sporeForceSummary
                        sporeInfDays += sporeForceSummary[t][ward];
                    }
                }
            }
            
            int nonAdmissionDuration = infTimes[patient] - 1 - admissionDuration + 1;
            
            double logProbAvoidInf = -((parm.betaBgroundHosp * admissionDuration) //hospital background
                                         + (parm.betaWard * wardInfDays) //ward infection force
                                         + (parm.betaHosp * hospInfDays) //across hospital infection force
                                         + (parm.betaComm * nonAdmissionDuration) //community background
                                         + (parm.betaWard * getSporeMultiplier(parm) * sporeInfDays ) //ward-based spore force
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
                for(int nonWard : hospitalWards[ward]) {
                    hospIt += wardI[t][nonWard];
                }
                
                //get number of spore day equivalents - from sporeForceSummary
                double sporeIt = sporeForceSummary[t][ward];
                
                betaI = parm.betaBgroundHosp + (parm.betaWard * wardIt) + (parm.betaHosp * hospIt) + (parm.betaWard * getSporeMultiplier(parm) * sporeIt);
            }
            
            double probInf = 1 - exp(-betaI);
            //if(patient==338) printf("%0.4f\t%0.4f\n", probInf, betaI);
            
            //Pr(infected by specific source)
            double probInfSource = 0;
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
                double specificSporeLevel = 0;
                for(SporeEvent sporeEvent: sporePatientI[ward][infSources[patient]]) { //for each time spore left on this ward
                    if( sporeEvent.start <= t & sporeEvent.end >=t) { //find the correct spore event
                        int specificSporeDuration = t - sporeEvent.start + 1; //get number of days since spore set, including day set
                        specificSporeLevel += pow((1-getSporeP(parm)), specificSporeDuration);
                    }
                }
                
                // prob(infected) * prob(infected by spore) * prob(infected by specific spore), total spore force cancels top and bottom
                probInfSource = probInf * ((getSporeMultiplier(parm) * parm.betaWard)  / betaI ) * (specificSporeLevel) ;
                //if(patient==338) printf("%0.4f\t%0.4f\t%0.4f\n", probInf, betaI, specificSporeLevel);
                
                if(specificSporeLevel==0) {
                    printf("--------------------------\n");
                    printf("spore duration error\n\n");
                    printf("victim: %d\n", patient);
                    printf("infected t=%d via spore on ward %d\n", infTimes[patient], ward);
                    printf("victim ward stays on ward %d: ", ward);
                    for(int tt=0; tt<=maxTime; tt++) {
                        if(ptLocation[patient][tt]==ward) {
                            printf("%d, ", tt);
                        }
                    }
                    printf("\n\n");
                    printf("source: %d\n", infSources[patient]);
                    printf("infected t=%d\n", infTimes[infSources[patient]]);
                    printf("recovered t=%d\n", recTimes[infSources[patient]]);
                    printf("source ward stays on ward %d: ", ward);
                    for(int tt=0; tt<=maxTime; tt++) {
                        if(ptLocation[infSources[patient]][tt]==ward) {
                            printf("%d, ", tt);
                        }
                    }
                    
                    printf("\n--------------------------\n\n");
                }
                
            }
            
            ll += logProbAvoidInf + log(probInfSource) + log(probNotStartInf);

        }
    }
    
    return ll;
}



//log likelihood contribution from sampling times - p(S | I, parm)
double llSample(int nInfPatients, vector<int> &infTimes, vector<int> &sampleTimes, Parm &parm) {
   
    // neg. binomial distributed time between infection and sampling
    double sampleSize = getSampleSize(parm);
    double sampleMu = parm.sampleMu;
    double sampleProb = sampleSize / (sampleSize + sampleMu);
    
    double ll = 0;
    for(int pt=0; pt<nInfPatients; pt++) {
        //for all infected patients
        int dd = sampleTimes[pt] - infTimes[pt];
        ll += dnbinom(dd, sampleSize, sampleProb, 1);
    }
    
    return ll;
}



//log likelihood contribution for recovery, p(R | I, parm)
double llRecover(int nInfPatients, vector<int> &sampleTimes, vector<int> &recTimes, Parm &parm) {
    //neg binomial distributed time between sampling and recovery
    //"recSize", "recMu" = parm.recSize and parm.recMu
    double recSize = parm.recSize;
    double recMu = parm.recMu;
    double recProb = recSize / (recSize + recMu);

    double ll = 0;
    for(int pt=0; pt<nInfPatients; pt++) {
        //for those patients with infections
        int dd = recTimes[pt] - sampleTimes[pt];
        ll += dnbinom(dd, recSize, recProb, 1);
        }
    
    return ll;
}

#ifndef no_genetic

//genetic log likelihood for single patient pair
double llGeneticSingle(vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm) {
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
        for (int nn=0; nn<nInfPatients; nn++) {
            if ((nn!=patient) & (infSourceType[nn]== SrcType::BGROUND_HOSP
                               | infSourceType[nn] == SrcType::BGROUND_COMM
                               | infSourceType[nn] == SrcType::START_POS)) {
                potentialNN.push_back(nn);
            }
        }
        
        //get average SNP difference to all potential NN
        double snpSum = 0;
        int sumNN = 0; //keep track of if have any nearest neighbours within 3000 SNPs
        
        for(int nn : potentialNN) {
            double snp = geneticDist[patient][nn];
            if(snp<=3000) { //restrict to only consider nearest neighbours within 3000 snps
                snpSum += snp;
                sumNN ++;
            }
        }
        
        //use mean SNPs to calculate prob of NN
        double meanSnp;
        if (sumNN>0) {
            meanSnp = snpSum / sumNN;
        }
        else {
            meanSnp = 3000; // if no nearest neighbours set the meanSnp to be the upper limit of 3000
        }
        
        double logProbNN = log(1-p) + meanSnp*log(p);
        ll = logProbNN;
    }
    else {
        //likelihood for direct transmission
        double snp = geneticDist[patient][transmissionSource];
        double time = fabs(sampleTimes[patient] - sampleTimes[transmissionSource]);
        double Ne = parm.directNe;
        double mu = parm.mu;
        ll = time/2/Ne + logIgamma(1+snp,mu*time+time/2/Ne) + snp*log(2*mu) - logFactorial(snp) - log(Ne) - (snp+1)*log(2*mu+1/Ne);
    }
    return ll;
}






//log likelihood contribution from the genetic distance matrix - p(G | I, S, parm)
double llGenetic(vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType,
                 vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm) {
    double ll = 0.00;
    
    //determine all possible nearest neighbours for imported cases, i.e. all other imported cases
    vector<int> nnList;
    nnList.reserve(nInfPatients);
    for (int nn =0; nn<nInfPatients; nn++) {
        if (infSources[nn]==-1) {
            nnList.push_back(nn);
        }
    }
    
    
    for (int patient =0; patient<nInfPatients; patient++) {
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
            
            //get average SNP difference to all potential NN
            double snpSum = 0.0;
            int sumNN = 0; //keep track of if have any nearest neighbours within 3000 SNPs
            
            for (int nn : nnList) {
                if (nn!=patient) {
                    double snp = geneticDist[patient][nn];
                    if(snp<=3000) { //restrict to only consider nearest neighbours within 3000 snps
                        snpSum += snp;
                        sumNN ++;
                    }
                }
            }

            //use mean SNPs to calculate prob of NN
            double meanSnp;
            if (sumNN>0) {
                meanSnp = snpSum / sumNN;
            }
            else {
                meanSnp = 3000;
            }
            double logProbNN = log(1-p) + meanSnp*log(p);
            ll += logProbNN;
        }
        
        else {
            //likelihood for direct transmission - ward, hospital or spore
            
            //tight bootleneck forces early coalescence with a host and therefore forces Ne lower
            // - i.e. bottleneck puts upper limit on Ne (most marked effect when generation time small)
            
            //sampling of source long after transmission, provides diversity within 2 hosts and tends to inflate Ne
            
            int sourcePatient = infSources[patient];
            double snp = geneticDist[patient][sourcePatient];
            double time = fabs(sampleTimes[patient] - sampleTimes[sourcePatient]);
            double Ne = parm.directNe;
            double mu = parm.mu;
            ll += time/2/Ne + logIgamma(1+snp,mu*time+time/2/Ne) + snp*log(2*mu) - logFactorial(snp) - log(Ne) - (snp+1)*log(2*mu+1/Ne);
            
        }
    }
    
    return ll;
}

#else
double llGenetic(vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &infSources, vector<int> &infSourceType,
                 vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm) {
    double ll = 0.00;
    return ll;
}

double llGeneticSingle(vector<int> &sampleTimes, int patient, int transmissionSource, vector<int> infSourceType,
                       vector<vector<double>> &geneticDist, int nInfPatients, Parm &parm) {
    double ll = 0.00;
    return ll;
}
    

#endif



//prior
double getPrior(Parm &parm) {
    //set gamma priors on beta values
    double priorBeta0 = dgamma(parm.betaBgroundHosp, 2, 0.002, 1);
    double priorBeta1 = dgamma(parm.betaWard, 2, 0.002, 1);
    double priorBeta2 = dgamma(parm.betaHosp, 2, 0.002, 1);
    double priorBetaComm = dgamma(parm.betaComm, 2, 0.002, 1);
    
    //interval from infection to sampling - set to be relatively imformative
    //double priorSampleSize = dgamma(parm.sampleSize, 3, 1/0.5, 1); //favours lower values of size parameter
    double priorSampleSize = dnorm(parm.sampleSize, 0.7, 0.5, 1); //most weight between 0 and 10
    //double priorSampleSize = dnorm(parm.sampleSize, -1, 0.01, 1); //tight to improve convergence as test
    double priorSampleMu = dgamma(parm.sampleMu, 3, 1/0.1, 1); //up to 100, but most mass around 10-30
    
    //genetic priors - base directNe on within host diversity estimates, and my on previous SNPs per year
    double priorDirectNe = dnorm(parm.directNe, 74.03717, 18.50929, 1); //within host diversity = 2.mu.Ne
    // try to match nejm range = 0.30 (0.13-0.42), actually 0.3 (0.15 - 0.45) below
    // 2.mu.Ne = 0.3, mu = 0.00202601, Ne = 0.3 / 0.004052019 = 74.03717 sd = 0.075/0.004052019 = 18.50929
    
    double priorIntroNe = dgamma(parm.introNe, 2, 10000, 1); //dexp(parm.introNe, 100, 1);
    
    double priorMu = dnorm(parm.mu, 0.8/365.25, 0.295/365.25, 1);
    // match nejm rate 0.74 (0.22 - 1.40) as close as possible mean = 0.8, 95% CI 0.21 - 1.39 per year, and convert to be per day
    
    //set logit transformed probabilties to be relatively uniform over 0 to 1, set priors on transformed scale to avoid Jacobian
    double priorSporeProbLogit = dnorm(parm.sporeProbLogit, 0, 1.7, 1); //relatively uniform over 0 to 1
    double priorSporeMultiplierLogit = dnorm(parm.sporeMultiplier, 0, 1.7, 1); //relatively uniform over 0 to 1

    //probability of starting infected, set to be relatively informative, i.e. that starting infected is rare (with CDI, not just carriage)
    double priorStartInfLogit = dnorm(parm.probStartInfLogit, -8, 1, 1); //mostly <1/2000 with some mass up to around 1/200
    
    //set recovery prior to be moderately informative
    double priorRecSize = dnorm(parm.recSize, 3, 0.5, 1); //relatively tight prior around 3, i.e. likely between 2 and 4
    double priorrecMu = dnorm(parm.recMu, 90, 15, 1); // recovery, mean 90 days, range 60-120 for mean
    
    double prior = priorBeta0 + priorBeta1 + priorBeta2 + priorSampleSize + priorSampleMu +
                    priorDirectNe + priorIntroNe + priorMu + priorStartInfLogit + priorBetaComm +
                    priorSporeProbLogit + priorSporeMultiplierLogit + priorRecSize + priorrecMu;
    return(prior);
}

//target distribution, i.e. non-normalised posterior
double targetDist (vector<vector<int>> &wardEver, vector<vector<int>> &hospitalWards, vector<int> &ward2Hospital, vector<vector<int>> &hospitalWardList,
                   vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes,
                   vector<int> &infSources, vector<int> &infSourceType,
                   vector<vector<vector<SporeEvent>>> &sporePatientI, vector<vector<double>> &sporeForceSummary,
                   vector<vector<vector<int>>> &wardLogInf, vector<vector<int>> &wardLogNeverInf,
                   vector<vector<vector<int>>> &inPtDays,
                   vector<vector<int>> &ptLocation,
                   vector<vector<int>> &wardI, int nInfPatients, int nNeverInfPatients, int nWards, int maxTime, vector<vector<double>> &geneticDist, Parm &parm) {

    
    double td = llTrans(wardEver, hospitalWards, ward2Hospital, hospitalWardList, infTimes, recoverTimes, infSourceType, infSources, sporePatientI, sporeForceSummary, wardLogInf, wardLogNeverInf, inPtDays, ptLocation, wardI, nInfPatients, nNeverInfPatients, nWards, maxTime, parm) +
                    llSample(nInfPatients, infTimes, sampleTimes, parm) +
                    llRecover(nInfPatients, sampleTimes, recoverTimes, parm) +
                    llGenetic(infTimes, sampleTimes, infSources, infSourceType, geneticDist, nInfPatients, parm) +
                    getPrior(parm);
    
    return td;
}
