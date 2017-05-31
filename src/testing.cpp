//
//  testing.cpp
//  inference_test
//
//  Created by David Eyre on 06/01/2017.
//  Copyright © 2017 David Eyre. All rights reserved.
//

#include "testing.hpp"

void runTest(vector<int> &infectedPatients, vector<int> &uninfectedPatients, vector<int> &infTimes, vector<int> &sampleTimes, vector<int> &recoverTimes, vector<double> &parm,
             vector<int> &infSources, vector<int> &infSourceType,
             vector<vector<vector<int>>> &sporeI, vector<vector<vector<double>>> &sporeForce, vector<vector<double>> &sporeForceSummary,
             vector<vector<double>> &geneticDist, unordered_map<int,int> geneticMap, vector<vector<vector<int>>> &wardLog,
             vector<vector<vector<int>>> &inPtDays,
             vector<vector<int>> &ptLocation,
             vector<vector<int>> &wardI, int nPatients, int nWards, int maxTime) {
    
    //timer start
    clock_t start,end;
    start = clock();
    
    //sample log likelihood
    double testSampleLL = llSample(infectedPatients, infTimes, sampleTimes, parm);
    printf("TEST: SampleLL = %0.4f (Expected value: -235.9437)\n\n", testSampleLL);
    
    //genetic log likelihood
    double testGeneticLL = llGenetic(infectedPatients, infTimes, sampleTimes, infSources, infSourceType, geneticDist, geneticMap, nPatients, parm);
    printf("TEST: GeneticLL = %0.4f (Expected value: -30344.3)\n\n", testGeneticLL);
    
    //test of single log genetic likelihood over all pairs
    double testGeneticSingleLL = llGeneticAlt(infectedPatients, infTimes, sampleTimes, infSources, infSourceType,
                                              geneticDist, geneticMap, nPatients, parm);
    printf("TEST: GeneticSingleLL = %0.4f (Expected value: -30344.3)\n\n", testGeneticSingleLL);
    
    //test of priors
    double testPrior = getPrior(parm);
    printf("TEST: Prior = %0.6f (Expected value: -9.559376)\n\n", testPrior);
    
    //test of transmission LL
    double testTransLL = llTrans(infectedPatients, uninfectedPatients, infTimes, infSourceType, infSources, sporeI, sporeForceSummary, wardLog,
                                 inPtDays, ptLocation, wardI, nPatients, nWards, maxTime, parm);
    printf("TEST: TransLL = %0.6f (Expected value: -1957.525)\n\n", testTransLL);
    
    //test of target distribution
    double testTarget = targetDist(infectedPatients, uninfectedPatients, infTimes, sampleTimes, recoverTimes, infSources, infSourceType, sporeI, sporeForceSummary,
                                   wardLog, inPtDays, ptLocation, wardI,
                                   nPatients, nWards, maxTime, geneticDist, geneticMap,
                                   parm);
    printf("TEST: TargetDist = %0.6f (Expected value: -32547.32)\n\n", testTarget);
    
    
    //report timing
    end = clock();			//Stop the clock.
    double timer = (double(end)-double(start))/CLOCKS_PER_SEC;
    printf("TEST Processing time: %0.6f seconds\n\n",timer);
    
}
