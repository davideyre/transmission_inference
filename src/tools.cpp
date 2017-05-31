
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
        for (int t = infTimes[sporePt]; t<=maxTime; t++) { //can only add spore after infected, hence start from there
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

//function to return log factorial
double logFactorial(int n) {
    if (1==1) { //temp over-ride for non-integer values
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*3.141592653589793238462643383279502884) + 1.0/(12.0*x);
    }
    if (n > 254)
    {
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*3.141592653589793238462643383279502884) + 1.0/(12.0*x);
    }
    else
    {
        vector<double> lf =
        {
            0.000000000000000,
            0.000000000000000,
            0.693147180559945,
            1.791759469228055,
            3.178053830347946,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.604602902745251,
            12.801827480081469,
            15.104412573075516,
            17.502307845873887,
            19.987214495661885,
            22.552163853123421,
            25.191221182738683,
            27.899271383840894,
            30.671860106080675,
            33.505073450136891,
            36.395445208033053,
            39.339884187199495,
            42.335616460753485,
            45.380138898476908,
            48.471181351835227,
            51.606675567764377,
            54.784729398112319,
            58.003605222980518,
            61.261701761002001,
            64.557538627006323,
            67.889743137181526,
            71.257038967168000,
            74.658236348830158,
            78.092223553315307,
            81.557959456115029,
            85.054467017581516,
            88.580827542197682,
            92.136175603687079,
            95.719694542143202,
            99.330612454787428,
            102.968198614513810,
            106.631760260643450,
            110.320639714757390,
            114.034211781461690,
            117.771881399745060,
            121.533081515438640,
            125.317271149356880,
            129.123933639127240,
            132.952575035616290,
            136.802722637326350,
            140.673923648234250,
            144.565743946344900,
            148.477766951773020,
            152.409592584497350,
            156.360836303078800,
            160.331128216630930,
            164.320112263195170,
            168.327445448427650,
            172.352797139162820,
            176.395848406997370,
            180.456291417543780,
            184.533828861449510,
            188.628173423671600,
            192.739047287844900,
            196.866181672889980,
            201.009316399281570,
            205.168199482641200,
            209.342586752536820,
            213.532241494563270,
            217.736934113954250,
            221.956441819130360,
            226.190548323727570,
            230.439043565776930,
            234.701723442818260,
            238.978389561834350,
            243.268849002982730,
            247.572914096186910,
            251.890402209723190,
            256.221135550009480,
            260.564940971863220,
            264.921649798552780,
            269.291097651019810,
            273.673124285693690,
            278.067573440366120,
            282.474292687630400,
            286.893133295426990,
            291.323950094270290,
            295.766601350760600,
            300.220948647014100,
            304.686856765668720,
            309.164193580146900,
            313.652829949878990,
            318.152639620209300,
            322.663499126726210,
            327.185287703775200,
            331.717887196928470,
            336.261181979198450,
            340.815058870798960,
            345.379407062266860,
            349.954118040770250,
            354.539085519440790,
            359.134205369575340,
            363.739375555563470,
            368.354496072404690,
            372.979468885689020,
            377.614197873918670,
            382.258588773060010,
            386.912549123217560,
            391.575988217329610,
            396.248817051791490,
            400.930948278915760,
            405.622296161144900,
            410.322776526937280,
            415.032306728249580,
            419.750805599544780,
            424.478193418257090,
            429.214391866651570,
            433.959323995014870,
            438.712914186121170,
            443.475088120918940,
            448.245772745384610,
            453.024896238496130,
            457.812387981278110,
            462.608178526874890,
            467.412199571608080,
            472.224383926980520,
            477.044665492585580,
            481.872979229887900,
            486.709261136839360,
            491.553448223298010,
            496.405478487217580,
            501.265290891579240,
            506.132825342034830,
            511.008022665236070,
            515.890824587822520,
            520.781173716044240,
            525.679013515995050,
            530.584288294433580,
            535.496943180169520,
            540.416924105997740,
            545.344177791154950,
            550.278651724285620,
            555.220294146894960,
            560.169054037273100,
            565.124881094874350,
            570.087725725134190,
            575.057539024710200,
            580.034272767130800,
            585.017879388839220,
            590.008311975617860,
            595.005524249382010,
            600.009470555327430,
            605.020105849423770,
            610.037385686238740,
            615.061266207084940,
            620.091704128477430,
            625.128656730891070,
            630.172081847810200,
            635.221937855059760,
            640.278183660408100,
            645.340778693435030,
            650.409682895655240,
            655.484856710889060,
            660.566261075873510,
            665.653857411105950,
            670.747607611912710,
            675.847474039736880,
            680.953419513637530,
            686.065407301994010,
            691.183401114410800,
            696.307365093814040,
            701.437263808737160,
            706.573062245787470,
            711.714725802289990,
            716.862220279103440,
            722.015511873601330,
            727.174567172815840,
            732.339353146739310,
            737.509837141777440,
            742.685986874351220,
            747.867770424643370,
            753.055156230484160,
            758.248113081374300,
            763.446610112640200,
            768.650616799717000,
            773.860102952558460,
            779.075038710167410,
            784.295394535245690,
            789.521141208958970,
            794.752249825813460,
            799.988691788643450,
            805.230438803703120,
            810.477462875863580,
            815.729736303910160,
            820.987231675937890,
            826.249921864842800,
            831.517780023906310,
            836.790779582469900,
            842.068894241700490,
            847.352097970438420,
            852.640365001133090,
            857.933669825857460,
            863.231987192405430,
            868.535292100464630,
            873.843559797865740,
            879.156765776907600,
            884.474885770751830,
            889.797895749890240,
            895.125771918679900,
            900.458490711945270,
            905.796028791646340,
            911.138363043611210,
            916.485470574328820,
            921.837328707804890,
            927.193914982476710,
            932.555207148186240,
            937.921183163208070,
            943.291821191335660,
            948.667099599019820,
            954.046996952560450,
            959.431492015349480,
            964.820563745165940,
            970.214191291518320,
            975.612353993036210,
            981.015031374908400,
            986.422203146368590,
            991.833849198223450,
            997.249949600427840,
            1002.670484599700300,
            1008.095434617181700,
            1013.524780246136200,
            1018.958502249690200,
            1024.396581558613400,
            1029.838999269135500,
            1035.285736640801600,
            1040.736775094367400,
            1046.192096209724900,
            1051.651681723869200,
            1057.115513528895000,
            1062.583573670030100,
            1068.055844343701400,
            1073.532307895632800,
            1079.012946818975000,
            1084.497743752465600,
            1089.986681478622400,
            1095.479742921962700,
            1100.976911147256000,
            1106.478169357800900,
            1111.983500893733000,
            1117.492889230361000,
            1123.006317976526100,
            1128.523770872990800,
            1134.045231790853000,
            1139.570684729984800,
            1145.100113817496100,
            1150.633503306223700,
            1156.170837573242400,
        };
        return lf[n];
    }
}




