//
//  struct.cpp
//  inference_test
//
//  Created by David Eyre on 27/01/2017.
//  Copyright © 2017 David Eyre. All rights reserved.
//

#include "struct.hpp"


//operator to allow comparison of transmission events
bool operator==(const trans& lhs, const trans& rhs)
{
    return lhs.t == rhs.t && lhs.victim == rhs.victim;
}



//convert negative/positive number to 0 to 1 scale
double logistic(double x) {
    return 1/(1+exp(-x));
}

//convert 0 to 1 to any real number
double logit(double x) {
    return log(x/(1-x));
}


//function to return spore prob
double getSporeP(Parm &parm) {
    double pSpore = logistic(parm.sporeProbLogit);
    return(pSpore);
}

//function to return start infected prob
double getStartInfP(Parm &parm) {
    double pStartInf = logistic(parm.probStartInfLogit);
    return(pStartInf);
}
