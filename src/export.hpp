//
//  export.hpp
//  inference_test
//
//  Created by David Eyre on 08/09/2016.
//  Copyright © 2016 David Eyre. All rights reserved.
//

#ifndef export_hpp
#define export_hpp

#include <fstream>
#include <vector>

using namespace std;

void exportChain(vector<vector<double>> &chain, vector<vector<int>> &chainInfTimes, vector<vector<int>> &chainInfSources, vector<vector<int>> &chainInfSourceTypes,
                 vector<vector<int>> &chainRecTimes, int steps, string filePath);

#endif /* export_hpp */
