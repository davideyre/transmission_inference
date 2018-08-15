# To do

## priorities
1. Review behaviour testing with two STs - problem is 3 tiers of diversity - whole dataset, within ST, transmission events , e.g. ST42 + ST10 or ST42+ST11

*** Review the need for the starting infected paramters



## Testing
1. Suite of test scenarios to repeat including all possibilities - initially try 1 of each


## Code
1. For now all input times are assumed to be numbered from 1, allow for conversion from date format
2. assume for now genetic distances provided for all samples - not explicitly checked - could add this check


## Improvements
1. Refine genetic likelihood


## Application
1. Apply to WGS1400 dataset


## Possibles
1. Allow for false negative tests, i.e. augmentation of sampled, but negative cases
2. Display ESS with each block of chain outputted to screen
3. The t_samp - t_inf interval follows the same distribtution in and out of hospital, whereas may be more likely to get abx and therefore symptoms, as well as being more likely to be sampled in hospital - this means that for patients diagnosed in the community the time infection occured may be biased to be too near the point of sampling, potentially
