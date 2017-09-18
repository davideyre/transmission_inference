To do
=====

Update spore handling
 - when readmitted reset spore, i.e. no spore while admitted and only one spore running at any one time

THIS IS PARTIALLY DONE - NEED TO FINISH UPDATING THE PART THAT SETS SPORES

Testing
-------

1. Suite of test scenarios to repeat including all possibilities - initially try 1 of each


Code
----
1. For now all input times are assumed to be numbered from 1, allow for conversion from date format
2. assume for now genetic distances provided for all samples - not explicitly checked - could add this check


Improvements
------------
1. Refine genetic likelihood


Application
-----------
1. Apply to WGS1400 dataset


Possibles
---------
1. Allow for false negative tests, i.e. augmentation of sampled, but negative cases
2. Display ESS with each block of chain outputted to screen
