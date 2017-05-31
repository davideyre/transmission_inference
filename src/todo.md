To do
=====

Code
----
1. Switch parm to be a struct for ease of reading code
2. Switch infSrcType to be enum, and move out of struct.h
3. Convert handling of never infected cases to be ward totals by day


Inference
---------

1. Update duration from infection to sampling to follow Gamma distn rather than Poisson
2. Update handling of directNe and introNe to reflect recent discussion


Testing
-------
1. Testing starting positive simualtions, and related MCMC updates


Application
-----------
1. Apply to WGS1400 dataset


Possibles
---------
1. Allow for false negative tests, i.e. augmentation of sampled, but negative cases
