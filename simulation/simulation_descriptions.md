Simulation scenarios
====================

Scenario types include:

* Background and Ward
* Background, Ward and Spore (short spore)
* Background, Ward and Spore (long spore)
* Background, Ward and Hospital
* Background, Ward, Hospital and Spore (short spore)


Common parameters
-----------------

Common parameters:

* Population=2000
* Duration=300
* losMean=5
* nWards=4
* probAdmission=0.005

These allow for approx 20-30 patients on a ward at once.

For now fix genetic parameters to produce something fairly similar to C diff, ~2 SNPs / yr, directNe=1, and introNe = 10000, with tight bottleneck (=1000).

Recovery fixed to follow negative binomial distribution (size=3, mean=30).


Scenario parameters
-------------------

Beta parameters for each scenario:

* Background and Ward
 - bground = 0.002
 - ward = 0.01
 - hosp = 0
 - p.spore = 1
* Background, Ward and Spore (short spore)
 - bground = 0.002
 - ward = 0.01
 - hosp = 0
 - p.spore = 0.4
* Background, Ward and Spore (long spore)
 - bground = 0.002
 - ward = 0.005
 - hosp = 0
 - p.spore = 0.1
* Background, Ward and Hospital
 - bground = 0.002
 - ward = 0.01
 - hosp = 0.001
 - p.spore = 1
* Background, Ward, Hospital and Spore (short spore)
 - bground = 0.002
 - ward = 0.007
 - hosp = 0.001
 - p.spore = 0.4
 
Test strategy
-------------
Generate 5 simulations for each scenario. (For now generate only one MCMC chain for each simulation, but this can be returned to, once can combine chains.)

Run each MCMC chain for 20,000 iterations initially.

Compare simulated and estimated (point estimated and 95% HPD):

* values of beta coefficients and p.spore
* values of epsilon (ignoring other estimates for now)
* transmission inference - infection time in 95% HPD, infection source type - point estimate and 95% HPD, infection source - point estimate and 95% HPD

May be helpful to generate detailed outputs from one simulation, e.g. histogram of most likely sources and source types (for colour) for a given case.

Also consider comparing performance to SNP based cut-off heuristic for most likely source:

* all previously sampled cases within 2 SNPs and 90 days considered as possible donors
* favour ward < spore < hospital
* when multiple possibilities pick the source with the smallest sampling date difference, and at random for ties 

