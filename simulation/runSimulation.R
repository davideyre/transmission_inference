rm(list = ls())
Sys.setenv(PATH="/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")
source('simulation.R')

#dependencies
## R - ape, phangorn, RColorBrewer
## Python - version 2.7, biopython
## seq-gen from https://github.com/rambaut/Seq-Gen, needs to be in path
## may need to comment out over-writing of $PATH above

##output
outDirRoot = "/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/"
randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))


### PARAMETERS ###
#transmission parameters
bgroundBeta = 0.002
wardBeta = 0.005
hospBeta = 0.000
commBeta = 0.001

pStartPos= 0.0

#sampling delay - negative bionmial
sample.size = 5; sample.mu = 10 

#recovery - neg binom distribution
rec.size = 3; rec.mu = 90

#spore decay - geometric distn
spore.p = 0.2
spore.multiplier = 0.5 #the relatively infectiousness of spores

#genetic parameters
directNe = 22.5 #within host population size
introNe = 20000
bottleneck = 1000 #set to 1000 for high bottleneck

#ward set up
nWardsPerHospital = 4; nHospitals = 1

#population set up
nPopulation = 2000
pAdmit = 0.005 #chance of being admitted at any time-step
losMean = 5 #length of stay - poisson distributed, assuming that diseases doesn't impact length of stay

#max time to run for
maxTime = 300

#seq-gen parameters
alnLength = 1000 #(arbitary); 
mutationRate = 1 #mutations per genome per year

simSettings = list(bgroundBeta=bgroundBeta, wardBeta=wardBeta, hospBeta=hospBeta, commBeta=commBeta,
                   pStartPos=pStartPos, sample.size=sample.size, sample.mu=sample.mu, 
                   rec.size=rec.size, rec.mu=rec.mu, 
                   spore.p=spore.p, spore.multiplier=spore.multiplier, 
                   directNe=directNe, introNe=introNe, bottleneck=bottleneck,
                   nWardsPerHospital=nWardsPerHospital, nHospitals=nHospitals, nPopulation=nPopulation, 
                   pAdmit=pAdmit, losMean=losMean,
                   maxTime=maxTime, randomSeed=randomSeed, outDirRoot=outDirRoot,
                   alnLength=alnLength, mutationRate=mutationRate)

runSim(simSettings)