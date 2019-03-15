rm(list = ls())
Sys.setenv(PATH="/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin")
source('simulation.R')

#dependencies
## R - ape, phangorn, RColorBrewer
## Python - version 2.7, biopython
## seq-gen from https://github.com/rambaut/Seq-Gen, needs to be in path
## may need to comment out over-writing of $PATH above


##output
outDirRoot = "/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/5_scenarios/"
nRep = 10


### PARAMETERS ###
#fixed parameters
pStartPos= 0.0
sample.size = 5; sample.mu = 10 #sampling delay - negative bionmial
rec.size = 3; rec.mu = 90 #recovery - neg binom distribution
directNe = 22.5; introNe = 2000; bottleneck = 1000 #set to 1000 for high bottleneck and to match typical within host diversity
nWardsPerHospital = 4; nHospitals = 1 #ward set up
nPopulation = 2000; pAdmit = 0.005; losMean = 5#population set up
maxTime = 300 #max time to run for
alnLength = 1000; mutationRate = 1 #mutations per genome per year, aln-length - arbitary


#scenario 1 - ward transmission only
bgroundBeta = 0.002; wardBeta = 0.005; hospBeta = 0.000; commBeta = 0.000
spore.p = 1; spore.multiplier = 0 #spore decay - geometric distn; the relatively infectiousness of spores

for (i in 1:nRep) {
  randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))
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
}

#scenario 2 - ward and hospital wide transmission
bgroundBeta = 0.002; wardBeta = 0.005; hospBeta = 0.001; commBeta = 0.000
spore.p = 1; spore.multiplier = 0 #spore decay - geometric distn; the relatively infectiousness of spores

for (i in 1:nRep) {
  randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))
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
}

#scenario 3 - ward and hospital wide transmission + short spores
bgroundBeta = 0.002; wardBeta = 0.005; hospBeta = 0.001; commBeta = 0.000
spore.p = 0.3; spore.multiplier = 0.5 #spore decay - geometric distn; the relatively infectiousness of spores

for (i in 1:nRep) {
  randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))
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
}

#scenario 4 - ward and hospital wide transmission + long spores
bgroundBeta = 0.002; wardBeta = 0.004; hospBeta = 0.001; commBeta = 0.000
spore.p = 0.1; spore.multiplier = 0.5 #spore decay - geometric distn; the relatively infectiousness of spores

for (i in 1:nRep) {
  randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))
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
}

#scenario 5 - ward and hospital wide transmission + short spores + community
bgroundBeta = 0.001; wardBeta = 0.003; hospBeta = 0.001; commBeta = 0.0001
spore.p = 0.3; spore.multiplier = 0.5 #spore decay - geometric distn; the relatively infectiousness of spores

for (i in 1:nRep) {
  randomSeed = as.numeric(gsub("0.", "", as.character(round(runif(1),8))))
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
}

