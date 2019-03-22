library(ape)
library(phangorn)
library(RColorBrewer)

# SI simulation with ward based transmission and background
# Also simulates genetic distances
# for now assume infection times known

### ASSUMPTIONS 
#  discharge times not linked to diagnosis

# MORE KEY ASSUMPTIONS
# Infectious between t_inf +1 to t_rec -1
# Can only shed spores while infectious
# Hence if infected on day of discharge, no spore shed

# Can only infect while infectious

# Can infect on day of discharge


# Can infect on day of discharge, but not after
# 
# 
# Spores set from t_rec onwards
# No onward transmission from spore on day of removal
# 
# 
# 
# IF discharged spore set from day after discharge
# 
# 
# Must be sampled while infectious, or on day of recovery




### NOTES
# see ml_simplified_attemp_fix.R in scatchpad for limitations of discrete time for densely infected outbreaks
# need to ensure relatively low proportion of patients infected
# i.e. that cases relatively sparse compared to admissions


runSim = function(simSettings) {
  
  #set random seed based on settings and import simSettings into namespace
  attach(simSettings)
  set.seed(randomSeed) 
  
  #set up output folders
  simulation = paste("simulation_", randomSeed, sep="")
  outDirBase = paste(outDirRoot, simulation, sep="")
  if(!dir.exists(outDirBase)) {
    dir.create(outDirBase)
    dir.create(paste(outDirBase, "/inference", sep=""))
    outDir = paste(outDirBase, "/input", sep="")
    dir.create(outDir)
  } else {
    print("folder already exists, please try again\n")
    stop()
  }
  
  #save settings
  capture.output(simSettings, file = paste(outDir, "/simulation_settings.txt", sep=""))
  save(simSettings, file = paste(outDir, "/simulation_settings.Rdata", sep=""))
  
  #set up hospital wards
  nWards = nWardsPerHospital * nHospitals
  hospitalWards = matrix(NA, nrow=nHospitals, ncol=nWardsPerHospital)
  for (h in 1:nHospitals) {
    for (w in 1:nWardsPerHospital) {
      wardIndex = (h-1)*nWardsPerHospital + w
      hospitalWards[h,w] = wardIndex
    }
  }
  
  ### RUN SIMULATION ###
  
  #simulate ward admissions, assuming that diseases doesn't impact length of stay
  
  #start with population all in community, rows=time, columns=patients (values = 0 for community, and then 1..nWards for wards)
  ptLocation = matrix(data=0, nrow=maxTime, ncol=nPopulation)
  
  #create level of spore sporeLevel[t][pt][ward]
  sporeLevel = array(0, c(maxTime, nPopulation, nWards))
  
  #at each time step iterate over all patients in community and choose who to admit
  for (t in 2:maxTime) { 
    for (pt in which(ptLocation[t-1,]==0))
      if(runif(1)<pAdmit) {
        #admit patient
        ward = sample(1:nWards,1) #choose a ward
        los = rpois(1, losMean) #determine los
        los_limit = min(maxTime, t+los) #find los
        ptLocation[t:los_limit,pt] = ward
      }
  }
  
  #set up matrix to log infection events
  infections = matrix(NA, nrow=nPopulation, ncol=5) 
  colnames(infections) = c("t_inf", "source", "type", "t_sample", "t_rec")
  
  #assign initial infections
  for(pt in 1:nPopulation) {
    if(runif(1)<pStartPos) {
      #patient starts infected
      print(paste("Patient", pt, "is alredy infected at t=1"))
      #simulate recovery time, accepting for now is approximation (as is sample time for these cases below)
      t.sample = 1 + rnbinom(1, sample.size, mu=sample.mu)
      t.rec = t.sample + rnbinom(1, rec.size, mu=rec.mu)
      infections[pt,] = c(1, -1, 4, t.sample, t.rec)
      
      ## set spores
      #determine if discharged during infectious period
      #get location when infected
      currentLocation = ptLocation[1,pt]
      for(tt in 2:min(t.rec,maxTime)) {
        if(ptLocation[tt,pt]!=currentLocation & currentLocation!=0) {
          #discharge event has occured from a ward, currentLocation, i.e. location at last iteration
          spore.start = tt
          spore.end = maxTime 
          #check if re-admitted to the same ward, if so end the spore at this point if still infectious
          sporeWard = currentLocation
          for(t_chk in spore.start:min(maxTime, t.rec-1)) {
            if(ptLocation[t_chk,pt]==sporeWard) {
              spore.end = t_chk - 1
              break
            }
          }
          if(spore.end>spore.start) {
            spore.ages = 1:(maxTime-tt+1)
            sporeLevel[spore.start:spore.end, pt, currentLocation] = sporeLevel[spore.start:spore.end, pt, currentLocation] + (1-spore.p)^spore.ages
            print(paste("Patient", pt, "has contaminated ward", currentLocation, 
                        "after discharge with spores from", spore.start, "until", spore.end))
          }
          
        }
        currentLocation = ptLocation[tt,pt]
      }
      #determine if recovers while still inpatient
      if(t.rec<=maxTime) {
        recoveryLocation = ptLocation[t.rec, pt]
        if(recoveryLocation>0) {
          spore.start = t.rec
          spore.end = maxTime 
          spore.ages = 1:(spore.end-t.rec+1)
          sporeLevel[spore.start:spore.end, pt, recoveryLocation] = sporeLevel[spore.start:spore.end, pt, recoveryLocation] + (1-spore.p)^spore.ages
          print(paste("Patient", pt, "has contaminated ward", recoveryLocation, 
                      "after recovery with spores from", spore.start, "until", spore.end))
        } 
      }
      
      #end of spore section
    }
  }
  
  #assign other infections
  for(t in 2:maxTime) {
    vectorS = which(is.na(infections[,1]))
    for (pt in vectorS) {
      #patient is still susceptible
      if(ptLocation[t,pt]==0) {
        #patient in community
        if(runif(1)<(1-exp(-commBeta))) {
          #background infection has occured
          srcType = 3 #community background
          src = -1
          t.sample = t + rnbinom(1, sample.size, mu=sample.mu)
          t.rec = t.sample + rnbinom(1, rec.size, mu=rec.mu)
          infections[pt,] = c(t, src, srcType, t.sample, t.rec)
          print(paste("Patient", pt, "is infected in the community at time", t))
          
          #note as set up no spores are shed by patients testing positive in the community
          #can review this
          
        }
      } else {
        #patient in hospital
        
        #get vector of infectious patients on the same ward
        ptWard = ptLocation[t,pt]
        vectorIWard = which(ptLocation[t,]==ptWard & infections[,1]<t & infections[,5]>t) #patients on the same ward infected before current time step
        nIWard = length(vectorIWard)
        
        #get vector of infectious patients on other wards
        #get vector of other wards in same hospital
        hosp = which(hospitalWards==ptWard, arr.ind=T)[1]
        hospWardsPt = hospitalWards[hosp,]
        vectorIHosp = which(ptLocation[t,] %in% hospWardsPt & ptLocation[t,]!=ptWard & ptLocation[t,]!=0 & infections[,1]<t & infections[,5]>t) 
        #patients on a different ward infected before current time step
        nIHosp = length(vectorIHosp)
        
        #get vector of spores that are infectious on the same ward
        vectorISpore = which(sporeLevel[t,,ptWard]>0)
        vectorSporeLevel = sporeLevel[t,vectorISpore,ptWard]
        sumSporeLevel = sum(vectorSporeLevel)
        
        betaI = c(bgroundBeta, (wardBeta*nIWard), (hospBeta*nIHosp), (wardBeta*spore.multiplier*sumSporeLevel)) #overall rate of infection
        totalBetaI = sum(betaI) 
        probSource = betaI/totalBetaI #vector of probabilities for each source type
        
        if(runif(1)<(1-exp(-totalBetaI))) {
          #an infection has occured
          srcType = sample(x=c(0,1,2,5), size=1, prob=probSource) #sample type as ward Bground, ward, hospital, spore = 0,1,2,5
          t.sample = t + rnbinom(1, sample.size, mu=sample.mu)
          #get recovery time
          t.rec = t.sample + rnbinom(1, rec.size, mu=rec.mu)
          
          if(srcType==0) {
            #background infection has occured
            src = -1
            srcType = 0
            print(paste("Patient", pt, "is infected via background in hospital at time", t))
          } else if (srcType==1) {
            #acquisition on the ward
            if(length(vectorIWard)>1) {
              src = sample(vectorIWard, 1) 
            } else {
              src = vectorIWard
            }
            srcType = 1
            print(paste("Patient", pt, "is infected on the same ward,",ptWard,", at time", t, "by", src))
          } else if (srcType==2) {
            #acqusition via hopsital
            if(length(vectorIHosp)>1) {
              src = sample(vectorIHosp, 1) 
            } else {
              src = vectorIHosp
            }
            srcType = 2
            print(paste("Patient", pt, "is infected via the hospital while on,",ptWard,", at time", t, "by", src))
          } else {
            #acquisition from spore
            if(length(vectorISpore)>1) {
              src = sample(vectorISpore, 1, prob=vectorSporeLevel) 
            } else {
              src = vectorISpore
            }
            srcType = 5
            print(paste("Patient", pt, "is infected via spore while on,",ptWard,", at time", t, "by", src))
          }
          infections[pt,] = c(t, src, srcType, t.sample, t.rec)
          
          ## set spores
          if(t<maxTime) {
            currentLocation = ptLocation[t+1,pt] #location on day first become infectious
            
            #determine if discharged during infectious period
            #get location when infected
            for(tt in (t+1):min(t.rec, maxTime)) { #can only drop spores until day before recovery, but go to day of recovery to find discharge in prior time interval
              if(ptLocation[tt,pt]!=currentLocation & currentLocation!=0) {
                #discharge event has occured from a ward, currentLocation, i.e. location at last iteration
                spore.start = tt
                spore.end = maxTime
                #check if re-admitted to the same ward, if so end the spore at this point if still infectious
                sporeWard = currentLocation
                for(t_chk in spore.start:min(maxTime, t.rec-1)) {
                  if(ptLocation[t_chk,pt]==sporeWard) {
                    spore.end = t_chk - 1
                    break
                  }
                }
                if(spore.end>spore.start) {
                  spore.ages = 1:(spore.end-tt+1)
                  sporeLevel[spore.start:spore.end, pt, currentLocation] = sporeLevel[spore.start:spore.end, pt, currentLocation] + (1-spore.p)^spore.ages
                  print(paste("Patient", pt, "has contaminated ward", currentLocation, 
                              "after discharge with spores from", spore.start, "until", spore.end))
                }
              }
              currentLocation = ptLocation[tt,pt]
            }
            if(t.rec<=maxTime) {
              #determine if recovers while still inpatient
              recoveryLocation = ptLocation[t.rec, pt]
              if(recoveryLocation>0) {
                spore.start = t.rec
                spore.end = maxTime 
                spore.ages = 1:(maxTime-t.rec+1)
                sporeLevel[spore.start:spore.end, pt, recoveryLocation] = sporeLevel[spore.start:spore.end, pt, recoveryLocation] + (1-spore.p)^spore.ages
                print(paste("Patient", pt, "has contaminated ward", recoveryLocation, 
                            "after recovery with spores from", spore.start, "until", spore.end))
              }
            }
            #end of spore section
          }
        } #end of inpt infection
      } #end of inpt loop
    } #susceptible pt end
  } #time loop end
  
  
  #infected patients
  infectedPatients = which(!is.na(infections[,1]))
  
  #summarise
  print(paste("Total infected: ", length(infectedPatients), "/", nPopulation, 
              " (", round(100*length(infectedPatients)/nPopulation,1), "%)", sep=""))
  print("Sources:")
  
  #plot summary of simulation
  #ward occupancy
  wardCount = matrix(NA, nrow=maxTime, ncol=(nWards+1))
  colnames(wardCount) = 0:nWards
  for (t in 1:maxTime) {
    for (ward in 0:nWards) {
      wardCount[t,(ward+1)] = sum(ptLocation[t,]==ward)    
    }
  }
  
  #get infectious count
  infCount = matrix(NA, nrow=maxTime, ncol=(nWards+1))
  colnames(wardCount) = 0:nWards
  for (t in 1:maxTime) {
    for (ward in 0:nWards) {
      infCount[t,(ward+1)] = sum(ptLocation[t,]==ward & infections[,1]<t & infections[,5]>t, na.rm=T)    
    }
  }
  
  #percentage infected
  pctCount = round(infCount/wardCount,2)
  pctCount[which(is.na(pctCount))] = 0
  
  # #plot the outbreak - can make more beautiful, but fine for now
  # colNumber = (nWards+1)
  # if (colNumber<3) {colNumber=3}
  # cols = brewer.pal(colNumber,"Set1")
  # 
  # for (wd in 0:nWards) {
  #   if (wd==0) {
  #     plot(pctCount[,1], type="l", col=cols[(wd+1)], ylim=c(0, max(pctCount[2:maxTime,])))
  #   } else {
  #     lines(pctCount[,(wd+1)], col=cols[(wd+1)])
  #   }
  # }
  
  #write log of sources for each infection
  siLog = cbind(infectedPatients, infections[infectedPatients,1:4])
  colnames(siLog) = c("patientid", "t_inf", "source", "source_type", "t_sample")
  
  #write table with log of infected patients, for coalescent simulation
  write.table(siLog, file = paste(outDir, "/siLog.txt", sep=""), 
              row.names = F, sep="\t", quote=F)
  print("siLog file written")
  
  
  ## run coalescent simulationinfectedPatients
  cmd = paste("./simCoalescent.py", 
              " -p ", directNe,
              " -c ", introNe,
              " -b ", bottleneck,
              " -s ", randomSeed,
              " -a ", alnLength,
              " -m ", mutationRate, 
              " -f ", outDir, sep="")
  #system(cmd, ignore.stdout = T, ignore.stderr = T)
  print(cmd)
  system(cmd)
  print("coalescent simulation done")
  
  #read newick file
  treeFile = paste(outDir, "/sim_newick.tree", sep="")
  tree = read.tree(treeFile)
  #convert tip labels back to patient ids
  tree$tip.label = paste("patient_", gsub("_1", "", gsub("Samp_", "", tree$tip.label)), sep="")
  # get distance matrix direct from coalescent tree - improvement would be to estimate tree by ML 
  # tree distances are simply in coalescent time, which has been provided in days
  geneticDist = cophenetic(tree)
  #test
  #geneticDist["269","161"]
  #plot(tree)
  
  write.csv(geneticDist, file=paste(outDir, "/simDistances.csv", sep=""), 
            quote=F)
  write.table(geneticDist, file=paste(outDir, "/simDistances.txt", sep=""),
              quote=F)
  print("genetic distances written")
  
  #write log of all patients, infected and non-infected
  infectionsStr = infections
  infectionsStr[,2] = ifelse(infections[,2]==-1, -1, paste("patient_", infections[,2], sep=""))
  patientLog = cbind(paste("patient_", 1:nPopulation, sep=""), infectionsStr)
  colnames(patientLog) = c("patient_id", "t_inf", "source", "source_type", "t_sample", "t_recover")
  write.csv(patientLog, file=paste(outDir, "/patientLog.csv", sep=""), row.names=FALSE, quote=FALSE)
  
  #admission log - convert back into patient, location, admit_date, discharge_date
  wardLog = c(NA, NA, NA, NA, NA)
  for(pt in 1:nPopulation) {
    prevLocation = 0
    location = NA
    t_admit = NA
    t_discharge = NA
    for (t in 1:maxTime) {
      if(ptLocation[t,pt]!=prevLocation) {
        if(prevLocation==0) {
          #an admission
          t_admit = t
          location = ptLocation[t,pt]
          hosp = which(hospitalWards==location, arr.ind=T)[1]
          prevLocation = location
        } 
        if(ptLocation[t,pt]==0) {
          #a discharge
          t_discharge = t-1
          wardLog = rbind(wardLog, c(paste("patient_", pt, sep=""), location, hosp, t_admit, t_discharge))
          prevLocation = 0
        }
      }
      if(t==maxTime & ptLocation[t,pt]!=0) {
        #record as discharge
        t_discharge = t
        wardLog = rbind(wardLog, c(paste("patient_", pt, sep=""), location, hosp, t_admit, t_discharge))
      }
    }
  }
  #remove top line
  wardLog = wardLog[2:nrow(wardLog),]
  colnames(wardLog) = c("patient_id", "ward", "hospital", "t_admit", "t_discharge")
  write.csv(wardLog, file=paste(outDir, "/wardLog.csv", sep=""), row.names=FALSE, quote=FALSE)
  
  print("simulation log written")
  print("done")
  
  print("infection sources")
  print(table(infections[,3]))
  
  sum(ptLocation!=0) #total inpt daysp
  
  dna.aln = read.dna(paste(outDir, "/sim_alignment.fa", sep=""), format="fasta")
  dna.aln.phyDat = phyDat(dna.aln, type = "DNA", levels = NULL)
  dist.aln = dist.ml(dna.aln.phyDat, model="JC69")
  tree.nj = NJ(dist.aln)
  tree.nj$tip.label = paste("patient_", gsub("_1", "", gsub("Samp_", "", tree.nj$tip.label)), sep="")
  geneticDist = round(cophenetic(tree.nj)*alnLength,0) #option integer distances based on alignment length of 10000
  
  write.csv(geneticDist, file=paste(outDir, "/geneticDistances_snps.csv", sep=""), 
            quote=F)
  write.table(geneticDist, file=paste(outDir, "/geneticDistances_snps.txt", sep=""),
              quote=F)
  tree.nj$edge.length = abs(round(tree.nj$edge.length*1000,1))
  write.tree(tree.nj, file=paste(outDir, "/sim_newick_snps.tree", sep=""))
  print("genetic distances from snps written")
  
  
  ##plot the true transmission tree
  cmd = paste('python plot_transmission_trees.py -f ', 
              outDirBase)
  print(cmd)
  system(cmd)
  print("true transmission tree written")
  
  print("all done")
  
  return(table(infections[,3]))
}













