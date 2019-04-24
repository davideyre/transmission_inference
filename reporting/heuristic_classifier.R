#!/usr/bin/env Rscript

rm(list = ls())
library(optparse)

#allow path to be hard coded, but also allow this to be changed at run time
path = "/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenarios/simulation_1723"
path = "/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/5_scenarios/simulation_911145"

#parse command line options
option_list = list(
  make_option( c("-d", "--directory"), type="character", default=path, 
               metavar="character") )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
path = paste(opt$directory, "/input", sep="")

setwd(path)

#read in simulated data
patientLog = read.csv("patientLog.csv", stringsAsFactors=F)
infectedPatients = which(!is.na(patientLog$t_inf))
geneticDist = read.table("geneticDistances_snps.txt")
wardLog = read.csv("wardLog.csv", stringsAsFactors=F)
nPopulation = nrow(patientLog)

#create ptLocation
max.time = max(wardLog[,5])
ptLocation = matrix(data=0, nrow=max.time, ncol=nPopulation)
for (i in 1:nrow(wardLog)) {
  row = wardLog[i,]
  pt = as.numeric(gsub("patient_", "", row[1]))
  loc = paste(row[3], ": ", row[2], sep="")
  adm = as.numeric(row[4])
  dis = as.numeric(row[5])
  ptLocation[adm:dis, pt] = loc
}

#create log of inferred sources
srcLog = matrix(NA, nrow=length(infectedPatients), ncol=3)
colnames(srcLog) = c("Patient", "Source Type", "Source")

# iterate over all patients finding the rule-based source    
for (pt.index in 1:length(infectedPatients)) {
    pt = infectedPatients[pt.index]
    
    #get all prior cases within 2 SNPs
    within.2.snps = gsub("patient_", "", 
                         colnames(geneticDist)[which(geneticDist[paste("patient_", pt, sep=""),]<=2)])
    within.2.snps = within.2.snps[which(within.2.snps != pt)]
      
    #restrict to only prior cases
    possible.sources = within.2.snps[which(patientLog[within.2.snps,"t_inf"] < patientLog[pt,"t_inf"])]
    sample.date.diff = patientLog[pt,"t_inf"] - patientLog[possible.sources,"t_inf"]
      
    #classify based on those with ward contact, then spore contact, then hospital at same time
    if(length(possible.sources)==0) {
        #no genetic match
        #check if admission in 84 days (12 weeks) prior to testing
        window.start = max(1, patientLog[pt, "t_sample"]-84)
        window.end = min(patientLog[pt, "t_sample"], max.time)
        inpt.days = sum(ptLocation[window.start:window.end, pt]!="0")
        if(inpt.days==0) {
          # no inpatient stay for 12 weeks - assign to community
          srcLog[pt.index,] = c(pt, "Community", -1)
        } else {
          # inpatient stay within <=12 weeks - assign to hospital
          srcLog[pt.index,] = c(pt, "Background", -1)
        }
        
    } else {
        src.list = matrix(NA, nrow=length(possible.sources), ncol=5)
        colnames(src.list) = c("src", "sample_dd", "ward_match", "hospita_match", "spore_match")
        src.index = 1
        #look for ward overlap
        for (src in possible.sources) {
          print(src)
          src=as.numeric(src)
          overlap.min = max(1, patientLog[src,"t_inf"]-14) #sampling date of source
          overlap.max = min(max.time, patientLog[pt,"t_inf"]) #sampling date of victim
          possible.days = overlap.min:overlap.max #possible days of overlap
          
          match.days = which(ptLocation[possible.days,src] != 0 & 
                               ptLocation[possible.days,src]==ptLocation[possible.days,pt]) #find days in same place, when same place not community (0)
          ward.overlaps = possible.days[match.days]
          print(ward.overlaps)
          
          #look for hospital overlap - for now assume that only one hospital simulated
          match.days = which(ptLocation[possible.days,src] != 0 & ptLocation[possible.days,pt] !=0 
                             & ptLocation[possible.days,src] != ptLocation[possible.days,pt])
          hospital.overlaps = possible.days[match.days]
          print(hospital.overlaps)
          
          #look for spore overlap within 28 days
          #first find any discharges of the src during the window of interest
          discharges = matrix(NA, nrow=1, ncol=2)
          colnames(discharges) = c("discharge_date", "ward")
          previous.location = ptLocation[overlap.min,src]
          i=1
          for (location in ptLocation[possible.days,src]) {
            #iterate through current locations for the source
            if(location!=previous.location & previous.location!=0) {
              # a discharge has occured
              discharges = rbind(discharges, c(possible.days[i], previous.location))
            }
            i=i+1
            previous.location = location
          }
          
          
          spore.overlaps = vector()
          if(nrow(discharges)>1) {
            if(nrow(discharges)>2) {
              discharges = discharges[2:nrow(discharges),] #discharges - remove blank top line
              for (i in 1:nrow(discharges)) {
                discharge = discharges[i,]
                spore.start = which(possible.days==discharge[1])+1
                spore.end = min(length(possible.days), spore.start+28)
                spores.present = ifelse(1:length(possible.days)>=spore.start & 1:length(possible.days)<=spore.end, discharge[2], 0)
                spores.compare = cbind(possible.days, ptLocation[possible.days,pt], spores.present)
                colnames(spores.compare) = c("possible.days", "victim_location", "spore_location")
                match.days = which(spores.compare[,2]!= 0 & spores.compare[,2] == spores.compare[,3])
                spore.overlaps = c(spore.overlaps, possible.days[match.days])
                print(spore.overlaps)
              }
            }
            else {
              discharge = discharges[2,] #remove blank top line
              spore.start = which(possible.days==discharge[1])+1
              spore.end = min(length(possible.days), spore.start+28)
              spores.present = ifelse(1:length(possible.days)>=spore.start & 1:length(possible.days)<=spore.end, discharge[2], 0)
              spores.compare = cbind(possible.days, ptLocation[possible.days,pt], spores.present)
              colnames(spores.compare) = c("possible.days", "victim_location", "spore_location")
              match.days = which(spores.compare[,2]!= 0 & spores.compare[,2] == spores.compare[,3])
              spore.overlaps =  possible.days[match.days]
              print(spore.overlaps)
            }
          } #end of discharges list
          print(spore.overlaps)
          
          src.list[src.index,] = c(src, sample.date.diff[src.index], ifelse(length(ward.overlaps)>0,1,0),
                                   ifelse(length(hospital.overlaps)>0,1,0),
                                   ifelse(length(spore.overlaps)>0,1,0))
          src.index = src.index +1
        } #end of src list
        
        
        ward.sources = which(src.list[,3]==1)
        hospital.sources = which(src.list[,4]==1)
        spore.sources = which(src.list[,5]==1)
        
        if(length(ward.sources)>=1) {
          if(length(ward.sources)==1) {
            srcLog[pt.index,] = c(pt, "Ward", src.list[ward.sources,1])
          } else {
            min.date = min(src.list[ward.sources,2])
            srcPt = which(src.list[ward.sources,2]==min.date)[1]
            srcLog[pt.index,] = c(pt, "Ward", src.list[srcPt,1])
          }
        } else if(length(spore.sources)>=1) {
          if(length(spore.sources)==1) {
            srcLog[pt.index,] = c(pt, "Spore", src.list[spore.sources,1])
          } else {
            min.date = min(src.list[spore.sources,2])
            srcPt = which(src.list[spore.sources,2]==min.date)[1]
            srcLog[pt.index,] = c(pt, "Spore", src.list[srcPt,1])
          }
        } else if(length(hospital.sources)>=1) {
          if(length(hospital.sources)==1) {
            srcLog[pt.index,] = c(pt, "Hospital", src.list[hospital.sources])
          } else {
            min.date = min(src.list[hospital.sources,2])
            srcPt = which(src.list[hospital.sources,2]==min.date)[1]
            srcLog[pt.index,] = c(pt, "Hospital", src.list[srcPt,1])
          }
        } else {
          #check if admission in 84 days (12 weeks) prior to testing
          window.start = max(1, patientLog[pt, "t_sample"]-84)
          window.end = min(patientLog[pt, "t_sample"], max.time)
          inpt.days = sum(ptLocation[window.start:window.end, pt]!="0")
          if(inpt.days==0) {
            # no inpatient stay for 12 weeks - assign to community
            srcLog[pt.index,] = c(pt, "Community", -1)
          } else {
            # inpatient stay within <=12 weeks - assign to hospital
            srcLog[pt.index,] = c(pt, "Background", -1)
          }
        }
      } #end of check for sources or no sources
      
    }
    



true.src.type = as.character(factor(patientLog[infectedPatients,"source_type"], levels=c(0,1,2,3,4,5), labels=c("Background", "Ward", "Hospital", "Community", "Start Positive", "Spore")))

    #srcLog[,1] = as.numeric(srcLog[,1])-1
    #srcLog[,3] = ifelse(srcLog[,3]=="-1", -1, as.numeric(srcLog[,3])-1)
    srcLogMerge = as.data.frame(cbind(srcLog, gsub("patient_", "", patientLog[infectedPatients,"source"]), 
                                      ifelse(srcLog[,3]==gsub("patient_", "", patientLog[infectedPatients,"source"]),1,0),
                           true.src.type, ifelse(as.character(true.src.type)==as.character(srcLog[,2]),1,0)))
    colnames(srcLogMerge) = c(colnames(srcLog)[1:3], "True Source", "Match Source", "True Source Type", "Match Source Type")

    srcLogMerge

write.csv(srcLogMerge, "../inference/heurstic_trans_stats.csv", row.names = F)
