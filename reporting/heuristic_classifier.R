rm(list = ls())

sims = read.csv("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/simulation_list.csv")
simLog = matrix(NA, nrow=nrow(sims), ncol=6)
sim.index=1
max.time = 300
for (simulation in sims[,1]) {
    print(simulation)
    load(paste("~/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/", simulation, "/input/simulation.RData", sep=""))
    
    srcLog = matrix(NA, nrow=length(infectedPatients), ncol=3)
    colnames(srcLog) = c("Patient", "Source", "Route")
    
    for (pt.index in 1:length(infectedPatients)) {
      pt = infectedPatients[pt.index]
      
      #get all prior cases within 2 SNPs
      within.2.snps = gsub("patient_", "", names(which(geneticDist[paste("patient_", pt, sep=""),]<=2)))
      within.2.snps = within.2.snps[which(within.2.snps != pt)]
      
      #restrict to only prior cases
      possible.sources = as.numeric(within.2.snps[which(patientLog[as.numeric(within.2.snps),5]<patientLog[pt,5])])
      sample.date.diff = patientLog[pt,5] - patientLog[possible.sources,5]
      
      
      
      #classify based on those with ward contact, then spore contact, then hospital at same time
      if(length(possible.sources)==0) {
        srcLog[pt.index,] = c(pt, "Background", -1)
      } else {
        src.list = matrix(NA, nrow=length(possible.sources), ncol=5)
        colnames(src.list) = c("src", "sample_dd", "ward_match", "hospita_match", "spore_match")
        src.index = 1
        #look for ward overlap
        for (src in possible.sources) {
          print(src)
          overlap.min = max(1, patientLog[src,5]-14) #sampling date of source
          overlap.max = min(max.time, patientLog[pt,5]) #sampling date of victim
          possible.days = overlap.min:overlap.max #possible days of overlap
          
          match.days = which(ptLocation[possible.days,src] != 0 & 
                               ptLocation[possible.days,src]==ptLocation[possible.days,pt]) #find days in same place, when same place not community (0)
          ward.overlaps = possible.days[match.days]
          print(ward.overlaps)
          #look for hospital overlap
          match.days = which(ptLocation[possible.days,src] != 0 & ptLocation[possible.days,pt] !=0 
                             & ptLocation[possible.days,src] != ptLocation[possible.days,pt])
          hospital.overlaps = possible.days[match.days]
          print(hospital.overlaps)
          #look for spore overlap within 28 days
          #first find any discharges of the src during the window of interest
          discharges = matrix(NA, nrow=1, ncol=2)
          previous.location = ptLocation[overlap.min,src]
          i=1
          for (location in ptLocation[possible.days,src]) {
            #iterate through current locations for the source
            if(location!=previous.location & previous.location!=0) {
              # a discharge has occured
              discharges = rbind(discharges, c(i, previous.location))
            }
            i=i+1
            previous.location = location
          }
          spore.overlaps = vector()
          if(nrow(discharges)>1) {
            discharges[,1] = possible.days[discharges[,1]]
            colnames(discharges) = c("discharge_date", "ward")
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
              discharge = discharges[2,]
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
          srcLog[pt.index,] = c(pt, "Background", -1)
        }
      } #end of check for sources or no sources
      
    }
    true.src.type = as.character(factor(patientLog[infectedPatients,4], levels=c(0,1,2,3,4,5), labels=c("Background", "Ward", "Hospital", "Community", "Start Positive", "Spore")))
    
    srcLog[,1] = as.numeric(srcLog[,1])-1
    srcLog[,3] = ifelse(srcLog[,3]=="-1", -1, as.numeric(srcLog[,3])-1)
    srcLog = as.data.frame(cbind(srcLog, patientLog[infectedPatients,3], ifelse(srcLog[,3]==patientLog[infectedPatients,3],1,0),
                           true.src.type, ifelse(as.character(true.src.type)==as.character(srcLog[,2]),1,0)))
    colnames(srcLog) = c(colnames(srcLog)[1:3], "True Source", "Match Source", "True Source Type", "Match Source Type")
    
    srcLog
    simLog[sim.index,] = c(simulation, length(infectedPatients), sum(as.numeric(as.character(srcLog[,5]))), sum(as.numeric(as.character(srcLog[,5])))/length(infectedPatients),
      sum(as.numeric(as.character(srcLog[,7]))), sum(as.numeric(as.character(srcLog[,7])))/length(infectedPatients))
    sim.index = sim.index +1
}

colnames(simLog) = c("simulation", "n", "src_match", "prop_src_match", "src_type_match", "prop_src_type_match")
write.csv(simLog, "/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/heurstic_trans_stats.csv", row.names = F)
