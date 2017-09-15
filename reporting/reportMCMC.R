#!/usr/bin/env Rscript

rm(list = ls())
library(coda)
library(optparse)
library(ggplot2)
library(gridExtra)
library(reshape)
library(RColorBrewer)


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#allow path to be hard coded, but also allow this to be changed at run time
path = "/Users/davideyre/Dropbox/Transmission_Inference/xcode_project/nejm_data_st42/"

#read in genetic distances
geneticDist = read.table(file=paste(path, "input/geneticDistances_snps.txt", sep=""))

#read in sample times
patientLog = read.csv(file = paste(path, "input/patientLog.csv", sep=""), stringsAsFactors=F)
infectedPatients = which(!is.na(patientLog$t_sample))
sampleTimes = patientLog$t_sample[infectedPatients]

#read in parameters
mcmcLog = paste(path, "inference/chain_parameters.txt", sep="")
chain = read.table(mcmcLog, header=T)

burnIn = nrow(chain) * 0.2

logistic = function(x) {return(1/(1+exp(-x)))}
finalChain = chain[burnIn:nrow(chain),]
finalChain = cbind(finalChain, logistic(finalChain$spore_prob_logit), logistic(finalChain$p_start_inf_logit))
#plot(as.mcmc(finalChain))

parmChainFile = paste(path, "inference/parm_plots.pdf", sep="")
pdf(parmChainFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

plot(as.mcmc(finalChain[,]),ask=FALSE)
dev.off()

#get mean
mean = apply( finalChain , 2 , mean )

#get 95% HPD interval
hpd = HPDinterval(as.mcmc(finalChain), prob=0.95)

#get ESS
ess = effectiveSize(as.mcmc(finalChain))

#print parameter summary and save to file
parmSummary = cbind(mean, hpd, ess)
print(parmSummary)
parmSummary = cbind(rownames(parmSummary), parmSummary)
colnames(parmSummary) = c("parameter", "mean", "lower", "upper", "ess")
parmSummary = as.data.frame(parmSummary, row.names=1:nrow(parmSummary))
parmFile = paste(path, "inference/parm_summary.csv", sep="")
write.csv(parmSummary, parmFile, row.names=F)

#analyse infection times
infTimeLog = paste(path, "inference/chain_inf_times.txt", sep="")
inf.infTimes = read.table(infTimeLog, header=T)
inf.infTimes = inf.infTimes[burnIn:nrow(inf.infTimes),]
inf.meanInfTimes = apply( inf.infTimes , 2 , mean )
inf.meanInfTimes = inf.meanInfTimes +1 #convert back to numbering from 1 rather than zero

#get ESS for infection times
ess.infTimes = effectiveSize(as.mcmc(inf.infTimes))
#get HPD
inf.infTimes.hpd = HPDinterval(as.mcmc(inf.infTimes))  +1 #convert back to numbering from 1 rather than zero
infTimes.summary = cbind(inf.meanInfTimes, inf.infTimes.hpd)

#analyse recovery times
recTimeLog = paste(path, "inference/chain_rec_times.txt", sep="")
rec.recTimes = read.table(recTimeLog, header=T)
rec.recTimes = rec.recTimes[burnIn:nrow(rec.recTimes),]
rec.meanrecTimes = apply( rec.recTimes , 2 , mean )
rec.meanrecTimes = rec.meanrecTimes +1 #convert back to numbering from 1 rather than zero

#get ESS for recection times
ess.recTimes = effectiveSize(as.mcmc(rec.recTimes))
#get HPD
rec.recTimes.hpd = HPDinterval(as.mcmc(rec.recTimes))  +1 #convert back to numbering from 1 rather than zero
recTimes.summary = cbind(rec.meanrecTimes, rec.recTimes.hpd)



##plot ESS and differences for infection times and recovery times
infRecFile = paste(path, "inference/inf_rec_plots.pdf", sep="")
pdf(infRecFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
par(mfrow = c(2,2))
hist(ess.infTimes)
hist(sampleTimes-inf.meanInfTimes)
hist(ess.recTimes)
hist(rec.meanrecTimes-sampleTimes)
dev.off()



#### source type analysis ####
infSourceTypesLog = paste(path, "inference//chain_inf_source_types.txt", sep="")
all.inf.infSourceTypes = read.table(infSourceTypesLog, header=T)
rowN = nrow(all.inf.infSourceTypes)
inf.infSourceTypes = all.inf.infSourceTypes[burnIn:rowN,]


#a - row per patient, column per transmission type
a = matrix(0, nrow=ncol(inf.infSourceTypes), ncol=6)
for (i in 1:ncol(inf.infSourceTypes)) {
  a[i,] = cbind(length(which(inf.infSourceTypes[,i]==0)), 
                length(which(inf.infSourceTypes[,i]==1)), 
                length(which(inf.infSourceTypes[,i]==2)),
                length(which(inf.infSourceTypes[,i]==3)),
                length(which(inf.infSourceTypes[,i]==4)),
                length(which(inf.infSourceTypes[,i]==5))
  )
}

#MAP esitmate
infSourceType.map = max.col(a)-1
infSourceType.summary = cbind(infSourceType.map, a)

#HPD estimate - function not used for now
rowCount = nrow(inf.infSourceTypes)
typeHPD = function(l, rowCount) {
  #get 95% cutoff
  cutoff = rowCount*0.95
  #get number of columsn to include
  colKeep = which(cumsum(sort(l, decreasing=T))<cutoff) + 1
  if (length(colKeep) == 0) {
    colKeep = 1
  }
  colKeep.max = max(colKeep)
  types = which(l %in% sort(l, decreasing=T)[1:colKeep.max])
  return (types-1)
}


#generate summary of source types
factoriseSrc = function(srcList) {
   o = factor(srcList, levels=c(0, 1,5,2,3,4), 
              labels=c("Hospital background", "Ward", "Spore", "Hospital-wide", "Community background", "Start positive"))
   return(o)
}

src.compare.map = factoriseSrc(infSourceType.map)

d = matrix(NA, ncol=6, nrow=nrow(inf.infSourceTypes))
for (i in 1:nrow(inf.infSourceTypes)) {
   d[i,] = table(factoriseSrc(inf.infSourceTypes[i,]))
}

summary.srctype = cbind(table(src.compare.map), apply( d , 2 , mean ), HPDinterval(as.mcmc(d), prob=0.95))
summary.srctype = cbind(rownames(summary.srctype), summary.srctype)
colnames(summary.srctype) = c("source_type", "map", "mean", "lower", "upper")
srcTypeFile = paste(path, "inference/source_type_summary.csv", sep="")
write.csv(summary.srctype, srcTypeFile, row.names=F)



### EXACT SOURCE MATCHES
infSourceLog = paste(path, "inference/chain_inf_sources.txt", sep="")
all.inf.infSources = read.table(infSourceLog, header=T)
rowN = nrow(all.inf.infSources)
inf.infSources = all.inf.infSources[burnIn:rowN,]
inf.infSources.mode = apply(inf.infSources, 2, Mode)

rowCount = nrow(inf.infSources)
sourceHPD = function(x, rowCount) {
  #get 95% cutoff
  cutoff = rowCount*0.95

  ux = unique(x)
  freq = tabulate(match(x, ux))
  colKeep = which(cumsum(sort(freq, decreasing=T))<cutoff) + 1
  if (length(colKeep) == 0) {
    colKeep = 1
  }
  colKeep.max = max(colKeep)
  sources = ux[match(sort(freq, decreasing=T)[1:colKeep.max], freq)]
  return (sources)
}

infSrc.summary = cbind(names(inf.infSources.mode), inf.infSources.mode)
colnames(infSrc.summary) = c("recipient", "mode_source")
rownames(infSrc.summary) = NULL

nCases = length(inf.infSources.mode)

snpList = rep(NA, nCases)

for (i in 1:nCases) {
   if(inf.infSources.mode[i]=="-1") snpList[i] = NA
   else snpList[i] = geneticDist[infSrc.summary[i,1], infSrc.summary[i,2]]
}


#get the minimum number of SNPs to each case
minSNPs = rep(NA, nCases)
for (i in 1:nCases) {
   minSNPs[i] = min(geneticDist[as.character(infSrc.summary[i,1]), which(colnames(geneticDist)!=as.character(infSrc.summary[i,1]))])
}

infSrc.summary = as.data.frame(cbind(infSrc.summary, infSourceType.map, snpList, minSNPs, ess.infTimes))
colnames(infSrc.summary) = c("recipient", "mode_source", "mode_source_type","SNPs_to_mode", "minSNPs", "ESS_inf_times")

srcFile = paste(path, "inference/source_summary.csv", sep="")
write.csv(infSrc.summary, srcFile, row.names=F)








## create plots of infection sources and times


getRoute = function(srcType) {
  o = factor(srcType, levels=c(0, 1,2,3,4,5),
             labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore"))
  return(as.character(o))
}

colors = brewer.pal(6,"Spectral")

#save plots as list
plot=list()
i=1
for (pt in infectedPatients) {
  d = as.data.frame(table(inf.infSources[,pt], inf.infSourceTypes[,pt]))
  colnames(d)[1:2] = c("Source", "SourceType")
  d = melt(d, id=c("Source", "SourceType"))
  d = cbind(d[,c(1,2,4)], d[,4]/sum(d[,4]))
  colnames(d) = c("Source", "SourceType", "Frequency", "Probability")
  d = d[which(d$Frequency>0),]
  plot[[i]] = ggplot(d, aes(x=factor(Source), y=Probability,
                            fill=factor(SourceType, levels=c(0, 1,2,3,4,5),
                                        labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore")))) +
    geom_bar(stat="identity") +
    scale_y_continuous(limits=c(0,1)) +
    labs(y="Posterior probability", x="Source",
         title=paste("Infection sources for patient ", patientLog$patient_id[pt], sep="")) +
    scale_fill_manual(values=colors, drop=F, name="Source Type") +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=10),
          axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10))
  i=i+1

  #plot inferred infection times
  inf.times = inf.infTimes[,pt] +1 #convert back to number from 1
  d = table(inf.times, inf.infSourceTypes[,pt])
  dMin = min(inf.times)
  dMax = max(inf.times)
  d = as.data.frame(cbind(as.numeric(rownames(d)), d))
  colnames(d)[1] = "inf_times"
  rownames(d) = 1:nrow(d)
  d = melt(d, id=c("inf_times"))
  d = cbind(d, d[,3]/sum(d[,3]))
  colnames(d) = c("InfectionTimes", "SourceType", "Frequency", "Probability")
  d = d[which(d$Frequency>0),]
  plot[[i]] = ggplot(d, aes(x=InfectionTimes, y=Probability,
                            fill=factor(SourceType, levels=c(0, 1,2,3,4,5),
                                        labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore")))) +
    geom_bar(stat="identity") +
    scale_y_continuous(limits=c(0,1)) +
    labs(y="Posterior probability", x="Estimated infection time",
         title=paste("Infection times for patient ", patientLog$patient_id[pt], sep="")) +
    scale_fill_manual(values=colors, drop=F, name="Source Type") +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=10),
          axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10)) +
    scale_x_continuous(breaks=dMin:dMax)
  i = i+1
}

#save plots to pdf - 4 per page
srcPlotFile = paste(path, "inference/source_plots.pdf", sep="")
glist = lapply(plot, ggplotGrob)
ggsave(srcPlotFile, marrangeGrob(glist, nrow = 2, ncol = 2, top=""), width=29.7/2.54, height=21/2.54, useDingbats=FALSE)



#read in wardLog
wardLog = read.csv(file = paste(path, "input/wardLog.csv", sep=""), stringsAsFactors=F)

#read in SNPs
snpList

#get wardLog for list of patients
plotWard = function(searchString, wardLog, wardMinMax, expand, patientLog) {
   wardLogSelectPatient = which(wardLog$patient_id %in% searchString)
   wardLogSubset = wardLog[wardLogSelectPatient,]
   
   patientLogSelectPatient = which(patientLog$patient_id %in% searchString)
   patientLogSubset = patientLog[patientLogSelectPatient,]
   
   wardMinMax = c(wardMinMax[1]-expand,wardMinMax[2]+expand)
   #restrict to min max
   df.ward = wardLogSubset[which(wardLogSubset$t_discharge>=wardMinMax[1] & wardLogSubset$t_admit<=wardMinMax[2]),]
   
   df.ward.sample = cbind(df.ward, patientLog$t_sample[match(df.ward$patient_id, patientLog$patient_id)])
   colnames(df.ward.sample) = c(colnames(df.ward), "t_sample")
   
   
   
   #df.ward.sample$t_sample[which(df.ward.sample$t_sample>df.ward.sample$t_discharge | df.ward.sample$t_sample<df.ward.sample$t_admit)] = NA
   
   print(df.ward.sample)
   print(patientLogSubset)
   
   p = ggplot(data=df.ward.sample) +
      geom_errorbarh(mapping=aes(y=factor(ward), x=t_admit, xmin=t_admit, xmax=t_discharge, 
                                               color=factor(patient_id)), height=0.4, size=1) +
      geom_point(mapping=aes(y=factor(ward), x=t_sample, color=factor(patient_id)), size=3, shape=4, stroke=1.5)
   print(p)
}

plotWardWrapper = function(ptString, wardLog, patientLog, inf.infSources, inf.infTimes) {
   expand = 50
   pt = which(patientLog$patient_id==ptString)
   sources = table(inf.infSources[,pt])
   ptList = c(ptString, names(sources))
   inf.times = inf.infTimes[,pt] +1
   wardMinMax = c(min(inf.times), max(inf.times))
   plotWard(ptList, wardLog, wardMinMax, expand, patientLog)
}


plotWardWrapper(ptString = "C00006232", wardLog, patientLog, inf.infSources, inf.infTimes)
