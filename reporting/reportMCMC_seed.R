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

logistic = function(x) {return(1/(1+exp(-x)))}


runReport = function(pathRoot, seed, burnIn.factor, thin.factor) {
  
  #read in genetic distances
  geneticDist = read.table(file=paste(pathRoot, "input/geneticDistances_snps.txt", sep=""))
  
  #read in sample times
  patientLog = read.csv(file = paste(pathRoot, "input/patientLog.csv", sep=""), stringsAsFactors=F)
  infectedPatients = which(!is.na(patientLog$t_sample))
  sampleTimes = patientLog$t_sample[infectedPatients]
  
  #read in wardLog
  wardLog = read.csv(file = paste(pathRoot, "input/wardLog.csv", sep=""), stringsAsFactors=F)
  
  #read in parameters
  mcmcLog = paste(pathRoot, "inference/", seed, "/chain_parameters.txt", sep="")
  chain = read.table(mcmcLog, header=T)
  
  burnIn = (nrow(chain) * burnIn.factor) +1
  
  finalChain = chain[burnIn:nrow(chain),]
  finalChain = cbind(finalChain, exp(finalChain$sampleSize), logistic(finalChain$spore_prob_logit), logistic(finalChain$p_start_inf_logit))
  finalChain = mcmc(finalChain[seq(1, nrow(finalChain), by = thin.factor),], start=burnIn, end=nrow(chain), thin=thin.factor)
  
  parmChainFile = paste(pathRoot, "inference/", seed, "/parm_plots.pdf", sep="")
  pdf(parmChainFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
  
  plot(finalChain,ask=FALSE)
  dev.off()
  
  #get mean
  mean = apply( finalChain , 2 , mean )
  
  #get 95% HPD interval
  hpd = HPDinterval(finalChain, prob=0.95)
  
  #get ESS
  ess = effectiveSize(finalChain)
  
  #print parameter summary and save to file
  parmSummary = cbind(mean, hpd, ess)
  print(parmSummary)
  parmSummary = cbind(rownames(parmSummary), parmSummary)
  colnames(parmSummary) = c("parameter", "mean", "lower", "upper", "ess")
  parmSummary = as.data.frame(parmSummary, row.names=1:nrow(parmSummary))
  parmFile = paste(pathRoot, "inference/", seed, "/parm_summary.csv", sep="")
  write.csv(parmSummary, parmFile, row.names=F)
  
  #analyse infection times
  infTimeLog = paste(pathRoot, "inference/", seed, "/chain_inf_times.txt", sep="")
  inf.infTimes = read.table(infTimeLog, header=T)
  #remove burn in and thin
  inf.infTimes = mcmc(inf.infTimes[seq(1, nrow(inf.infTimes), by = thin.factor),], 
                      start=burnIn, end=nrow(inf.infTimes), thin=thin.factor)
  
  #mean infection times
  inf.meanInfTimes = apply( inf.infTimes , 2 , mean )
  inf.meanInfTimes = inf.meanInfTimes +1 #convert back to numbering from 1 rather than zero
  
  #get ESS for infection times
  ess.infTimes = effectiveSize(inf.infTimes)
  #get HPD
  inf.infTimes.hpd = HPDinterval(inf.infTimes)  +1 #convert back to numbering from 1 rather than zero
  infTimes.summary = cbind(inf.meanInfTimes, inf.infTimes.hpd)
  
  #analyse recovery times
  recTimeLog = paste(pathRoot, "inference/", seed, "/chain_rec_times.txt", sep="")
  rec.recTimes = read.table(recTimeLog, header=T)
  # remove burn in and thin
  rec.recTimes = mcmc(rec.recTimes[seq(1, nrow(rec.recTimes), by = thin.factor),], 
                      start=burnIn, end=nrow(rec.recTimes), thin=thin.factor)
  
  rec.meanrecTimes = apply( rec.recTimes , 2 , mean )
  rec.meanrecTimes = rec.meanrecTimes +1 #convert back to numbering from 1 rather than zero
  
  #get ESS for recection times
  ess.recTimes = effectiveSize(rec.recTimes)
  #get HPD
  rec.recTimes.hpd = HPDinterval(rec.recTimes)  +1 #convert back to numbering from 1 rather than zero
  recTimes.summary = cbind(rec.meanrecTimes, rec.recTimes.hpd)
  
  ##plot ESS and differences for infection times and recovery times
  infRecFile = paste(pathRoot, "inference/", seed, "/inf_rec_plots.pdf", sep="")
  pdf(infRecFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
  par(mfrow = c(2,2))
  hist(ess.infTimes)
  hist(sampleTimes-inf.meanInfTimes)
  hist(ess.recTimes)
  hist(rec.meanrecTimes-sampleTimes)
  dev.off()
  
  
  
  #### source type analysis ####
  infSourceTypesLog = paste(pathRoot, "inference/", seed, "/chain_inf_source_types.txt", sep="")
  all.inf.infSourceTypes = read.table(infSourceTypesLog, header=T)
  rowN = nrow(all.inf.infSourceTypes)
  #remove burn in and thin
  inf.infSourceTypes = all.inf.infSourceTypes[burnIn:rowN,]
  inf.infSourceTypes = inf.infSourceTypes[seq(1, nrow(inf.infSourceTypes), by = thin.factor),]
  
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
  srcTypeFile = paste(pathRoot, "inference/", seed, "/source_type_summary.csv", sep="")
  write.csv(summary.srctype, srcTypeFile, row.names=F)
  
  
  
  ### EXACT SOURCE MATCHES
  infSourceLog = paste(pathRoot, "inference/", seed, "/chain_inf_sources.txt", sep="")
  all.inf.infSources = read.table(infSourceLog, header=T)
  rowN = nrow(all.inf.infSources)
  inf.infSources = all.inf.infSources[burnIn:rowN,]
  inf.infSources = inf.infSources[seq(1, nrow(inf.infSources), by = thin.factor),]
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
  
  srcFile = paste(pathRoot, "inference/", seed, "/source_summary.csv", sep="")
  write.csv(infSrc.summary, srcFile, row.names=F)
  
  
  
  
  
  
  
  
  ## create plots of infection sources and times
  
  
  ## functions for ward overlap plotting
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
    
    p = ggplot(data=df.ward.sample) +
      geom_errorbarh(mapping=aes(y=factor(ward), x=t_admit, xmin=t_admit, xmax=t_discharge, 
                                 color=factor(patient_id)), height=0.4, size=1, alpha=0.5) +
      geom_point(mapping=aes(y=factor(ward), x=t_sample, color=factor(patient_id)), size=3, shape=4, stroke=1.5) +
      labs(y="Ward", x="Time",
           title=paste("Ward overlaps")) +
      scale_colour_discrete(name="Patient")
    return(p)
  }
  
  plotWardWrapper = function(ptString, wardLog, patientLog, inf.infSources, inf.infTimes) {
    expand = 100
    pt = which(patientLog$patient_id==ptString)
    sources = table(inf.infSources[,pt])
    ptList = c(ptString, names(sources))
    inf.times = inf.infTimes[,pt] +1
    wardMinMax = c(min(inf.times), max(inf.times))
    p = plotWard(ptList, wardLog, wardMinMax, expand, patientLog)
    return(p)
  }
  
  
  ## start loop through all infection source and time plotting
  
  
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
            axis.text=element_text(size=10), axis.text.x=element_text(angle=45,hjust=1),
            axis.title=element_text(size=10), plot.title=element_text(size=10))
    
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
    
    dLen = dMax-dMin
    x.ticks = seq(dMin, dMax, ceiling(dLen/10))
    
    plot[[i]] = ggplot(d, aes(x=InfectionTimes, y=Probability,
                              fill=factor(SourceType, levels=c(0, 1,2,3,4,5),
                                          labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore")))) +
      geom_bar(stat="identity") +
      scale_y_continuous(limits=c(0,1)) +
      scale_x_discrete(limits=x.ticks) +
      labs(y="Posterior probability", x="Estimated infection time",
           title=paste("Infection times for patient ", patientLog$patient_id[pt], sep="")) +
      scale_fill_manual(values=colors, drop=F, name="Source Type") +
      theme(legend.text=element_text(size=10), legend.title=element_text(size=10),
            axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10))
    i = i+1
    
    plot[[i]] = plotWardWrapper(ptString = patientLog$patient_id[pt], wardLog, patientLog, inf.infSources, inf.infTimes)
    
    i = i+1
  }
  
  #save plots to pdf - 4 per page
  srcPlotFile = paste(pathRoot, "inference/", seed, "/source_plots.pdf", sep="")
  glist = lapply(plot, ggplotGrob)
  ggsave(srcPlotFile, marrangeGrob(glist, layout_matrix=matrix(c(1,2,3,3), 2), nrow=2, ncol=2), width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
  
  
  #create a graph of the mode sources
  library(igraph)
  g = make_empty_graph() 
  # add all the cases
  g = g + vertices(as.character(infSrc.summary$recipient), size=5, label.cex=0.5)
  
  #add edges for genetic links
  geneticList = melt(cbind(rownames(geneticDist), geneticDist))
  geneticMatches = geneticList[which(geneticList$value<=2),]
  for (i in 1:nrow(geneticMatches)) {
    if (gsub("C","",geneticMatches[i,1]) < gsub("C","",geneticMatches[i,2])) {
      g = g + edge(geneticMatches[i,1], geneticMatches[i,2], 
                   color="grey", size=2, arrow.width=0.1, arrow.size=0.1)
    }
  }
  
  
  #add edges for transmitted cases
  for (i in 1:nrow(infSrc.summary)) {
    if (infSrc.summary$mode_source[i]!=-1) {
      g = g + edge(as.character(infSrc.summary$mode_source[i]), as.character(infSrc.summary$recipient[i]), 
                   color="red", size=2, arrow.width=1, arrow.size=10)
    }
  }
  
  #graph summary
  pdf(file=paste(pathRoot, "inference/", seed, "/transmission_snp_network.pdf", sep=""),  width=297/25.4, height=210/25.4)
  plot(g, layout=layout_with_kk(g, epsilon=0, maxiter=50000, kkconst=200))
  dev.off()
 
}

#allow path to be hard coded, but also allow this to be changed at run time
pathRoot = "/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenarios/simulation_22789326"
burnIn.factor = 0.2
thin.factor = 10

#parse command line options 
option_list = list(
  make_option( c("-d", "--directory"), type="character", default=pathRoot, 
               metavar="character"),
  make_option( c("-b", "--burnin"), type="double", default=burnIn.factor, 
               help="burn in"),
  make_option( c("-t", "--thin"), type="double", default=thin.factor, 
               help="thin chains by saving 1 in x iterations")
  )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pathRoot = paste(opt$directory, "/", sep="")
burnIn.factor = opt$burnin
thin.factor = opt$thin

simId = tail(strsplit(pathRoot, '/')[[1]], n=1)
seedList = list.dirs(path = paste(pathRoot,"/inference", sep=""), full.names = F, recursive = F)
seedList = grep("seed*", seedList, value = TRUE) #require to include seed to allow to ignore some chains
seedList = seedList[(which(seedList!="seed_combined"))]


#create folder with merged data
dir.create(paste(pathRoot, "inference/seed_combined", sep=""))

chainList = mcmc.list()
chainListParm = mcmc.list()

i=1
for (seed in seedList) {
  mcmcLog = paste(pathRoot, "inference/", seed, "/chain_parameters.txt", sep="")
  infTimeLog = paste(pathRoot, "inference/", seed, "/chain_inf_times.txt", sep="")
  recTimeLog = paste(pathRoot, "inference/", seed, "/chain_rec_times.txt", sep="")
  infSourceLog = paste(pathRoot, "inference/", seed, "/chain_inf_sources.txt", sep="")
  infSourceTypesLog = paste(pathRoot, "inference/", seed, "/chain_inf_source_types.txt", sep="")
  
  chain = read.table(mcmcLog, header=T)
  all.inf.infTimes = read.table(infTimeLog, header=T)
  all.rec.recTimes = read.table(recTimeLog, header=T)
  all.inf.infSources = read.table(infSourceLog, header=T)
  all.inf.infSourceTypes = read.table(infSourceTypesLog, header=T)
  rowN = nrow(chain)
  burnIn = (rowN * burnIn.factor)+1
  
  finalChain = chain[burnIn:nrow(chain),]
  finalChain = cbind(finalChain, logistic(finalChain$spore_prob_logit), logistic(finalChain$p_start_inf_logit))
  thinChain = finalChain[seq(1, nrow(finalChain), by = thin.factor),]
  infTimes = all.inf.infTimes[burnIn:rowN,][seq(1, rowN-burnIn+1, by = thin.factor),]
  recTimes = all.inf.infTimes[burnIn:rowN,][seq(1, rowN-burnIn+1, by = thin.factor),]
  infSources = all.inf.infSources[burnIn:rowN,][seq(1, rowN-burnIn+1, by = thin.factor),]
  infSourceTypes = all.inf.infSourceTypes[burnIn:rowN,][seq(1, rowN-burnIn+1, by = thin.factor),]
  
  chainList[[i]] = mcmc(thinChain, thin=thin.factor)
  chainListParm[[i]] = mcmc(thinChain, thin=thin.factor)[,1:14]
  
  if (i==1) {
    mergedChain = thinChain
    mergedInfTimes = infTimes
    mergedRecTimes = recTimes
    mergedInfSources = infSources
    mergedInfSourceTypes = infSourceTypes
  } else {
    mergedChain = rbind(mergedChain, thinChain)
    mergedInfTimes = rbind(mergedInfTimes, infTimes)
    mergedRecTimes = rbind(mergedRecTimes, recTimes)
    mergedInfSources = rbind(mergedInfSources, infSources)
    mergedInfSourceTypes = rbind(mergedInfSourceTypes, infSourceTypes)
  }
  i = i +1
}

#save chains
write.table(mergedChain[,1:19], paste(pathRoot, "inference/seed_combined/chain_parameters.txt", sep=""), row.names = F, sep="\t")
write.table(mergedInfTimes, paste(pathRoot, "inference/seed_combined/chain_inf_times.txt", sep=""), row.names = F, sep="\t")
write.table(mergedRecTimes, paste(pathRoot, "inference/seed_combined/chain_rec_times.txt", sep=""), row.names = F, sep="\t")
write.table(mergedInfSources, paste(pathRoot, "inference/seed_combined/chain_inf_sources.txt", sep=""), row.names = F, sep="\t")
write.table(mergedInfSourceTypes, paste(pathRoot, "inference/seed_combined/chain_inf_source_types.txt", sep=""), row.names = F, sep="\t")

#save converagance plots for parameters, posterior and likelihoods
parmChainFile = paste(pathRoot, "inference/seed_combined/covergence_plots.pdf", sep="")
pdf(parmChainFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
plot(chainList)
gelman.plot(chainList)
dev.off()

#save converagance plots for parameters only
parmChainFile = paste(pathRoot, "inference/seed_combined/covergence_plots_parm.pdf", sep="")
pdf(parmChainFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
plot(chainListParm)
gelman.plot(chainListParm)
dev.off()

#save gelman diagonal to file
gelman = gelman.diag(chainListParm)
capture.output(gelman, file = paste(pathRoot, "inference/seed_combined/covergence_factors.txt", sep=""))

#run reporting on merged data
runReport(pathRoot, seed="seed_combined", burnIn.factor=0, thin.factor=1)

#report each seed separately
for (seed in seedList) {
  runReport(pathRoot, seed, burnIn.factor, thin.factor)
}

