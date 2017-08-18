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
path = "/Users/davideyre/Dropbox/Transmission_Inference/xcode_project/sim_data/simulation_23574379"

#parse command line options - more here - https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
option_list = list(
  make_option( c("-d", "--directory"), type="character", default=path, 
               metavar="character") )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
path = paste(opt$directory, "/", sep="")

simId = strsplit(opt$directory, '/')[[1]]
simId = simId[length(simId)]


## read in true values
patientLog = read.csv(file = paste(path, "input/patientLog.csv", sep=""))
wardLog = read.csv(file = paste(path, "input/wardLog.csv", sep=""))
geneticDist = read.table(file=paste(path, "input/simDistances.txt", sep=""))

#reformat simulation log for ease of reading code
true.infTimes = patientLog$t_inf
true.infSources = cbind(patientLog$source, patientLog$source_type)
sampleTimes = patientLog$t_sample
true.recTimes = patientLog$t_recover
true.sporeDurations = patientLog$spore_duration
infectedPatients = which(!is.na(true.infTimes))
n = length(infectedPatients)
maxTime = max(wardLog$t_discharge)
nWards = max(wardLog$ward)
nPatients = max(patientLog$patient_id)

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

#checks for logistic transform of spore.p
# logistic = function(x) {return(1/(1+exp(-x)))}
# plot(as.mcmc(logistic(finalChain$spore_prob)))
# plot(as.mcmc((finalChain$beta1)))
# mean(as.mcmc(logistic(finalChain$spore_prob)))
# mean(as.mcmc((finalChain$beta1)))
# HPDinterval(as.mcmc(logistic(finalChain$spore_prob)))
# HPDinterval(as.mcmc((finalChain$beta1)))
# effectiveSize(as.mcmc(logistic(finalChain$spore_prob)))
# effectiveSize(as.mcmc((finalChain$beta1)))



#print parameter summary and save to file
parmSummary = cbind(mean, hpd, ess)
print(parmSummary)
parmSummary = cbind(rep(simId, nrow(parmSummary)), rownames(parmSummary), parmSummary)
colnames(parmSummary) = c("simulation", "parameter", "mean", "lower", "upper", "ess")
parmSummary = as.data.frame(parmSummary, row.names=1:nrow(parmSummary))
parmFile = paste(path, "inference/parm_summary.csv", sep="")
write.csv(parmSummary, parmFile, row.names=F)


#epsilon check
#getEll = function(lambda, dd) {
#  ll = sum(dpois(dd, lambda, log=T))
#  return(ll)
#}
#dd = sampleTimes[infectedPatients]-true.infTimes[infectedPatients]
#o.epsilon = optimise(getEll, c(0.1,10), dd, maximum = T)
#print("ML estimate of epsilon")
#print(o.epsilon$maximum)


#pairs(data.frame(cbind(finalChain$spore_prob, finalChain$beta1)))


#compare infection times to true infection times
infTimeLog = paste(path, "inference/chain_inf_times.txt", sep="")
inf.infTimes = read.table(infTimeLog, header=T)
inf.infTimes = inf.infTimes[burnIn:nrow(inf.infTimes),]
inf.meanInfTimes = apply( inf.infTimes , 2 , mean )
inf.meanInfTimes = inf.meanInfTimes +1 #convert back to numbering from 1 rather than zero

#get ESS for infection times
ess.infTimes = effectiveSize(as.mcmc(inf.infTimes))


inf.infTimes.hpd = HPDinterval(as.mcmc(inf.infTimes))  +1 #convert back to numbering from 1 rather than zero
inf.infTimes.hpd.match = true.infTimes[infectedPatients] >= inf.infTimes.hpd[,1] & true.infTimes[infectedPatients] <= inf.infTimes.hpd[,2]

inf.diffInfTimes = inf.meanInfTimes - true.infTimes[infectedPatients]
infTimes.compare = cbind(inf.meanInfTimes, true.infTimes[infectedPatients], inf.diffInfTimes, inf.infTimes.hpd, inf.infTimes.hpd.match)
#print(infTimes.compare)


#inf times with exact match
n.infTime.exact = sum(round(inf.diffInfTimes)==0)
# inf times with hpd match
n.infTime.hpd = sum(infTimes.compare[,6])

# create transmission summary 
transSummary = c("Infection Times", n, n.infTime.exact, n.infTime.exact/n, n.infTime.hpd, n.infTime.hpd/n)
names(transSummary) = c("parameter", "n", "n_exact_match", "proportion_exact_match", "n_hpd_match", "proportion_hpd_match")


#compare recovery times to true recovery times
recTimeLog = paste(path, "inference/chain_rec_times.txt", sep="")
rec.recTimes = read.table(recTimeLog, header=T)
rec.recTimes = rec.recTimes[burnIn:nrow(rec.recTimes),]
rec.meanrecTimes = apply( rec.recTimes , 2 , mean )
rec.meanrecTimes = rec.meanrecTimes +1 #convert back to numbering from 1 rather than zero

#get ESS for recection times
ess.recTimes = effectiveSize(as.mcmc(rec.recTimes))


rec.recTimes.hpd = HPDinterval(as.mcmc(rec.recTimes))  +1 #convert back to numbering from 1 rather than zero
rec.recTimes.hpd.match = true.recTimes[infectedPatients] >= rec.recTimes.hpd[,1] & true.recTimes[infectedPatients] <= rec.recTimes.hpd[,2]

rec.diffrecTimes = rec.meanrecTimes - true.recTimes[infectedPatients]
recTimes.compare = cbind(rec.meanrecTimes, true.recTimes[infectedPatients], rec.diffrecTimes, rec.recTimes.hpd, rec.recTimes.hpd.match)
#print(recTimes.compare)


#rec times with exact match
n.recTime.exact = sum(round(rec.diffrecTimes)==0)
# rec times with hpd match
n.recTime.hpd = sum(recTimes.compare[,6])

# create transmission summary 
transSummary = rbind(transSummary, c("Recovery Times", n, n.recTime.exact, n.recTime.exact/n, n.recTime.hpd, n.recTime.hpd/n))






##plot ESS and differences for infection times and recovery times

infRecFile = paste(path, "inference/inf_rec_plots.pdf", sep="")
pdf(infRecFile, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)
par(mfrow = c(2,2))
hist(ess.infTimes)
hist(inf.diffInfTimes)
hist(ess.recTimes)
hist(rec.diffrecTimes)
dev.off()










#### compare source types to true source types ####
infSourceTypesLog = paste(path, "inference//chain_inf_source_types.txt", sep="")
all.inf.infSourceTypes = read.table(infSourceTypesLog, header=T)
rowN = nrow(all.inf.infSourceTypes)
inf.infSourceTypes = all.inf.infSourceTypes[burnIn:rowN,]

##ESS for infection source types
#print ( effectiveSize(as.mcmc(inf.infSourceTypes)) )

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
true.infSources = true.infSources[infectedPatients,]

#MAP esitmate
infSourceType.map = max.col(a)-1
infSourceType.match = (infSourceType.map == true.infSources[,2])
infSourceType.compare = cbind(a, true.infSources[,2], infSourceType.match)
#print(infSourceType.compare)



#HPD estimate
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

infSourceType.hpd = rep(0, ncol(inf.infSourceTypes))
for (i in 1:ncol(inf.infSourceTypes)) {
  if( true.infSources[i,2] %in% typeHPD(a[i,], rowCount) ) {
    infSourceType.hpd[i] = 1
  }
}


#inf source types with exact match
n.srcType.exact = sum(infSourceType.compare[,8])
# inf source types with hpd match
n.srcType.hpd = sum(infSourceType.hpd)

# create transmission summary 
transSummary = rbind(transSummary, c("Source Types", n, n.srcType.exact, n.srcType.exact/n, n.srcType.hpd, n.srcType.hpd/n))



### EXACT SOURCE MATCHES
infSourceLog = paste(path, "inference/chain_inf_sources.txt", sep="")
all.inf.infSources = read.table(infSourceLog, header=T)
rowN = nrow(all.inf.infSources)
inf.infSources = all.inf.infSources[burnIn:rowN,]

inf.infSources.mode = apply(inf.infSources, 2, Mode)

infSources.compare = cbind(true.infSources[,1], inf.infSources.mode, inf.infSources.mode==true.infSources[,1])



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

infSources.hpd = rep(0, ncol(inf.infSources))
for (i in 1:ncol(inf.infSources)) {
  if( true.infSources[i,1] %in% sourceHPD(inf.infSources[,i] ,rowCount) ) {
    infSources.hpd[i] = 1
  }
}


cbind(infSources.compare,infSources.hpd)


#inf source types with exact match
n.src.exact = sum(infSources.compare[,3])
# inf source types with hpd match
n.src.hpd = sum(infSources.hpd)

# create transmission summary 
transSummary = rbind(transSummary, c("Sources", n, n.src.exact, n.src.exact/n, n.src.hpd, n.src.hpd/n))
transSummary = cbind(rep(simId, nrow(transSummary)), transSummary)
colnames(transSummary) = c("simulation", "parameter","n","n_exact_match","proportion_exact_match","n_hpd_match","proportion_hpd_match")


#write transmission summary to file
transSummary = as.data.frame(transSummary, row.names=1:nrow(transSummary))
transSummaryFile = paste(path, "inference/trans_summary.csv", sep="")
write.csv(transSummary, transSummaryFile, row.names=F)
print(transSummary)
##export summary of key parameters





##THINGS GET WRONG
print("Patients with no match in source HPD")
errorPt = which(infSources.hpd==0)
if(length(errorPt)>1) {
  errorTable = cbind(inf.infSources.mode[errorPt], infSourceType.map[errorPt], true.infSources[errorPt,])
  colnames(errorTable) = c("src_mode", "srcType_mode", "src_true", "srcType_true")
  print(errorTable)
} else if (length(errorPt)==1) {
  errorTable = c(inf.infSources.mode[errorPt], infSourceType.map[errorPt], true.infSources[errorPt,])
  names(errorTable) = c("src_mode", "srcType_mode", "src_true", "srcType_true")
  print(errorTable)
}



ixPatient = function(pt) {
  ptIndex = which(infectedPatients==208)
  print(paste("t_inf:", true.infTimes[pt]))
  print(paste("t_sample:", sampleTimes[pt]))
  print(paste("true_src:", true.infSources[ptIndex,1]))
  print(paste("true_srcType:", true.infSources[ptIndex,2]))
  print("admissions for patient:")
  print(wardLog[which(wardLog[,1]==pt),])
  print("--------")
  
  print("admissons for source")
  src = true.infSources[ptIndex,1]
  srcIndex = which(infectedPatients==src)
  print(wardLog[which(wardLog[,1]==src),])
  print(paste("src_t_inf:", true.infTimes[src]))
  print(paste("src_t_sample:", sampleTimes[src]))
  print("--------")
  
  print("inferred sources:")
  print(table(inf.infSources[,paste("patient_", pt-1, sep="")]))
  print("inferred source types:")
  print(table(inf.infSourceTypes[,paste("patient_", pt-1, sep="")]))
  
}

#call using patient numbering starting from one
#ixPatient(208)

#ix
# table(inf.infSources[,"patient_2"]+1)
# table(inf.infSourceTypes[,"patient_2"])
cbind(infectedPatients, infSourceType.compare)[which(infSourceType.hpd==0),] 
#right vs wrong, by true source type

##table of src type results with HPD match or not

## table of accuracy of source type prediction, by true infection source type
srcComparison = matrix(NA, nrow=3, ncol=6)
colnames(srcComparison) = c("hosp_bground", "ward", "hosp", "comm_bground", "start_inf", "spore")
rownames(srcComparison) = c("true", "hpd_match", "hpd_no_match")
for (i in 1:3) {
  for (j in 1:6) {
    srcIndex = j-1
    if (i==1) {
      #true answer
      srcComparison[i,j] = sum(true.infSources==srcIndex)
    } else if (i==2) {
      srcComparison[i,j] = sum(infSourceType.compare[which(infSourceType.hpd==1),7]==srcIndex)
    } else {
      srcComparison[i,j] = sum(infSourceType.compare[which(infSourceType.hpd==0),7]==srcIndex)
    }
  }
}
print(t(srcComparison))


#compare true (simulated) and rec/inf times from MCMC
#scatter.smooth((true.recTimes-true.infTimes)[infectedPatients], rec.meanrecTimes-inf.meanInfTimes)
#tend to underestimate slow recoveries, as little information in data?



#compare recovery times - recTimes.compare[which(rec.recTimes.hpd.match==0),]
# check ward log - wardLog[which(wardLog[,1]==289),]


## check estimate of Ne

# geneticDist = read.table(file=paste(path, "input/simDistances_snps.txt", sep=""))
# #get all direct transmission pairs
# trans = cbind(infectedPatients, true.infSources[infectedPatients,1])
# trans = trans[which(trans[,2]>-1),]
# t.diff = rep(NA, nrow(trans))
# snps = rep(NA, nrow(trans))
# 
# for (i in 1:nrow(trans)) {
#   snps[i] = geneticDist[paste("patient_", trans[i,1], sep=""),paste("patient_", trans[i,2], sep="")]
#   t.diff[i] = abs(sampleTimes[trans[i,1]]-sampleTimes[trans[i,2]])
# }
# 
# 
# igamma = function(a,z) pgamma(z,a,lower=FALSE)*gamma(a)
# # Likelihood for an individual observation
# liki = Vectorize(function(snp,time,theta,mu) {
#   Ne = theta/2/mu
#   if(time==0 | mu==0) return(theta^snp/(theta+1)^(snp+1))
#   return(    exp(time/2/Ne+log(igamma(1+snp,mu*time+time/2/Ne)))     *    (2*mu)^snp   /  factorial(snp)/Ne/(2*mu+1/Ne)^(snp+1))
# })
# # Generate a likelihood across all observations with parameter vector log(theta,mu)
# get_loglik = function(snp,time) {
#   loglik = function(param) sum(log(liki(snp,time,exp(param[1]),exp(param[2]))))
# }
# # Generate likelihood function using all observations
# loglik = get_loglik(snps,t.diff)
# init_param = log(c(0.5,0.5))
# opt = optim(init_param,loglik,control=list(fnscale=-1))
# # Parameter estimates
# (theta_hat = exp(opt$par[1]))
# (mu_hat = exp(opt$par[2]))
# (Ne_hat = theta_hat/2/mu_hat)


getRoute = function(srcType) {
  o = factor(srcType, levels=c(0, 1,2,3,4,5), 
             labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore"))
  return(as.character(o))
}

colors = brewer.pal(6,"Spectral")

#save plots as list
plot=list()
i=1
for (pt.R in infectedPatients) {
  pt = paste("patient_", pt.R, sep="")
  pt.index = which(infectedPatients==pt.R)
  d = table(inf.infSources[,pt], inf.infSourceTypes[,pt])
  d = as.data.frame(cbind(as.numeric(rownames(d)), d))
  colnames(d)[1] = "source"
  rownames(d) = 1:nrow(d)
  d = melt(d, id=c("source"))
  d = cbind(d, d[,3]/sum(d[,3]))
  colnames(d) = c("Source", "SourceType", "Frequency", "Probability")
  d = d[which(d$Frequency>0),]
  plot[[i]] = ggplot(d, aes(x=factor(Source), y=Probability, 
                            fill=factor(SourceType, levels=c(0, 1,2,3,4,5), 
                                        labels=c("Hospital background", "Ward", "Hospital-wide", "Community background", "Start positive", "Spore")))) +   
    geom_bar(stat="identity") +
    scale_y_continuous(limits=c(0,1)) + 
    labs(y="Posterior probability", x="Source", 
         title=paste("Infection sources for patient ", pt.R, 
                     " (true source: ", true.infSources[pt.index,1], 
                     " via route: ", getRoute(true.infSources[pt.index,2]), ")", sep="")) +
    scale_fill_manual(values=colors, drop=F, name="Source Type") +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=10), 
          axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10)) 
  i=i+1
  
  #plot difference in inferred infection time and true time
  inf.diffs = inf.infTimes[,pt.index] - (true.infTimes[infectedPatients[pt.index]] - 1)
  d = table(inf.diffs, inf.infSourceTypes[,pt])
  dMin = min(inf.diffs)
  dMax = max(inf.diffs)
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
    labs(y="Posterior probability", x="Estimated - true infection time", 
         title=paste("Infection times for patient ", pt.R, sep="")) +
    scale_fill_manual(values=colors, drop=F, name="Source Type") +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=10), 
          axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10)) +
    scale_x_continuous(breaks=dMin:dMax)
  
  # d = as.data.frame(inf.infTimes[,pt.index] - true.infTimes[infectedPatients[pt.index]])
  # dMin = min(d)
  # dMax = max(d)
  # colnames(d) = c("infTimes")
  # plot[[i]] = ggplot(d, aes(x=infTimes)) + 
  #   geom_histogram(binwidth = 1, col="orange", fill="orange", alpha=0.7) +
  #   labs(x="Estimated - true infection time", y="Frequency", title=paste("Infection times for patient ", pt.R, sep="")) +
  #   theme_bw() +
  #   theme(axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10)) +
  #   scale_x_continuous(breaks=dMin:dMax)
  
  i = i+1
}

#save plots to pdf - 4 per page
srcPlotFile = paste(path, "inference/source_plots.pdf", sep="")
glist = lapply(plot, ggplotGrob)
ggsave(srcPlotFile, marrangeGrob(glist, nrow = 2, ncol = 2, top=""), width=29.7/2.54, height=21/2.54, useDingbats=FALSE)



#generate comparison of true sources
factoriseSrc = function(srcList) {
  o = factor(srcList, levels=c(0, 1,5,2,3,4), 
         labels=c("Hospital background", "Ward", "Spore", "Hospital-wide", "Community background", "Start positive"))
  return(o)
}

src.compare.true = factoriseSrc(true.infSources[,2])
src.compare.map = factoriseSrc(infSourceType.map)
d = matrix(NA, ncol=6, nrow=nrow(inf.infSourceTypes))
for (i in 1:nrow(inf.infSourceTypes)) {
  d[i,] = table(factoriseSrc(inf.infSourceTypes[i,]))
}



summary.srctype = cbind(table(src.compare.true), table(src.compare.map), apply( d , 2 , mean ), HPDinterval(as.mcmc(d), prob=0.95))
summary.srctype = cbind(rep(simId, nrow(summary.srctype)), rownames(summary.srctype), summary.srctype)
colnames(summary.srctype) = c("simulation", "source_type", "simulated_counts", "map", "mean", "lower", "upper")
srcTypeFile = paste(path, "inference/source_type_summary.csv", sep="")
write.csv(summary.srctype, srcTypeFile, row.names=F)
