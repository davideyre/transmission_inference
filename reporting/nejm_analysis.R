rm(list = ls())
library(ggplot2)
library(grid)
library(gridExtra)

logistic = function(x) {return(1/(1+exp(-x)))}

setwd("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/nejm/")


parms = c("beta0", "beta1", "beta2", "betacomm", 
          "logistic(finalChain$spore_prob_logit)",
          "spore_multiplier", "exp(finalChain$sampleSize)",
          "sampleMu", "rec_size", "rec_mu", "directNe", "introNe", "mu")

prior = matrix(NA, ncol=4, nrow=length(parms))
prior[1,] = c(parms[1], qgamma(c(0.5,0.025,0.975), shape=2, rate=1/0.002))
prior[2,] = c(parms[2], qgamma(c(0.5,0.025,0.975), shape=2, rate=1/0.002))
prior[3,] = c(parms[3], qgamma(c(0.5,0.025,0.975), shape=2, rate=1/0.002))
prior[4,] = c(parms[4], qgamma(c(0.5,0.025,0.975), shape=2, rate=1/0.002))

prior[5,] = c(parms[5], logistic(qnorm(c(0.5,0.025,0.975), mean=0, sd=1.7)))
prior[6,] = c(parms[6], logistic(qnorm(c(0.5,0.025,0.975), mean=0, sd=1.7)))

prior[7,] = c(parms[7], exp(qnorm(c(0.5,0.025,0.975), mean=0.7, sd=0.5)))
prior[8,] = c(parms[8], qgamma(c(0.5,0.025,0.975), shape=3, rate=0.1))

prior[9,] = c(parms[9], qnorm(c(0.5,0.025,0.975), mean=3, sd=0.5))
prior[10,] = c(parms[10], qnorm(c(0.5,0.025,0.975), mean=90, sd=3))

prior[11,] = c(parms[11], qnorm(c(0.5,0.025,0.975), mean=22.5, sd=0.1))
prior[12,] = c(parms[12], qgamma(c(0.5,0.025,0.975), shape=2, rate=1/10000))
prior[13,] = c(parms[13], qnorm(c(0.5,0.025,0.975), mean=1/365.25, sd=0.05/365.25))

colnames(prior) = c("parameter", "mean", "lower", "upper")
df.prior = as.data.frame(cbind(st="Prior", prior, ess=NA))

#based on simulation from prior - sampling delay
set.seed(42)
sampling.shape = exp(rnorm(100000, mean=0.7, sd=0.5))
sampling.mu = rgamma(100000, shape=3, rate=0.1)
sampling.rate = sampling.shape/sampling.mu
df.sampling = as.data.frame(rgamma(100000, shape=sampling.shape, rate=sampling.rate))
colnames(df.sampling) = c("sampling")

#based on simulation from prior - recovery time
recovery.shape = rnorm(100000, mean=3, sd=0.5)
recovery.mu = rnorm(100000, mean=90, sd=3)
recovery.rate = recovery.shape/recovery.mu
df.recovery = as.data.frame(rgamma(100000, shape=recovery.shape, rate=recovery.rate))
colnames(df.recovery) = c("recovery")

#spore persistence
spore.p = logistic(rnorm(100000, mean=0, sd=1.7))
df.spore = as.data.frame(rgeom(100000, spore.p))
colnames(df.spore) = c("spore")

p.s = ggplot(df.sampling, aes(x=sampling)) +
  geom_density(alpha=.65, fill="red") +
  labs(y="Density", x="Time from infection to first postive test", title="Prior distribution: Sampling")+
  xlim(c(0,400))
  
p.r = ggplot(df.recovery, aes(x=recovery)) +
  geom_density(alpha=.65, fill="lightblue") +
  labs(y="Density", x="Time from first postive test to recovery", title="Prior distribution: Recovery") +
  xlim(c(0,400))

p.spore = ggplot(df.spore, aes(x=spore)) +
  geom_density(alpha=.65, fill="green") +
  labs(y="Density", x="Time from ward contamination to spore removal", title="Prior distribution: Spores") +
  xlim(c(0,250))


#plot posterior distributions for ST1
df.st1.sampling = read.csv("merged_summary/st1_sampling_intervals.csv")
st1.samplediff = data.frame(sampling=as.vector(as.matrix(df.st1.sampling)))
st1.s = ggplot(st1.samplediff, aes(x=sampling)) +
  geom_histogram(aes(y=..density..), bins=60, fill="red") +
  labs(y="Density", x="Time from infection to first postive test", title="Posterior distribution, ST1: Sampling")+
  xlim(c(0,400)) +
  ylim(c(0,0.03))

df.st1.rec = read.csv("merged_summary/st1_recovery_intervals.csv")
st1.recdiff = data.frame(rec=as.vector(as.matrix(df.st1.rec)))
st1.r = ggplot(st1.recdiff, aes(x=rec)) +
  geom_histogram(aes(y=..density..), bins=60, fill="lightblue") +
  labs(y="Density", x="Time from first postive test to recovery", title="Posterior distribution, ST1: Recovery")+
  xlim(c(0,400)) +
  ylim(c(0,0.01))

df.st1.spore.p = read.csv("merged_summary/st1_spore_p.csv")
df.st1.spore = data.frame(spore=rgeom(100000, df.st1.spore.p$spore.p))
st1.spore = ggplot(df.st1.spore, aes(x=spore)) +
  geom_histogram(aes(y=..density..), bins=60, fill="green") +
  labs(y="Density", x="Time from ward contamination to spore removal", title="Posterior distribution, ST1: Spores") +
  xlim(c(0,250))

#plot posterior distributions for ST2
df.st2.sampling = read.csv("merged_summary/st2_sampling_intervals.csv")
st2.samplediff = data.frame(sampling=as.vector(as.matrix(df.st2.sampling)))
st2.s = ggplot(st2.samplediff, aes(x=sampling)) +
  geom_histogram(aes(y=..density..), bins=60, fill="red") +
  labs(y="Density", x="Time from infection to first postive test", title="Posterior distribution, ST2: Sampling")+
  xlim(c(0,400)) +
  ylim(c(0,0.03))

df.st2.rec = read.csv("merged_summary/st2_recovery_intervals.csv")
st2.recdiff = data.frame(rec=as.vector(as.matrix(df.st2.rec)))
st2.r = ggplot(st2.recdiff, aes(x=rec)) +
  geom_histogram(aes(y=..density..), bins=60, fill="lightblue") +
  labs(y="Density", x="Time from first postive test to recovery", title="Posterior distribution, ST2: Recovery")+
  xlim(c(0,400)) +
  ylim(c(0,0.01))

df.st2.spore.p = read.csv("merged_summary/st2_spore_p.csv")
df.st2.spore = data.frame(spore=rgeom(100000, df.st2.spore.p$spore.p))
st2.spore = ggplot(df.st2.spore, aes(x=spore)) +
  geom_histogram(aes(y=..density..), bins=60, fill="green") +
  labs(y="Density", x="Time from ward contamination to spore removal", title="Posterior distribution, ST2: Spores") +
  xlim(c(0,250)) +
  ylim(c(0,0.06))

p.many = arrangeGrob(p.s, p.r, p.spore, st1.s, st1.r, st1.spore, st2.s, st2.r, st2.spore, nrow=3, ncol=3)
ggsave("merged_summary/plot_prior_sampling_recovery_spore.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)



## to do - save sampling intervals and recovery times, by ST
## for spore recovery - save values of spore.p from chain and simulate from these











### BY ST COMPARISONS ###

## Parameters
df = read.csv("merged_summary/parm_compare_by_st.csv")
parmList = unique(df$parameter)
parmList.keep = parmList[c(1,2,3,9,21,12,20,5,13,14,6,7,8)]

df = df[which(df$parameter %in% parmList.keep),]

df = rbind(df, df.prior)
df$mean = as.numeric(df$mean)
df$lower = as.numeric(df$lower)
df$upper = as.numeric(df$upper)

df$parm.lab = factor(df$parameter, levels=c("beta0", "beta1", "beta2", "betacomm", 
                                         "logistic(finalChain$spore_prob_logit)",
                                         "spore_multiplier", "exp(finalChain$sampleSize)",
                                         "sampleMu", "rec_size", "rec_mu", "directNe", "introNe", "mu"),
                  labels=c("beta_hospital_bg", "beta_ward", "beta_hospital", "beta_community",
                           "spore_decay_probability", "spore_beta_multiplier", "sample_delay_dist_size",
                           "sample_delay_dist_mu", "recovery_delay_dist_size", "recovery_delay_dist_mu",
                           "direct_Ne", "intro_Ne", "mutation_rate"))

parm.lab=c(expression(paste("Parameter estimates, by ST: ", beta["hospital background"])), 
           expression(paste("Parameter estimates, by ST: ", beta["ward"])), 
           expression(paste("Parameter estimates, by ST: ", beta["hospital-wide"])), 
           expression(paste("Parameter estimates, by ST: ", beta["community background"])), 
           "Parameter estimates, by ST: spore decay probability", 
           expression(paste("Parameter estimates, by ST: spore ", beta, " multiplier")), 
           "Parameter estimates, by ST: sample delay distribution, size",
           "Parameter estimates, by ST: sample delay distribution, mean",
           "Parameter estimates, by ST: recovery delay distribution, size",
           "Parameter estimates, by ST: recovery delay distribution, mean",
           expression("Parameter estimates, by ST: direct N"[e]),
           expression("Parameter estimates, by ST: introduction N"[e]),
           "Parameter estimates, by ST: mutation rate")

df$st = as.character(df$st)
df$st = sapply(df$st, toupper)
df$st[which(df$st=="ST_OTHER")] = "Other"
df$st[which(df$st=="PRIOR")] = "Prior"
df$st = factor(df$st, levels=c("ST1", "ST2", "ST3", "ST5", "ST6", "ST8", "ST10", "ST11", "ST42", "ST44", "Other", "Prior"))

parmList.plot = levels(df$parm.lab)
plots = list()
i=1
for (parm in parmList.plot) {
  df.p = df[which(df$parm.lab==parm),]
  if(i<5) {
    df.p = df.p[which(df.p$st!="Prior"),]
  }
  p = ggplot(df.p, aes(x=st, y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper, width=0.6)) +
    labs(title = parm.lab[i], 
         x="Sequence Type", y="Mean (95% HPD)")
  plots[[i]] = p
  i = i+1
}

p.many = arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2)
ggsave("merged_summary/plot_parm_beta.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

p.many = arrangeGrob(plots[[5]], plots[[6]], nrow=2, ncol=2)
ggsave("merged_summary/plot_parm_spore.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

p.many = arrangeGrob(plots[[7]], plots[[8]], plots[[9]], plots[[10]], nrow=2, ncol=2)
ggsave("merged_summary/plot_parm_sample_recovery.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

p.many = arrangeGrob(plots[[11]], plots[[12]], plots[[13]], nrow=2, ncol=2)
ggsave("merged_summary/plot_parm_genetic.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)


## Source type

out = read.csv("merged_summary/trans_src_type_compare_by_st.csv", stringsAsFactors = F)
out = out[which(out$source_type!="Start positive"),]
srcList = unique(out$source_type)
out$st[which(out$st!="st_other")] = sapply(out$st, toupper)
out$st[which(out$st=="st_other")] = "Other"
out$st = factor(out$st, levels=c("ST1", "ST2", "ST3", "ST5", "ST6", "ST8", "ST10", "ST11", "ST42", "ST44", "Other"))

plots = list()
i=1
y.limit = c(1, 0.4, 0.4, 0.4, 1)
for (src in srcList) {
  df.p = out[which(out$source_type==src),]
  p = ggplot(df.p, aes(x=st, y=prop_mean)) +
    geom_point() +
    ylim(0, y.limit[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_errorbar(aes(ymin=prop_lower, ymax=prop_upper, width=0.6)) +
    labs(title = src, 
         x="Sequence Type", y="Proportion of Cases Aquired via this Route\nMean (95% HPD)")
  plots[[i]] = p
  i = i+1
}
p.many = grid.arrange(plots[[2]], plots[[3]], plots[[4]], plots[[1]], plots[[5]], nrow=2)
ggsave("merged_summary/plot_route.pdf", p.many, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)


### Over time analysis

#create data frame
df.bd = read.csv("beddays/beddays_quarter.csv", header=T)
df.all = read.csv("merged_summary/trans_quarter_inpt_summary.csv", header=T)
df.whs = read.csv("merged_summary/trans_quarter_whs_summary.csv", header=T)
df.all = cbind(Route="Total\n(including hospital background)\n", df.all)[which(!df.all$quarter %in% c("2007 Q1", "2007 Q2", "2011 Q2")),]
df.whs = cbind(Route="Sampled source with CDI\n(via ward, spore, hospital-wide)\n", df.whs)[which(!df.whs$quarter %in% c("2007 Q1", "2007 Q2", "2011 Q2")),]
df.all = cbind(df.all, beddays=df.bd$beddays, prop_mean = df.all$mean/df.bd$beddays*10000,
               prop_lower = df.all$lower/df.bd$beddays*10000, prop_upper = df.all$upper/df.bd$beddays*10000)
df.whs = cbind(df.whs, beddays=df.bd$beddays, prop_mean = df.whs$mean/df.bd$beddays*10000,
               prop_lower = df.whs$lower/df.bd$beddays*10000, prop_upper = df.whs$upper/df.bd$beddays*10000)

df = rbind(df.all, df.whs)

p = ggplot(df, aes(x=quarter, y=prop_mean, color=Route, group=Route)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  geom_errorbar(aes(ymin=prop_lower, ymax=prop_upper, width=0.6)) +
  labs(x="Quarter", y="Inpatient aquisitions per 10,000 bed-days\nmean (95% HPD)")
ggsave("merged_summary/plot_over_time.pdf", p, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

### By hospital analysis

## overall by hospital
df.bd = read.csv("beddays/beddays_hosp.csv")
df.all = read.csv("merged_summary/trans_hospital_summary.csv")
df = merge(df.bd, df.all, by="location")
df$spec = "All"

## by hospital limited to specific specialty
df.bdspec = read.csv("beddays/beddays_spec_hosp.csv")
df.spec = read.csv("merged_summary/trans_hospital_spec_summary.csv", stringsAsFactors = F)
df.spec = merge(df.bdspec, df.spec, by="location")
df.spec$location = as.character(df.spec$location)
df.spec$spec = substr(df.spec$location,4,unlist(lapply(df.spec$location, nchar)))
df.spec$location = substr(df.spec$location,1,2)

## merge and label
df = rbind(df, df.spec)
df = cbind(df, prop_mean = df$mean/df$beddays*10000, prop_lower = df$lower/df$beddays*10000,
           prop_upper = df$upper/df$beddays*10000)
df = df[which(df$location!="CO"),]

df$loc.lab = factor(df$location, levels=c("CH", "HG", "JR", "NC"), 
                    labels=c("Churchill", "Horton\nGeneral", "John\nRadcliffe", "Nuffield\nOrthopedic"))
df$spec.lab = factor(df$spec, levels=c("All", "Acute Medicine & Geratology", "General Surgery"))


## plot over all specialities and within AGM and general surgery
df.p = df[which(df$spec %in% c("All", "General Surgery", "Acute Medicine & Geratology")),]



p = ggplot(df.p, aes(x=loc.lab, y=prop_mean, fill=location)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=prop_lower, ymax=prop_upper, width=0.6)) +
  geom_text(lineheight = 0.95, aes(x=loc.lab,y=1.1,label=paste(round(mean,0), "\nCDI\ncases", sep=""))) +
  labs(x="Hospital", y="Inpatient aquisitions per 10,000 bed-days\nmean (95% HPD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  facet_wrap(~spec.lab, nrow=1)
ggsave("merged_summary/plot_by_hospital.pdf", p, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)



### By specialty
df.bd = read.csv("beddays/beddays_spec.csv")
df = read.csv("merged_summary/trans_spec_summary.csv")
df = merge(df, df.bd, by="location")
df$Route = "Total\n(including hospital background)\n"
df.whs = read.csv("merged_summary/trans_spec_src_summary.csv")
df.whs = merge(df.whs, df.bd, by="location")
df.whs$Route = "Sampled source with CDI\n(via ward, spore, hospital-wide)\n"
df.p = rbind(df, df.whs)
df.p = df.p[which(df.p$location!="Community"),]
df.p = cbind(df.p, prop_mean = df.p$mean/df.p$beddays*10000,
               prop_lower = df.p$lower/df.p$beddays*10000, prop_upper = df.p$upper/df.p$beddays*10000)

df.p$loc.lab = factor(df.p$location, levels=c("Acute Medicine & Geratology", "Cardiac", "Gastro", "General Surgery", "Haematology & Oncology",
                                              "ICU", "O&G", "Paediatrics",  "Renal & Transplant", "Specialist Medicine", "Specialist Surgery", "T&O"),
                      labels = c("Acute Medicine\n& Geratology", "Cardiac", "Gastroenterology", "General Surgery", "Haematology\n& Oncology",
                                 "ICU", "Obstetrics\n& Gynaecology", "Paediatrics",  "Renal\n& Transplant", "Specialist Medicine", 
                                 "Specialist Surgery", "Trauma\n& Orthopaedics"))


p = ggplot(df.p, aes(x=loc.lab, y=prop_mean, group=Route, fill=Route)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=prop_lower, ymax=prop_upper, width=0.6), position=position_dodge(width=0.9)) +
  labs(x="Specialty", y="Inpatient aquisitions per 10,000 bed-days\nmean (95% HPD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  geom_text(lineheight = 0.95, aes(x=loc.lab,y=-0.4,label=round(mean,0)),
            position=position_dodge(width=0.9)) 

p = p + annotation_custom(
  grob = textGrob(label = "CDI cases:", hjust = 0, gp = gpar(cex = 0.9)),
  ymin = -0.4,      
  ymax = -0.4,
  xmin = -0.23,         
  xmax = -0.23)


# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
ggsave("merged_summary/plot_by_spec.pdf", gt, width=29.7/2.54, height=21/2.54, useDingbats=FALSE)

