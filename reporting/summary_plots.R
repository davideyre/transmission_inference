rm(list = ls())
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

sims = read.csv("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/simulation_list.csv")
trans = read.csv("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/combined_simulation_trans_summary.csv")
heur = read.csv("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/heurstic_trans_stats.csv")
colnames(heur) = paste("h_", colnames(heur), sep="")
trans=trans[which(as.character(trans[,1])!="simulation" & as.character(trans[,1])!="simulation_716365"),]
trans[, 3:7] = sapply(trans[, 3:7], as.character)
trans[, 3:7] = sapply(trans[, 3:7], as.numeric)
trans = cbind(trans, sims$type[match(trans$simulation, sims$simulation)])
colnames(trans)[ncol(trans)] = "type"
trans = cbind(trans, heur[match(trans$simulation, heur$h_simulation),])

parm = read.csv("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/combined_simulation_parm_summary.csv")
parm=parm[which(as.character(parm[,1])!="simulation" & as.character(parm[,1])!="simulation_716365"),]
parm[, 3:6] = sapply(parm[, 3:6], as.character)
parm[, 3:6] = sapply(parm[, 3:6], as.numeric)
parm = cbind(parm, sims$type[match(parm$simulation, sims$simulation)])
colnames(parm)[ncol(parm)] = "type"


#covert spore.prob
logistic = function(x) {return(1/(1+exp(-x)))}

parm[which(parm$parameter=="spore_prob"),3:6] = logistic(parm[which(parm$parameter=="spore_prob"),3:6])


srcs = "/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/combined_simulation_source_type_summary.csv"
srcs = read.csv(srcs)
srcs=srcs[which(as.character(srcs[,1])!="simulation" & as.character(srcs[,1])!="simulation_716365"),]
srcs[, 3:7] = sapply(srcs[, 3:7], as.character)
srcs[, 3:7] = sapply(srcs[, 3:7], as.numeric)
srcs = cbind(srcs, sims$type[match(srcs$simulation, sims$simulation)])
colnames(srcs)[ncol(srcs)] = "type"

pdf("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/summary_plots/src_plots.pdf", width=29.7/2.54, height=21/2.54)
p = list()
i=1
for (simType in levels(srcs$type)) {
  srcs.subset = srcs[which(srcs$type==simType),]
  
  typename = factor(simType, 
                      levels=c("ward", "ward_spore_short", "ward_spore_long", 
                               "ward_hospital", "ward_hospital_spore_short"),
                      labels=c("Ward", "Ward + Spore (short)", "Ward + Spore (long)",
                               "Ward + Hospital", "Ward + Spore (short) + Hospital"))
  
  p[[i]] = ggplot(data=srcs.subset, aes(x=factor(gsub("simulation_", "", simulation)), color=factor(source_type))) +
    geom_point(aes(y=simulated_counts), shape=4, position = position_dodge(width=0.4), size=4) +
    geom_point(aes(y=mean), shape=19, position = position_dodge(width=0.4), size=2.5) +
    geom_errorbar(aes(ymin=lower, ymax=upper) , width=0.3, position = position_dodge(width=0.4), size=.6) +
    guides(colour = guide_legend(title="Source type\nMean, 95% HPD\nSimulated value = X", override.aes = list(shape = 19))) + 
    labs(title=paste("Relative source attribution: ", typename, sep=""), x="Simulation", y="Source counts", legend="Source Type") +
    theme(legend.text=element_text(size=14), legend.title=element_text(size=14), axis.text=element_text(size=14), 
          axis.title=element_text(size=14), plot.title=element_text(size=16))
  print(p[[i]])
  i = i+1
}
dev.off()





p = list()
i=1
for (parmType in c("Sources", "Source Types", "Infection Times")) {
   p[[i]] = ggplot(data = trans[which(trans$parameter==parmType),]) +
      geom_dotplot(aes(y=proportion_exact_match, fill="Inference", x=factor(type, 
                                                          levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                   "ward_hospital", "ward_hospital_spore_short"),
                                                          labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                   "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital"))),
                   binaxis = "y", stackdir = "center") +
      scale_y_continuous(limits=c(0,1)) +
      labs(title=paste(parmType, ": match", sep=""), x="Transmission routes included",y="Proportion correctly identified") + 
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="bottom") +
      theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16))  
   
   # if(parmType=="Sources") {
   #   p[[i]] = p[[i]] + 
   #     geom_dotplot(aes(y=h_prop_src_match, fill="Transmission Heuristic", x=factor(type, 
   #                                                               levels=c("ward", "ward_spore_short", "ward_spore_long", 
   #                                                                        "ward_hospital", "ward_hospital_spore_short"),
   #                                                               labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
   #                                                                        "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital"))),
   #                  binaxis = "y", stackdir = "center")
   # } else if (parmType=="Source Types") {
   #   p[[i]] = p[[i]] + 
   #     geom_dotplot(aes(y=h_prop_src_type_match, fill="Transmission Heuristic", x=factor(type, 
   #                                                                     levels=c("ward", "ward_spore_short", "ward_spore_long", 
   #                                                                              "ward_hospital", "ward_hospital_spore_short"),
   #                                                                     labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
   #                                                                              "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital"))),
   #                  binaxis = "y", stackdir = "center", position="dodge")
   # }
   # 
   
   
   p[[i+1]] = ggplot(data = trans[which(trans$parameter==parmType),]) +
      geom_dotplot(aes(y=proportion_hpd_match, fill="Inference", x=factor(type, 
                                                        levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                 "ward_hospital", "ward_hospital_spore_short"),
                                                        labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                 "Ward + \nHospital", "Ward + \nSpore (short) +\nHospital"))),
                   binaxis = "y", stackdir = "center") +
      scale_y_continuous(limits=c(0,1)) +
      labs(title=paste(parmType, ": within 95% HPD", sep=""), x="Transmission Routes Included",y="Proportion within 95% HPD") + 
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position="bottom") +
      theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) 
   i = i+2
}


pdf("/Users/davideyre/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/summary_plots/plots.pdf", width=29.7/2.54, height=21/2.54)
grid.arrange( p[[1]], p[[2]], ncol=2)
grid.arrange( p[[3]], p[[4]], ncol=2)
grid.arrange( p[[5]], p[[6]], ncol=2)


colourCount = 25
getPalette = colorRampPalette(brewer.pal(5, "Blues"))


parmType = "beta0"
p = ggplot(data = parm[which(parm$parameter==parmType),], aes(fill=simulation, y=mean, x=factor(type, 
                                                                                            levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                                                     "ward_hospital", "ward_hospital_spore_short"),
                                                                                            labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                                                     "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital")))) +
   geom_bar(position="dodge", stat="identity") +
   geom_errorbar(aes(ymax = upper, ymin=lower), position="dodge") +
   labs(title="Beta, background", x="Transmission routes included",y="Mean (95% HPD)") + 
   scale_fill_manual(values=getPalette(colourCount)) +
   theme_bw() +
   theme(legend.position="none") +
   theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) 
print(p)

parmType = "beta1"
p=ggplot(data = parm[which(parm$parameter==parmType),], aes(fill=simulation, y=mean, x=factor(type, 
                                                                       levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                                "ward_hospital", "ward_hospital_spore_short"),
                                                                       labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                                "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital")))) +
   geom_bar(position="dodge", stat="identity") +
   geom_errorbar(aes(ymax = upper, ymin=lower), position="dodge") +
   labs(title="Beta, ward", x="Transmission routes included",y="Mean (95% HPD)") + 
   scale_fill_manual(values=getPalette(colourCount)) +
   theme_bw() +
   theme(legend.position="none") +
   theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) 
print(p)

parmType = "beta2"
p=ggplot(data = parm[which(parm$parameter==parmType),], aes(fill=simulation, y=mean, x=factor(type, 
                                                                                            levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                                                     "ward_hospital", "ward_hospital_spore_short"),
                                                                                            labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                                                     "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital")))) +
   geom_bar(position="dodge", stat="identity") +
   geom_errorbar(aes(ymax = upper, ymin=lower), position="dodge") +
   labs(title="Beta, hospital", x="Transmission routes included",y="Mean (95% HPD)") + 
   scale_fill_manual(values=getPalette(colourCount)) +
   theme_bw() +
   theme(legend.position="none") +
   theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) 
print(p)


parmType = "spore_prob"
p=ggplot(data = parm[which(parm$parameter==parmType),], aes(fill=simulation, y=mean, x=factor(type, 
                                                                                            levels=c("ward", "ward_spore_short", "ward_spore_long", 
                                                                                                     "ward_hospital", "ward_hospital_spore_short"),
                                                                                            labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
                                                                                                     "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital")))) +
   geom_bar(position="dodge", stat="identity") +
   geom_errorbar(aes(ymax = upper, ymin=lower), position="dodge") +
   labs(title="Spore probability", x="Transmission routes included",y="Mean (95% HPD)") + 
   scale_fill_manual(values=getPalette(colourCount)) +
   theme_bw() +
   theme(legend.position="none") +
   theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) +
   scale_y_continuous(breaks = 1:10/10)
print(p)

dev.off()


#check for all parm - look fine

# for (parmType in unique(parm$parameter)) {
#    p=ggplot(data = parm[which(parm$parameter==parmType),], aes(fill=simulation, y=mean, x=factor(type, 
#                                                                                                levels=c("ward", "ward_spore_short", "ward_spore_long", 
#                                                                                                         "ward_hospital", "ward_hospital_spore_short"),
#                                                                                                labels=c("Ward", "Ward + \nSpore (short)", "Ward + \nSpore (long)",
#                                                                                                         "Ward + \nHospital", "Ward + \nSpore (short) + \nHospital")))) +
#       geom_bar(position="dodge", stat="identity") +
#       geom_errorbar(aes(ymax = upper, ymin=lower), position="dodge") +
#       labs(title=parmType, x="Transmission routes included",y="Mean (95% HPD)") + 
#       scale_fill_manual(values=getPalette(colourCount)) +
#       theme_bw() +
#       theme(legend.position="none") +
#       theme(legend.text=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=16)) 
#    print(p)
# }
