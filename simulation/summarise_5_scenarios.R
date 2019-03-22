#!/usr/bin/env Rscript

rm(list = ls())
library(ggplot2)

#get parameter estimates and save into single file
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/5_scenarios")
parms = data.frame()
for (folder in list.files(".")) { 
	parms = rbind(parms, read.csv(paste(folder, "inference/parm_summary.csv", sep="/")))
	}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(parms, "5_scenario_parms.csv")

#get source types summary
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/5_scenarios")
srctype = data.frame()
for (folder in list.files(".")) { 
  srctype = rbind(srctype, read.csv(paste(folder, "inference/source_type_summary.csv", sep="/")))
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(srctype, "5_scenario_src_type.csv")

#get transmission times summary
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/5_scenarios")
trans = data.frame()
for (folder in list.files(".")) { 
  trans = rbind(trans, read.csv(paste(folder, "inference/trans_summary.csv", sep="/")))
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(trans, "5_scenario_trans.csv")