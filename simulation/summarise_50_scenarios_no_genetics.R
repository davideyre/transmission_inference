#!/usr/bin/env Rscript

rm(list = ls())
library(ggplot2)

#get parameter estimates and save into single file
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/50_scenarios_no_genetic")
parms = data.frame()
for (folder in list.files(".")) { 
  filename = paste(folder, "inference/parm_summary.csv", sep="/")
  if(file.exists(filename)) {
    parms = rbind(parms, read.csv(filename))
  }
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(parms, "50_scenario_parms_no_genetic.csv")

#get source types summary
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/50_scenarios_no_genetic")
srctype = data.frame()
for (folder in list.files(".")) { 
  filename = paste(folder, "inference/source_type_summary.csv", sep="/")
  if(file.exists(filename)) {
    srctype = rbind(srctype, read.csv(filename))
  }
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(srctype, "50_scenario_src_type_no_genetic.csv")

#get transmission times summary
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/50_scenarios_no_genetic")
trans = data.frame()
for (folder in list.files(".")) { 
  filename = paste(folder, "inference/trans_summary.csv", sep="/")
  if(file.exists(filename)) {
    trans = rbind(trans, read.csv(filename))
  }
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(trans, "50_scenario_trans_no_genetic.csv")