#!/usr/bin/env Rscript

rm(list = ls())
library(ggplot2)

setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/5_scenarios")
parms = data.frame()
for (folder in list.files(".")) { 
	parms = rbind(parms, read.csv(paste(folder, "inference/parm_summary.csv", sep="/")))
	}

setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")

write.csv(parms, "5_scenario_parms.csv")
