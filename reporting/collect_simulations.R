#!/usr/bin/env Rscript
rm(list = ls())
library(ggplot2)

df = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/5_scenario_parms.csv")
df$true_value = as.numeric(as.character(df$true_value))
df.bg = df[which(df$parameter=="bgroundBeta"),]
#df.bg = df[which(df$parameter=="spore.multiplier"),]
df.bg = cbind(df.bg[order(df.bg$true_value),], sim.num = 1:50)



ggplot(df.bg, aes(x=sim.num, y=mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(x=sim.num, y=true_value), shape=4)

