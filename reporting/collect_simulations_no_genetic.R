#!/usr/bin/env Rscript
rm(list = ls())
library(ggplot2)
library(reshape2)


## accuracy of parameter estimates
df = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenario_parms_no_genetic.csv")
df$true_value = as.numeric(as.character(df$true_value))
simlist = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenarios.csv")
simlist$simulation = paste("simulation_", simlist$simulation, sep="")

setwd("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/")

for (parm in c("bgroundBeta", "wardBeta", "hospBeta", "spore.multiplier", "spore.p", "commBeta")) {
  df.p = df[which(df$parameter==parm),]
  df.p = merge(df.p, simlist, "simulation")
  df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))
  
  p = ggplot(df.p, aes(x=sim.num, y=mean, color=scenario_type)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper)) +
    geom_point(aes(x=sim.num, y=true_value), shape=4) +
    labs(title = paste("Parameter inference performance on simulated data:", parm), x="Simulation identifier", y=parm)
  
  print(p)
  ggsave(paste("50_scenarios_plot_parm_",parm,"_no_genetic.pdf",sep=""), width = 29.7, height = 21.0, units = "cm")
}


## accuracy of source type attribution
df = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenario_src_type_no_genetic.csv")
df = merge(df, simlist, "simulation", all=T)
df = df[which(df$source_type!="Start positive"),]
for (scenario in unique(df$scenario_type)) {
  df.p = df[which(df$scenario_type==scenario),]
  df.p = cbind(df.p[order(df.p$simulation, df.p$source_type),], sim.num = 1:nrow(df.p))
  p = ggplot(df.p, aes(x=sim.num, y=mean, color=source_type)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper))+
    geom_point(aes(x=sim.num, y=simulated_counts), shape=4) +
    labs(title = paste("Transmission routes inference performance on simulated data, scenario:", scenario), x="Simulation and parameter identifier", y="Infections via this route")
  print(p)
  ggsave(paste("50_scenarios_plot_source_",scenario,"_no_genetic.pdf",sep=""), width = 29.7, height = 21.0, units = "cm")
}


## accuracy of overall source type attribution
df = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenario_trans_no_genetic.csv")
df = merge(df, simlist, "simulation", all=T)

df.p = df[which(df$parameter=="Source Types"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))
df.p = melt(df.p[,c("parameter", "sim.num", "scenario_type", "proportion_exact_match", "proportion_hpd_match")], id=c("parameter", "sim.num", "scenario_type"))
p = ggplot(df.p, aes(x=sim.num, y=value, shape=variable, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type inference, performance on simulated data"), x="Simulation identifier", y="Proportion correctly estimated",
       shape="", color="Scenario") +
  ylim(0,1) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  scale_shape_manual(values = c(16, 17), labels = c("Proportion with exact match", "Proportion within 95% HPD")) +
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background")) +
  theme(legend.key.size = unit(1.8, 'lines'))

print (p)
ggsave("50_scenarios_plot_source_types_no_genetic.pdf", width = 29.7, height = 21.0, units = "cm")

df.p = df[which(df$parameter=="Sources"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))
df.p = melt(df.p[,c("parameter", "sim.num", "scenario_type", "proportion_exact_match", "proportion_hpd_match")], id=c("parameter", "sim.num", "scenario_type"))
p = ggplot(df.p, aes(x=sim.num, y=value, shape=variable, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source inference, performance on simulated data"), x="Simulation identifier", y="Proportion correctly estimated",
       shape="", color="Scenario") +
  ylim(0,1) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  scale_shape_manual(values = c(16, 17), labels = c("Proportion with exact match", "Proportion within 95% HPD")) +
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background")) +
  theme(legend.key.size = unit(1.8, 'lines'))

print (p)
ggsave("50_scenarios_plot_sources_no_genetic.pdf", width = 29.7, height = 21.0, units = "cm")

df.p = df[which(df$parameter=="Infection Times"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))
df.p = melt(df.p[,c("parameter", "sim.num", "scenario_type", "proportion_exact_match", "proportion_hpd_match")], id=c("parameter", "sim.num", "scenario_type"))
p = ggplot(df.p, aes(x=sim.num, y=value, shape=variable, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Infection time inference, performance on simulated data"), x="Simulation identifier", y="Proportion correctly estimated",
       shape="", color="Scenario") +
  ylim(0,1) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  scale_shape_manual(values = c(16, 17), labels = c("Proportion with exact match", "Proportion within 95% HPD")) +
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background")) +
  theme(legend.key.size = unit(1.8, 'lines'))

print (p)
ggsave("50_scenarios_plot_inf_times_no_genetic.pdf", width = 29.7, height = 21.0, units = "cm")
