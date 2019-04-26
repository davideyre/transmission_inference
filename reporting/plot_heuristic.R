#!/usr/bin/env Rscript
rm(list = ls())
library(reshape2)


df = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenario_heuristic.csv",
              stringsAsFactors = T)

simlist = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/sim_data/50_scenarios.csv",
                   stringsAsFactors = T)
simlist$simulation = paste("simulation_", simlist$simulation, sep="")
simlist = simlist[,c("simulation", "scenario_type")]


df.all = merge(df, simlist, "simulation", all=T)


### PLOT BIAS BY TYPE

df = cbind(df.all, (df$bg-df$bg.true)/df$n, (df$ward-df$ward.true)/df$n, (df$spore-df$spore.true)/df$n, 
           (df$hosp-df$hosp.true)/df$n, (df$comm-df$comm.true)/df$n)

colnames(df) = c(colnames(df)[1:15], "bg.diff", "ward.diff", "spore.diff", "hosp.diff", "comm.diff")

df = melt(df[,c("simulation", "scenario_type", "bg.diff", "ward.diff", "spore.diff", "hosp.diff", "comm.diff")], 
          id=c("simulation", "scenario_type"))


df.p = df[which(df$variable=="ward.diff"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Sources attributed to ward (estimated-true) / n",
       shape="", color="Scenario") 
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic_bias_ward.pdf", width = 29.7, height = 21.0, units = "cm")

df.p = df[which(df$variable=="bg.diff"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Sources attributed to background (estimated-true) / n",
       shape="", color="Scenario") 
scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic_bias_bground.pdf", width = 29.7, height = 21.0, units = "cm")

df.p = df[which(df$variable=="spore.diff"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Sources attributed to spore (estimated-true) / n",
       shape="", color="Scenario") 
scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic_bias_spore.pdf", width = 29.7, height = 21.0, units = "cm")


df.p = df[which(df$variable=="hosp.diff"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Sources attributed to hospital (estimated-true) / n",
       shape="", color="Scenario") 
scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic_bias_hospital.pdf", width = 29.7, height = 21.0, units = "cm")


df.p = df[which(df$variable=="comm.diff"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Sources attributed to community (estimated-true) / n",
       shape="", color="Scenario") 
scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic_bias_community.pdf", width = 29.7, height = 21.0, units = "cm")


### PLOT SUMMARY

df = melt(df.all[,c("simulation", "src_correct", "src_type_correct", "scenario_type")], id=c("simulation", "scenario_type"))

df.p = df[which(df$variable=="src_type_correct"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source type heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Proportion correctly estimated",
       shape="", color="Scenario") +
  ylim(0,1) +
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_type_heuristic.pdf", width = 29.7, height = 21.0, units = "cm")




df.p = df[which(df$variable=="src_correct"),]
df.p = cbind(df.p[order(df.p$scenario_type),], sim.num = 1:nrow(df.p))

p = ggplot(df.p, aes(x=sim.num, y=value, color=scenario_type)) +
  geom_point() +
  labs(title = paste("Transmission source heurstic inference, performance on simulated data"), 
       x="Simulation identifier", y="Proportion correctly estimated",
       shape="", color="Scenario") +
  ylim(0,1) +
  scale_color_discrete(labels = c("Hospital background & ward", "Hospital background, ward\n& hospital-wide", 
                                  "Hospital background, ward,\nhospital-wide & long spore", "Hospital background, ward,\nhospital-wide & short spore",
                                  "Hospital background, ward, hospital-wide,\nspore & community background"))

print (p)
ggsave("50_scenarios_plot_source_heuristic.pdf", width = 29.7, height = 21.0, units = "cm")