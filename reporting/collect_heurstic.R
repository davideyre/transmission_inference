#!/usr/bin/env Rscript

rm(list = ls())

#get parameter estimates and save into single file
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data/50_scenarios")
nsim = length(list.files("."))
df = matrix("", nrow=nsim, ncol=14)
colnames(df) = c("simulation", "n", "src_correct", "src_type_correct", 
                 "bg", "bg.true", "ward", "ward.true", "spore", "spore.true",
                 "hosp", "hosp.true", "comm", "comm.true")
i = 1
for (folder in list.files(".")) { 
  filename = paste(folder, "inference/heurstic_trans_stats.csv", sep="/")
  if(file.exists(filename)) {
    data = read.csv(filename)
    n = nrow(data)
    src_correct = sum(data$`Match.Source`) / n
    src_type_correct = sum(data$`Match.Source.Type`) / n
    bg = sum(data$Source.Type=="Background")
    bg.true = sum(data$True.Source.Type=="Background")
    ward = sum(data$Source.Type=="Ward")
    ward.true = sum(data$True.Source.Type=="Ward")
    spore = sum(data$Source.Type=="Spore")
    spore.true = sum(data$True.Source.Type=="Spore")
    hospital = sum(data$Source.Type=="Hospital")
    hospital.true = sum(data$True.Source.Type=="Hospital")
    community = sum(data$Source.Type=="Community")
    community.true = sum(data$True.Source.Type=="Community")
    row = c(folder, n, round(src_correct,4), round(src_type_correct,4),
            bg, bg.true, ward, ward.true, spore, spore.true, hospital, hospital.true, community, community.true)
    df[i,] = row
    i = i+1
  }
}
setwd("/users/bag/deyre/analysis/transmission_inference/sim_data")
write.csv(df, "50_scenario_heuristic.csv", row.names = F)
