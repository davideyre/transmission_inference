rm(list = ls())
library(ggplot2)

load("~/Dropbox/Epi_WGS_MCMC/ward_hosp_comm_sim/simulation_87628393/input/simulation.RData")

idx = 1:length(infectedPatients)


dat = cbind(infectedPatients, infections[infectedPatients,c(1,3)])
dat = cbind(idx, dat[order(infections[infectedPatients,1]),])
colnames(dat) = c("pt_index", "pt", "t_inf", "inf_type")
dat = as.data.frame(dat)


ward.dat = wardLog[which(wardLog[,1] %in% infectedPatients),]
ward.dat = cbind(dat$pt_index[match(ward.dat[,1], dat$pt)], ward.dat)
colnames(ward.dat) = c("pt_index", "pt", "ward", "t_admit", "t_discharge")
ward.dat = as.data.frame(ward.dat, row.names = 1:nrow(ward.dat))

all.dat = cbind(ward.dat, dat$t_inf[match(ward.dat[,2], dat$pt)], dat$inf_type[match(ward.dat[,2], dat$pt)])
colnames(all.dat) = c("pt_index", "pt", "ward", "t_admit", "t_discharge", "t_inf", "inf_type")
all.dat$inf_type[which(all.dat$inf_type==5)] = 2


ggplot(all.dat, aes(x=t_inf, y=pt_index, color=inf_type)) +
   geom_point(shape=1) +
   geom_errorbarh(aes(xmax=t_discharge, xmin=t_admit, y=pt_index))

