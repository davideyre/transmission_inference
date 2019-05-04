#!/usr/bin/env Rscript
rm(list = ls())

library(ggplot2)
library(zoo)

stList = c(1, 2, 3, 5, 6, 8, 10, 11, 42, 44, "_other")
 
# ##collect parameter estimates
# for (st in stList) {
#   path = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/parm_summary.csv", sep="")
#   df = read.csv(path)
#   df = cbind(rep(paste("st", st, sep=""), nrow(df)), df)
#   colnames(df)[1] = "st"
#   if(st=="1") {
#     out = df
#   } else {
#     out = rbind(out, df)
#   }
# }
# outFile = "/home/davideyre/transmission_inference/nejm/parm_compare_by_st.csv"
# write.csv(out, outFile, row.names = F)
# 
# 
# ##collect transmission summaries
# for (st in stList) {
#   path = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/source_type_summary.csv", sep="")
#   df = read.csv(path)
#   df = cbind(rep(paste("st", st, sep=""), nrow(df)), df, df$mean/sum(df$map), df$lower/sum(df$map), df$upper/sum(df$map))
#   colnames(df)[1] = "st"
#   colnames(df)[7:9] = c("prop_mean", "prop_lower", "prop_upper")
#   if(st=="1") {
#     out = df
#   } else {
#     out = rbind(out, df)
#   }
# }
# outFile = "/home/davideyre/transmission_inference/nejm/trans_src_type_compare_by_st.csv"
# write.csv(out, outFile, row.names = F)

##collect transmission summaries
for (st in stList) {
  path.loc = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/chain_inf_locations.txt", sep="")
  df.loc = read.table(path.loc, header=T, sep="\t", stringsAsFactors = F)
  path.src = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/chain_inf_source_types.txt", sep="")
  df.src = read.table(path.src, header=T, sep="\t", stringsAsFactors = F)
  if(st=="1") {
    out.loc = df.loc[,1:ncol(df.loc)-1]
    out.src = df.src[,1:ncol(df.src)-1]
  } else {
    out.loc =cbind(out.loc, df.loc[,1:ncol(df.loc)-1])
    out.src =cbind(out.src, df.src[,1:ncol(df.src)-1])
  }
}
#discard burn-in
out.loc = out.loc[((nrow(out.loc)*.2)+1):nrow(out.loc),]
out.src = out.src[((nrow(out.src)*.2)+1):nrow(out.src),]
#thin by factor of 10
out.loc = out.loc[seq(1, nrow(out.loc), by = 10),]
out.src = out.src[seq(1, nrow(out.src), by = 10),]


out.loc.mat = as.matrix(out.loc)
out.src.mat = as.matrix(out.src)

out.comb.mat = matrix(paste(out.src.mat, out.loc.mat, sep="_"),
                      nrow=nrow(out.loc.mat), ncol=ncol(out.loc.mat))

wardList = sort(unique(as.vector(out.comb.mat)))

s.df = matrix(NA, ncol = length(wardList), nrow=nrow(out.loc))
colnames(s.df) = wardList

for (i in 1:nrow(out.comb.mat)) {
  iter = as.vector(out.comb.mat[i,])
  for (ward in wardList) {
    s.df[i,ward] = sum(iter==ward)
  }
}

#save file with one row per iteration and a count of each route/ward combination for each column
outFile = "/home/davideyre/transmission_inference/nejm/trans_location.csv"
write.csv(s.df, outFile, row.names = F)

#read in lookup table that can join specialty to ward
spec.df = read.csv("/home/davideyre/transmission_inference/nejm/spec_lookup.csv", stringsAsFactors = F)
wards = gsub("^[0-9]_", "", wardList)
spec.search = match(wards, spec.df$location)
spec = spec.df[spec.search, "spec"]
hosp.fix = spec.df[spec.search, "hospital"]

lookup = as.data.frame(cbind(wardList, gsub("^[0-9]_", "", wardList), 
      hosp.fix, 
      substr(wardList, 1, 1),
      spec,
      paste(hosp.fix, spec, sep="_")), stringsAsFactors=F)
colnames(lookup) = c("tag", "ward", "hospital", "src", "spec", "hosp_spec")

#create by specialty summary
u.spec = unique(spec)
spec.list = matrix(NA, ncol = length(u.spec), nrow=nrow(out.loc))
spec.list.trans = matrix(NA, ncol = length(u.spec), nrow=nrow(out.loc))
colnames(spec.list) = u.spec
colnames(spec.list.trans) = u.spec

for (i in 1:nrow(spec.list)) { #each iteration
  for (spec in u.spec) { #each specialty
    spec.list[i,spec] = sum(s.df[i,which(lookup$spec==spec)])
    spec.list.trans[i,spec] = sum(s.df[i,which(lookup$spec==spec & lookup$src!=0 &lookup$src!=3)]) #exclude hospital and community background
  }
}

spec.summary = cbind(location=u.spec, mean=apply(spec.list, 2, mean), 
                     lower=apply(spec.list, 2, quantile, probs=0.025), 
                     upper=apply(spec.list, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_spec_summary.csv"
write.csv(spec.summary, outFile, row.names = F)

spec.summary.trans = cbind(location=u.spec, mean=apply(spec.list.trans, 2, mean), 
                     lower=apply(spec.list.trans, 2, quantile, probs=0.025), 
                     upper=apply(spec.list.trans, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_spec_src_summary.csv"
write.csv(spec.summary.trans, outFile, row.names = F)




#create by hospital summary
u.hospital = unique(lookup$hospital)
hospital.list = matrix(NA, ncol = length(u.hospital), nrow=nrow(out.loc))
colnames(hospital.list) = u.hospital

for (i in 1:nrow(hospital.list)) { #each iteration
  for (hospital in u.hospital) { #each specialty
    hospital.list[i,hospital] = sum(s.df[i,which(lookup$hospital==hospital)])
  }
}

hospital.summary = cbind(location=u.hospital, mean=apply(hospital.list, 2, mean), 
                     lower=apply(hospital.list, 2, quantile, probs=0.025), 
                     upper=apply(hospital.list, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_hospital_summary.csv"
write.csv(hospital.summary, outFile, row.names = F)


#create hospital speciality summary
u.hosp_spec = unique(lookup$hosp_spec)
hosp_spec.list = matrix(NA, ncol = length(u.hosp_spec), nrow=nrow(out.loc))
colnames(hosp_spec.list) = u.hosp_spec

for (i in 1:nrow(hosp_spec.list)) { #each iteration
  for (hosp_spec in u.hosp_spec) { #each specialty
    hospital = unlist(strsplit(hosp_spec, "_"))[1]
    spec = unlist(strsplit(hosp_spec, "_"))[2]
    hosp_spec.list[i,hosp_spec] = sum(s.df[i,which(lookup$hospital==hospital & lookup$spec==spec)])
  }
}

hosp_spec.summary = cbind(location=u.hosp_spec, mean=apply(hosp_spec.list, 2, mean), 
                         lower=apply(hosp_spec.list, 2, quantile, probs=0.025), 
                         upper=apply(hosp_spec.list, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_hospital_spec_summary.csv"
write.csv(hosp_spec.summary, outFile, row.names = F)




## BY ST comparisons - parameters


out = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/nejm/parm_compare_by_st.csv")
parmList = unique(out$parameter)
for (parm in parmList) {
  df.p = out[which(out$parameter==parm),]
  p = ggplot(df.p, aes(x=st, y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper)) +
    labs(title = paste("Parameter inference, by ST:", parm), 
         x="Sequence Type", y="Mean (95% HPD)")
  print(p)
}


## src type by ST

out = read.csv("/Users/davideyre/Drive/academic/research/transmission_modelling/cdiff_transmission_inference/xcode_project/nejm/trans_src_type_compare_by_st.csv")
srcList = unique(out$source_type)
for (src in srcList) {
  df.p = out[which(out$source_type==src),]
  p = ggplot(df.p, aes(x=st, y=prop_mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=prop_lower, ymax=prop_upper)) +
    labs(title = paste("Parameter inference, by ST:", src), 
         x="Sequence Type", y="Proportion of cases aquired via this route, mean (95% HPD)")
  print(p)
}


## TRANSMISSION TIME ANALYSIS

# get all transmission times where t=1 (14-Jan-2007)
for (st in stList) {
  path = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/chain_inf_times.txt", sep="")
  df = read.table(path, header=T, sep="\t", stringsAsFactors = F)
  path.src = paste("/home/davideyre/transmission_inference/nejm/st", st, "/inference/chain_inf_source_types.txt", sep="")
  df.src = read.table(path.src, header=T, sep="\t", stringsAsFactors = F)
  if(st=="1") {
    out = df[,1:ncol(df)-1]
    out.src = df.src[,1:ncol(df.src)-1]
  } else {
    out = cbind(out, df[,1:ncol(df)-1])
    out.src =cbind(out.src, df.src[,1:ncol(df.src)-1])
  }
}
#discard burn-in
out = out[((nrow(out)*.2)+1):nrow(out),]
out.src = out.src[((nrow(out.src)*.2)+1):nrow(out.src),]
#thin by factor of 10
out = out[seq(1, nrow(out), by = 10),]
out.src = out.src[seq(1, nrow(out.src), by = 10),]

out.mat = as.matrix(out)
out.src.mat = as.matrix(out.src)

out.comb = ifelse(out.src.mat==1 | out.src.mat==2 | out.src.mat==5, out.mat, NA)




#now have iterations as rows and cases as columns
offset = as.Date("2007-01-14", "%Y-%m-%d")
out.dt = as.matrix(offset+out.mat, nrow = nrow(out.mat), ncol = ncol(out.mat))
out.dt.trans = as.matrix(offset+out.comb, nrow = nrow(out.comb), ncol = ncol(out.comb))
colnames(out.dt) = colnames(out.mat)
rownames(out.dt) = rownames(out.mat)
colnames(out.dt.trans) = colnames(out.mat)
rownames(out.dt.trans) = rownames(out.mat)
out.q = matrix(as.character(as.yearqtr(out.dt, format = "%Y-%m-%d")), nrow = nrow(out.mat), ncol = ncol(out.mat))
out.q.trans = matrix(as.character(as.yearqtr(out.dt.trans, format = "%Y-%m-%d")), nrow = nrow(out.mat), ncol = ncol(out.mat))
quarters = sort(unique(out.q[1,]))

q.list = matrix(NA, nrow=nrow(out.q), ncol = length(quarters))
q.list.trans = matrix(NA, nrow=nrow(out.q), ncol = length(quarters))
colnames(q.list) = quarters
colnames(q.list.trans) = quarters
for (i in 1:nrow(out.q)) {
  for(q in quarters) {
    q.list[i,q] = sum(out.q[i,]==q)
    q.list.trans[i,q] = sum(out.q.trans[i,]==q, na.rm=T)
  }
}


q.summary = cbind(quarter=quarters, mean=apply(q.list, 2, mean), 
                          lower=apply(q.list, 2, quantile, probs=0.025), 
                          upper=apply(q.list, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_quarter_summary.csv"
write.csv(q.summary, outFile, row.names = F)

q.summary.trans = cbind(quarter=quarters, mean=apply(q.list.trans, 2, mean), 
                  lower=apply(q.list.trans, 2, quantile, probs=0.025), 
                  upper=apply(q.list.trans, 2, quantile, probs=0.975))
outFile = "/home/davideyre/transmission_inference/nejm/trans_quarter_whs_summary.csv"
write.csv(q.summary.trans, outFile, row.names = F)

#add ward_spore_hosp to quarter summary