population = 50
true.p.geom.inf = 0.2
true.p.geom.sample = 0.5

#true dataset of infection times, geom with parameter 0.2
true.inf.times = rgeom(population, true.p.geom.inf)
sample.times = true.inf.times + rgeom(population, true.p.geom.sample) #delay observation 

#likelihood parm | infection times
ll.parm = function(inftimes, p.geom.inf, p.import) {
  ll.data = ifelse(inftimes<0, log(p.import), log(1-p.import) + dgeom(inftimes, p.geom.inf, log=T))
  return(sum(ll.data))
}

# likelihood of sampled data | infection times
ll.sampling = function(inftimes, samplingtimes, p.geom.sample) {
  ll.data = dgeom(samplingtimes-inftimes, p.geom.sample, log=T)
  return(sum(ll.data))
}

ll.total = function(inftimes, samplingtimes, p.geom.sample, p.geom.inf, p.import) {
  ll.parm(inftimes, p.geom.inf, p.import) + ll.sampling(inftimes, samplingtimes, p.geom.sample)
}

#true ll
ll.total(true.inf.times, sample.times, 0.5, 0.2, 0)

#alternative, with all inf.times set 5 days earlier
for (offset in 1:20) {
  print (ll.total(true.inf.times-offset, sample.times, p.geom.sample=0.06, p.geom.inf=0.1, p.import=0.94))
}


#the extreme
for (offset in 1:50) {
  print (ll.total(true.inf.times-offset, sample.times, p.geom.sample=0.05, p.geom.inf=0.1, p.import=0.999))
}

