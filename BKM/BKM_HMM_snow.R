library(snow)



mod <- jags.model('BKM_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)

coda.samples.wrapper <- function(j)
{ 
  temp.model = jags.model('BKM_Bugs_HMM.R', 
                          inits=list(.RNG.name="base::Wichmann-Hill",
                                     .RNG.seed=j), 
                          data=data, n.chains=cha, n.adapt=ada)
  coda.samples(temp.model, params, n.iter=iter, thin=th) 
}

snow.start.time = proc.time()
cl <- makeCluster(n.chains=cha, "SOCK")
##Make sure the rjags library is loaded in each worker
clusterEvalQ(cl, library(rjags))
##Send data to workers, then fit models. One disadvantage of this
##parallelization is that you lose the ability to watch the progress bar.
clusterExport(cl, list("x","y.obs","N","n.iter","n.thin"))
par.samples = clusterApply(cl, 1:n.chains, coda.samples.wrapper)
##Reorganize 'par.samples' so that it is recognizeable as an 'mcmc.list' object
for(i in 1:length(par.samples)) { par.samples[[i]] <- par.samples[[i]][[1]] }
class(par.samples) <- "mcmc.list"
stopCluster(cl)
snow.end.time = proc.time()
snow.dtime = snow.end.time - snow.start.time
