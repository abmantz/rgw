rGWsqrtdist = function(n, a=2.0) (runif(n) * 2.0*(sqrt(a)-sqrt(1.0/a))/2.0 + sqrt(1.0/a))^2

GoodmanWeare = function(ensemble, lnpost, Nsteps, current.lnP=NULL, mc.cores=getOption("mc.cores", 1L), ...) {
  ## ensemble should be an Nparam*Nwalkers matrix
  nwalkers = ncol(ensemble)
  if (nwalkers == 0 || nwalkers %% 2 == 1) stop('Number of walkers must be positive and even')
  nparam = nrow(ensemble)
  group1 = 1:(nwalkers/2)
  group2 = nwalkers/2 + group1
  if (length(current.lnP) != nwalkers)
    current.lnP = simplify2array(mclapply(1:nwalkers, function(k) lnpost(ensemble[,k], ...), mc.cores=mc.cores, mc.allow.recursive=FALSE))
  for (i in 1:Nsteps) {
    stretches = rGWsqrtdist(nwalkers)
    lnr = log(runif(nwalkers))
    js = c(sample(group2, length(group1), replace=TRUE), sample(group1, length(group2), replace=TRUE))
    newens = simplify2array(mclapply(group1, function(k){
      prop = ensemble[,js[k]] + stretches[k] * (ensemble[,k] - ensemble[,js[k]])
      trial.logP = lnpost(prop, ...)
      lnq = (nparam - 1.0) * log(stretches[k]) + trial.logP - current.lnP[k]
      if (lnr[k] <= lnq) {
        ensemble[,k] = prop
        current.lnP[k] = trial.logP
      }
      c(ensemble[,k], current.lnP[k])
    }, mc.cores=mc.cores, mc.allow.recursive=FALSE))
    ensemble[,group1] = newens[1:nparam,]
    current.lnP[group1] = newens[nparam+1,]
    newens = simplify2array(mclapply(group2, function(k){
      prop = ensemble[,js[k]] + stretches[k] * (ensemble[,k] - ensemble[,js[k]])
      trial.logP = lnpost(prop, ...)
      lnq = (nparam - 1.0) * log(stretches[k]) + trial.logP - current.lnP[k]
      if (lnr[k] <= lnq) {
        ensemble[,k] = prop
        current.lnP[k] = trial.logP
      }
      c(ensemble[,k], current.lnP[k])
    }, mc.cores=mc.cores, mc.allow.recursive=FALSE))
    ensemble[,group2] = newens[1:nparam,]
    current.lnP[group2] = newens[nparam+1,]
  }
  list(ensemble=ensemble, current.lnP=current.lnP)
}

GoodmanWeare.rem = function(post, lnpost, thin=1, mention.every=NA, save.every=NA, save.file=NA, ...) {
  ## post should be an Nparam*Nwalkers*Nsteps array with the initial walker positions
  ## set in post[,,1]. The ensemble will be saved in post[,,i] every thin
  ## iterations until post is filled.
  res = list()
  for (i in 2:dim(post)[3]) {
    res = GoodmanWeare(post[,,i-1], lnpost, thin, res$current.lnP, ...)
    post[,,i] = res$ensemble
    if (!is.na(mention.every) & i %% mention.every == 0) message(paste('Finished iteration', i))
    if (!is.na(save.every) & !is.na(save.file) & i %% save.every == 0) save(post, file=save.file)
  }
  post
}
