rm(list=ls())
options(stringsAsFactors=FALSE)
setwd('~/downloads/data/denogi/')

library(rjags)
library(RColorBrewer)
library(maptools)
library(countrycode)

source('denogi_functions.r')
source('denogi_data.r')

mar <- c(4, 3, 3.5, 1)

mp <- mp[mp$century != 'overall', ]
mp <- mp[order(mp$century, mp$country), ]

a <- array(NA, dim=c(length(unique(mp$century)), 3, length(unique(countries))))
for (i in 1:length(countries)){
    mat <- mp[mp$country==countries[i], grep('mean', names(mp))]
    a[, , i] <- as.matrix(mat)
}
colnames(a) <- grepr('mean', names(mp))
rownames(a) <- unique(mp$century)

datlist <- list(y=a, N=dim(a)[3], K=dim(a)[2], T=dim(a)[1], 
                g0=c(0, 0), G0=diag(1e-7, 2), drift=0.1)
inits <- FA_inits(mp, colnames(a))
initslist <- list(beta=inits$gamma, tau=1/inits$omega^2, x=matrix(inits$xi, ncol=datlist$N))
mdl <- jags.model('denogi_dyn.bug', data=datlist, inits=initslist)
update(mdl, 1e4)
smpls <- coda.samples(mdl, n.iter=5e4, thin=10, 
    variable.names=c('country', 'load', 'intercept'))

countries <- extract_from_coda(smpls, 'country')
par(mfrow=c(2, 3))
for (i in 1:ncol(countries)){
    plot(countries[, i], type='l')
}
loads <- extract_from_coda(smpls, 'load')
par(mfrow=c(2, 3))
for (i in 1:ncol(loads)){
    plot(loads[, i], type='l')
}
intercepts <- extract_from_coda(smpls, 'intercept')
par(mfrow=c(2, 3))
for (i in 1:ncol(intercepts)){
    plot(intercepts[, i], type='l')
}