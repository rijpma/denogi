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
mp <- mp[order(mp$region, mp$century, mp$country), ]
mp$ncentury <- as.numeric(as.factor(mp$century))
levels <- as.factor(paste(mp$region, mp$century))

m <- factanal(mp[complete.cases(mp), grepl('mean', names(mp))], 
    factors=1, scores='Bartlett')

means <- mp[, grepl('mean', names(mp))]
rownames(means) <- paste0(mp$country, mp$century)
sds <- mp[, grepl('sd', names(mp))]
rownames(sds) <- paste0(mp$country, mp$century)

inits <- FA_inits(means, c('mean_celib', 'mean_smam', 'mean_cmplx'))
il <- list(beta=inits$gamma, tau=1/(inits$omega^2), x=inits$xi,
    .RNG.name="base::Mersenne-Twister", .RNG.seed=217563)

datlist <- list(y=means, N=nrow(means), K=3, # sigma=sds,
    g0=c(0, 0), G0=diag(1e-7, 2), 
    level=levels, M=length(unique(levels)))

fit <- jags.model('denogi_ml.bug', data=datlist, inits=il)
update(fit, 5e3)
smpls <- coda.samples(fit, n.iter=25e4, thin=100,
    variable.names=c('load', 'intercept', 'country', 'region', 'y'))

rd <- raftery.diag(smpls)
# hist(rd[[1]]$resmatrix[,2])
# 50k observations good enough for most par., some ctrs need 250k

countries <- extract_from_coda(smpls, 'country')
regions <- extract_from_coda(smpls, 'region')
intercepts <- extract_from_coda(smpls, 'intercept')
loads <- extract_from_coda(smpls, 'load')
ys <- extract_from_coda(smpls, '^y\\[')

regionshat <- sumstatsDF(regions)
regionshat$region <- vstrsplit(unique(as.character(levels)), ' c')[,1]
regionshat$century <- as.factor(vstrsplit(unique(as.character(levels)), ' c')[,2])

loadshat <- sumstatsDF(loads)
rownames(loadshat) <- c('Celibacy', 'SMAM', 'Complex HH')
interceptshat <- sumstatsDF(intercepts)
rownames(interceptshat) <- c('Celibacy', 'SMAM', 'Complex HH')

ctrshat <- data.frame(sumstatsDF(countries),
    region=mp$region, iso3=mp$iso3, century=mp$century,
    ncentury=as.numeric(as.factor(mp$century)), country=mp$country)

yhat <- sumstatsDF(ys)
yhat$country <- mp$country
yhat$century <- mp$century
yhat$ncentury <- as.numeric(as.factor(mp$century))
yhat$iso3 <- mp$iso3
yhat$vrb <- rep(colnames(datlist$y), each=200)

write.csv(format(interceptshat, digits=3), 'intercepts.csv')
write.csv(format(loadshat, digits=3), 'loadings.csv')
write.csv(ctrshat, 'countries.csv', row.names=F)

# matrixes with odds lv-score i > j
regionprobs <- matrix(NA, ncol=ncol(regions), nrow=ncol(regions))
for (i in 1:ncol(regions)){
    probs <- colSums(regions[, -i] > regions[, i]) / nrow(regions)
    regionprobs[-i, i] <- probs
}
colnames(regionprobs) <- rownames(regionprobs) <- gsub('[a-z ]', '', unique(levels))

ctrprobs <- matrix(NA, ncol=ncol(countries), nrow=ncol(countries))
for (i in 1:ncol(countries)){
    probs <- colSums(countries[, -i] > countries[, i]) / nrow(countries)
    ctrprobs[-i, i] <- probs
}
colnames(ctrprobs) <- rownames(ctrprobs) <- paste0(mp$country, mp$century)

write.csv(regionprobs, 'regionprobs.csv')
write.csv(ctrprobs, 'ctrprobs.csv')

pdf('mpbyregion.pdf')
par(mfrow=c(2,2 ), mar=mar, font.main=1)
for (reg in unique(regionshat$region)){
    d <- regionshat[regionshat$region==reg, ]
    plot(1:5, d$q50, lwd=2,
        type='l', ylim=c(-1.5, 1.5), bty='n',
        xlab='century', ylab='Marriage pattern', axes=F)
    lines(1:5, d$q05, lty=3)
    lines(1:5, d$q95, lty=3)
    axis(1, labels=unique(regionshat$century), at=1:5)
    axis(2)
    title(main=reg, line=-0.6)
}
dev.off()

pdf('countries.pdf')
par(mfrow=c(2, 2), mar=mar)
rg <- range(ctrshat[, c('q05', 'q50', 'q95')])
for (ctr in unique(ctrshat$country)){
    ds <- ctrshat[ctrshat$country==ctr, ]
    plot(q50 ~ ncentury, data=ds, type='l', lwd=2, ylim=rg, axes=F, ylab='Marriage pattern')
    lines(q95 ~ ncentury, data=ds, type='l', lty=3)
    lines(q05 ~ ncentury, data=ds, type='l', lty=3)
    title(main=ctr, line=-0.7)
    axis(1, labels=unique(centuries)[ds$ncentury], at=ds$ncentury)
    axis(2)
}
dev.off()

countrieshat <- ctrshat[order(ctrshat$q50), ]
par(mfrow=c(1, 1))
plot(countrieshat$q50, 1:datlist$N, , type='p')
segments(countrieshat$q05, 1:datlist$N, countrieshat$q95, 1:datlist$N)
place_text(countrieshat, countrieshat$iso3)

data(wrld_simpl)
c18 <- countrieshat[countrieshat$century=='c18', ]
c18 <- aggregate(q50 ~ iso3, data=c18, mean)

pdf('eumarpattern_c18.pdf')
cols = brewer.pal(5, "RdPu")
wrld_simpl@data <- data.frame(wrld_simpl@data, 
    c18[match(wrld_simpl@data$ISO3, c18$iso3), ])
eur_simpl <- wrld_simpl[!is.na(wrld_simpl@data$iso3), ]
eur_simpl@data$col <- as.character(factor(cut(eur_simpl@data$q50, breaks=5), labels=cols))
par(mfrow=c(1,1))
plot(eur_simpl, col=eur_simpl$col, xlim=c(-10, 40), ylim=c(40, 70))
dev.off()

pdf('impcheck.pdf')
par(mfrow=c(2, 2))
for (iso3 in unique(yhat$iso3)){
    dat <- yhat[yhat$iso3==iso3, ]
    plot(y=range(dat[, 1:4]), x=range(dat$ncentury), type='n', 
        main=paste('celib', unique(iso3)), bty='l')
    lines(mean ~ ncentury, data=dat[dat$vrb=='mean_celib', ])
    lines(q95 ~ ncentury, data=dat[dat$vrb=='mean_celib', ], lty=2)
    lines(q05 ~ ncentury, data=dat[dat$vrb=='mean_celib', ], lty=2)
    points(mean_celib ~ ncentury, data=mp[mp$iso3==iso3, ])

    plot(y=range(dat[, 1:4]), x=range(dat$ncentury), type='n', 
        main='smam', bty='l')
    lines(mean ~ ncentury, data=dat[dat$vrb=='mean_smam', ])
    lines(q95 ~ ncentury, data=dat[dat$vrb=='mean_smam', ], lty=2)
    lines(q05 ~ ncentury, data=dat[dat$vrb=='mean_smam', ], lty=2)
    points(mean_smam ~ ncentury, data=mp[mp$iso3==iso3, ])

    plot(y=range(dat[, 1:4]), x=range(dat$ncentury), type='n', 
        main='cmplx', bty='l')
    lines(mean ~ ncentury, data=dat[dat$vrb=='mean_cmplx', ])
    lines(q95 ~ ncentury, data=dat[dat$vrb=='mean_cmplx', ], lty=2)
    lines(q05 ~ ncentury, data=dat[dat$vrb=='mean_cmplx', ], lty=2)
    points(mean_cmplx ~ ncentury, data=mp[mp$iso3==iso3, ])

    plot(y=range(ctrshat[ctrshat$iso3==iso3, 1:4]),
         x=range(ctrshat$ncentury[ctrshat$iso3==iso3]), 
         type='n', main='', bty='l')
    lines(q50 ~ ncentury, data=ctrshat[ctrshat$iso3==iso3, ], type='l', main='')
    lines(q05 ~ ncentury, data=ctrshat[ctrshat$iso3==iso3, ], lty=2)
    lines(q95 ~ ncentury, data=ctrshat[ctrshat$iso3==iso3, ], lty=2)
}
dev.off()