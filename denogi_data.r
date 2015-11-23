options(stringsAsFactors=FALSE)

celib <- read.table('denogi_celib.csv', skip=1, na.string='n/a', sep=',', header=T)
smam <- read.table('denogi_smam.csv', skip=1, na.string='n/a', sep=',', header=T)
cmplx <- read.table('denogi_complex.csv', skip=1, na.string='n/a', sep=',', header=T)

celib$X[!celib$X %in% smam$X]
smam$X[!smam$X %in% celib$X]
smam$X[smam$X=='France (whole)'] <- 'France (Whole)'

celib$X[!celib$X %in% cmplx$X]
cmplx$X[!cmplx$X %in% celib$X]

smam$X[!smam$X %in% cmplx$X]
cmplx$X[!cmplx$X %in% smam$X]

centuries = rep(c('c10-5', 'c16', 'c17', 'c18', 'c19', 'overall'), each=3)
variables = c('country', rep(c('mean', 'sd', 'N'), length(centuries) /3))

celib <- tablefix(celib)
smam  <- tablefix(smam)
cmplx <- tablefix(cmplx)

mp <- merge(celib, smam, by=c('country', 'century'), all=T)
mp <- merge(mp, cmplx, by=c('country', 'century'), all=T)
names(mp)[-c(1:2)] <- paste0(rep(c('mean', 'sd', 'N'), 3), rep(c('_celib', '_smam', '_cmplx'), each=3))

mp$nobs <- apply(mp[-c(1, 2)], 1, function(x) sum(!is.na(x)))

mp <- mp[!mp$country=='', ]
countries <- unique(mp$country)

regions <- countrycode(countries, 'country.name', 'region')
mp$region <- regions[match(mp$country, countries)]
mp$region[mp$country=='England'] <- "Western Europe"
mp$region[mp$country=='Scotland'] <- "Western Europe"

iso3 <- countrycode(countries, 'country.name', 'iso3c')
mp$iso3 <- iso3[match(mp$country, countries)]
mp$iso3[mp$country=='England'] <- "GBR"
mp$iso3[mp$country=='Scotland'] <- "GBR"

### GDP data from clio-infra

# gdp <- read.csv('~/dropbox/cliodata/GDPperCapita.csv', skip=2)
# names(gdp)[1:2] <- c('ccode', 'country')
# gdp <- gdp[!is.na(gdp$ccode), ]
# yrs <- names(gdp)[-c(1:2)]
# gdp$iso3 <- countrycode(gdp$ccode, 'iso3n', 'iso3c')
# gdp[is.na(gdp$iso3), c('iso3', 'ccode', 'country')]
# cat(gdp[is.na(gdp$iso3), c('country')], sep='\n')
# gdp$iso3[gdp$country=="Federal Republic of Germany (until 1990)"] <- "DEU"
# gdp$iso3[gdp$country=="Kosovo"] <- "KOS"
# gdp$iso3[gdp$country=="Serbia and Montenegro (until 2006)"] <- "SCG"
# gdp$iso3[gdp$country=="Socialist Federal Republic of Yugoslavia (until 1992)"] <- "YUG" 
# gdp$iso3[gdp$country=="USSR (until 1991)"] <- "SUN"
# gdp$iso3[gdp$country=="Pacific Islands (Trust Territory) (1947-1986)"] <- "PCI"
# gdp$iso3[gdp$country=="Tibet"] <- "TIB"
# gdp$iso3[gdp$country=="Taiwan"] <- "TWN"

# x <- data.frame(gdp=unlist(gdp[, grep('X', names(gdp))]), 
#     year=as.numeric(gsub('X', '', rep(yrs, each=nrow(gdp)))),
#     ccode=rep(gdp$ccode, times=length(yrs)),
#     iso3=rep(gdp$iso3, times=length(yrs)))
# x <- x[!is.na(x$iso3), ]
# dim(x)
# write.csv(x, '~/dropbox/cliodata/cliogdp.csv', row.names=F)

x <- read.csv('~/dropbox/cliodata/cliogdp.csv')
x$century <- trunc(x$year / 100)
x$century[x$century >= 10 & x$century <=15] <- "10-5"
x$century <- paste0('c', x$century)

meangdp <- aggregate(gdp ~ century + iso3, data=x, mean)