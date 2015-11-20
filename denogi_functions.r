tablefix <- function(dat){
    fill <- NULL
    for (century in unique(centuries)){
        temp <- dat[, c(1, which(centuries %in% century) + 1)]
        names(temp) <- variables[1:4]
        temp$century <- century
        fill <- rbind(fill, temp)
    }
    return(fill)    
}

FA_inits <- function(y, vrs){
    y$init <- NA
    for (i in 3:length(vrs)){ # min 3 for fantanal
        combinations <- combn(vrs, i)
        for (j in 1:ncol(combinations)){
            variables <- combinations[, j]
            dat_fa <- y[complete.cases(y[, variables]), variables]
            try(fa <- factanal(dat_fa, 1, scores='Bartlet'))
            y$init[match(rownames(fa$scores), rownames(y))] <- fa$scores
        }
    }
    initxi <- y$init
    initgamma <- t(apply(y[, variables], 2, function(x) lm(x ~ y$init)$coef))
    initomega <- apply(y[, variables], 2, function(x) summary(lm(x ~ y$init))$sigma)
    return(list(xi=initxi, gamma=initgamma, omega=initomega))
}

extract_from_coda <- function(coda, varname){
    keep <- grep(varname, colnames(coda[[1]]))
    # out <- coda[[1]][, keep]
    out <- coda[, keep, drop=FALSE]
    out <- as.matrix(out)
    return(out)
}

sumstats <- function(x){
    out <- c(mean(x), quantile(x, c(0.05, 0.5, 0.95)))
    names(out) <- c('mean', 'q05', 'q50', 'q95')
    return(out)
}

sumstatsDF <- function(coda){
    out <- apply(coda, 2, sumstats)
    out <- t(out)
    out <- as.data.frame(out)
    return(out)
}

vstrsplit <- function(strs, split){
    spltlist <- strsplit(strs, split)
    spltvec <- do.call(rbind, spltlist)
    return(spltvec)
}

place_text <- function(d, countries, col=1, cex=0.6, clio=F){
    # d should be an ordered dataframe
    N <- nrow(d)

    y_left <- seq(from=1, to=N, by=2)
    x_left <- d$q05[y_left] 
    text_left <- countries[y_left]
    text(x_left, y_left, text_left, cex=cex, pos=2, col=col)

    y_right <- seq(from=2, to=N, by=2)
    x_right <- d$q95[y_right]
    text_right <- countries[y_right]
    text(x_right, y_right, text_right, cex=cex, pos=4, col=col)
}