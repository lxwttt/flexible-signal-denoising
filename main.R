library(ashr)
source("./utils.R")

smash.mu <- function(wc, data.var, tsum, ashparam, J, n) {
    wmean <- matrix(0, J, n)
    wvar <- matrix(0, J, n)
    y <- wc
    vtable <- sdtable(data.var)$sumtable
    for (j in 0:(J - 1)) {
        ind.nnull <- (vtable[j + 2, ] != 0)
        zdat.ash <- call.ash(y[j + 2, ind.nnull], sqrt(vtable[
            j + 2,
            ind.nnull
        ]), ashparam,
        df = NULL, SGD = FALSE
        )
        wmean[j + 1, ind.nnull] <- get_pm(zdat.ash) / 2
        wmean[j + 1, !ind.nnull] <- 0
    }
    mu.est <- rev.wave(tsum, wmean, -wmean)
    return(mu.est)
}


smash.var <- function(data, data.var, x.var.ini, ashparam,
                      weight, J, n) {
    wmean <- matrix(0, J, n)
    wvar <- matrix(0, J, n)
    vtable <- sdtable(data.var)$sumtable
    vdtable <- sdtable(data)$difftable
    for (j in 0:(J - 1)) {
        ind.nnull <- (vtable[j + 2, ] != 0)
        zdat.ash <- call.ash(vdtable[j + 2, ind.nnull],
            sqrt(vtable[j + 2, ind.nnull]),
            ashparam,
            df = min(50, 2^(j + 1)), SGD = TRUE
        )
        wmean[j + 1, ind.nnull] <- get_pm(zdat.ash) / 2
        wmean[j + 1, !ind.nnull] <- 0
    }
    var.est <- rev.wave(weight * sum(data) + (1 - weight) *
        sum(x.var.ini), wmean, -wmean)
    return(var.est)
}


call.ash <- function(x, x.sd, ashparam, df, SGD) {
    zdat.ash <- withCallingHandlers(do.call(
        ash,
        c(list(betahat = x, sebetahat = x.sd), ashparam)
    ))
    return(zdat.ash)
}


sdtable <- function(signal) {
    n <- length(signal)
    J <- as.integer(log2(n))

    tempvec <- rep(0, 2 * n)
    tempvec2 <- rep(0, 2 * n)
    sumtable <- matrix(0, nrow = J + 1, ncol = n)
    difftable <- matrix(0, nrow = J + 1, ncol = n)
    sumtable[1, ] <- signal
    difftable[1, ] <- signal

    for (D in 0:(J - 1)) {
        nD <- 2^(J - D)
        pD <- 2^D
        for (l in 0:(pD - 1)) {
            a <- l * nD + 1
            tempvec[1:(nD - 1)] <- sumtable[D + 1, a:(a + nD - 2)]
            tempvec[nD] <- sumtable[D + 1, a + nD - 1]
            tempvec[nD + 1] <- sumtable[D + 1, a + nD - 1]
            tempvec[(nD + 2):(2 * nD)] <- sumtable[D + 1, a:(a + nD - 2)]

            for (i in 0:(nD - 1)) {
                sumtable[D + 2, a + i] <- tempvec[2 * i + 1] + tempvec[2 * i + 2]
                difftable[D + 2, a + i] <- tempvec[2 * i + 1] - tempvec[2 * i + 2]
            }
        }
    }

    return(list(sumtable = sumtable, difftable = difftable))
}

rev.wave <- function(estimate, pmat, qmat) {
    J <- nrow(pmat)
    np <- ncol(pmat)
    est <- rep(estimate[1], np)

    for (D in 0:(J - 1)) {
        nD <- 2^(D + 1)
        pD <- 2^(J - D - 1)
        nDo2 <- nD / 2
        tempvecl <- numeric(nD)
        tempvecr <- numeric(nD)

        for (l in 0:(pD - 1)) {
            a <- l * nD + 1
            for (i in 0:(nDo2 - 1)) {
                dep <- est[a + i] / 2
                dp <- pmat[J - D, a + i]
                dq <- qmat[J - D, a + i]
                tempvecl[2 * i + 1] <- dep + dp
                tempvecl[2 * i + 2] <- dep + dq
            }
            if (nDo2 > 1) {
                for (i in nDo2:(nD - 2)) {
                    dep <- est[a + i + 1] / 2
                    dp <- pmat[J - D, a + i + 1]
                    deq <- est[a + i] / 2
                    dq <- qmat[J - D, a + i]
                    tempvecr[2 * (i - nDo2) + 1] <- deq + dq
                    tempvecr[2 * (i - nDo2) + 2] <- dep + dp
                }
            }
            dep <- est[a + nDo2] / 2
            dp <- pmat[J - D, a + nDo2]
            deq <- est[a + nD - 1] / 2
            dq <- qmat[J - D, a + nD - 1]
            tempvecr[nD - 1] <- deq + dq
            tempvecr[nD] <- dep + dp

            for (i in 0:(nD - 1)) {
                est[a + i] <- 0.5 * (tempvecl[i + 1] + tempvecr[i + 1])
            }
        }
    }
    return(est)
}

smash.gaus <- function(x, sigma = NULL, weight = 0.5, ashparam = list()) {
    if (!ispowerof2(length(x))) {
        reflect.res <- reflect(x)
        idx <- reflect.res$idx
        x <- reflect.res$x
    } else {
        idx <- 1:length(x)
    }
    J <- log2(length(x))
    n <- length(x)

    ashparam <- list(
        pointmass = TRUE, prior = "nullbiased",
        gridmult = 2, mixcompdist = "normal",
        nullweight = 10, outputlevel = 2, fixg = FALSE
    )


    tsum <- sum(x)
    x.w.d <- sdtable(x)$difftable


    var.est1.ini <- (x - lshift(x))^2 / 2
    var.est2.ini <- (rshift(x) - x)^2 / 2
    var.est.ini <- (var.est1.ini + var.est2.ini) / 2
    mu.est <- smash.mu(
        x.w.d, var.est.ini, tsum, ashparam, J, n
    )
    var.est <- (x - mu.est)^2
    var.var.est <- 2 / 3 * var.est^2
    var.est <- smash.var(
        var.est, var.var.est, var.est.ini,
        ashparam, weight, J, n
    )
    var.est[var.est <= 0] <- 1e-08
    sigma <- sqrt(var.est)


    ashparam.mean <- ashparam
    ashparam.mean$gridmult <- 64
    mu.res <- smash.mu(
        x.w.d, sigma^2, tsum, ashparam.mean, J, n
    )

    return(mu.res[idx])
}
