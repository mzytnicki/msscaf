.notSquareICE <- function(A) {
    #message(str(summary(A)))
    s               <- rowSums(A)
    sm              <- mean(s)
    bias            <- s / sm
    bias[bias == 0] <- 1.0
    A               <- A / bias
    #message(str(bias))
    #message(str(summary(A)))
    return(A)
}

notSquareICE <- function(A, nIter = 1000, maxDelta = 0.1) {
    Asav <- A
    for (iter in seq.int(nIter)) {
        A     <- .notSquareICE(A)
        A     <- t(.notSquareICE(t(A)))
        s     <- colSums(A)
        sm    <- mean(s)
        delta <- max(abs(s / sm - 1))
        if (any(is.na(delta))) {
          message("Warning, algorithm diverges")
          return(Asav)
        }
        if (delta < maxDelta) {
            return(A)
        }
    }
    message(paste0("Warning, algorithm did not converge in ", nIter, " iterations."))
    return(Asav)
}

normalizeNotSquareICE <- function(object) {
    data <- object@interactionMatrix
    if (nrow(data) == 0) {
        # Matrix is empty, do not normalize it
        return(object)
    }
    mat <- sparseMatrix(data$bin1 + 1, data$bin2 + 1, x = data$count)
    mat <- summary(notSquareICE(mat))
    object@interactionMatrix <- tibble::tibble(bin1  = mat$i - 1,
                                       bin2  = mat$j - 1,
                                       count = mat$x)
    return(object)
}

ICE <- function(A, nIter = 10000, maxDelta = 0.1) {
    for (iter in seq.int(nIter)) {
        s       <- colSums(A)
        sm      <- mean(s)
        if (log10(sm) > 300) {
          stop("Warning, algorithm diverges")
        }
        bias    <- s / sm
        A       <- A / bias
        A       <- A %*% Diagonal(x = 1 / bias)
        #message(paste0("max(s): ", max(s), ", min(sm): ", min(sm)))
        delta   <- max(abs(s / sm - 1))
        #message(paste0("delta: ", delta))
        if (delta < maxDelta) {
            return(A)
        }
    }
    message(paste0("Warning, algorithm did not converge in ", nIter, " iterations."))
    return(A)
}


KR <- function(A, tol = 1e-6, delta = 0.1, Delta = 3) {
    n <- nrow(A)
    e <- matrix(1, nrow = n, ncol = 1)
    x0 <- e
    # inner stopping criterior
    g <- 0.9
    etamax <- 0.1
    eta <- etamax
    stop_tol <- tol * .5
    x <- x0
    rt <- tol^2
    v <- x * (A %*% x)
    rk <- 1 - v
    rho_km1 <- drop(t(rk) %*% rk)
    rout <- rho_km1
    rold <- rout
    MVP <- 0 # We'll count matrix vector products.
    i <- 0 # Outer iteration count.
    while (rout > rt) { # Outer iteration
        i <- i + 1
        k <- 0
        y <- e
        innertol <- max(c(eta^2 * rout, rt))
        while (rho_km1 > innertol) { #Inner iteration by CG
            k <- k + 1
            if (k == 1) {
                Z <- rk / v
                p <- Z
                rho_km1 <- drop(t(rk) %*% Z)
            }
            else {
                beta <- rho_km1 / rho_km2
                p <- Z + beta * p
            }
            
            # update search direction efficiently
            w <- x * (A %*% (x * p)) + v * p
            alpha <- rho_km1 / drop(t(p) %*% w)
            ap <- alpha * p
            # test distance to boundary of cone
            ynew <- y + ap
            if (min(ynew) <= delta) {
                if (delta == 0) break()
                ind <- which(ap < 0)
                gamma <- min((delta - y[ind]) / ap[ind])
                y <- y + gamma * ap
                break()
            }
            if (max(ynew) >= Delta) {
                ind <- which(ynew > Delta)
                gamma <- min((Delta - y[ind]) / ap[ind])
                y <- y + gamma * ap
                break()
            }
            y <- ynew
            rk <- rk - alpha * w
            rho_km2 <- rho_km1
            Z <- rk / v
            rho_km1 <- drop(t(rk) %*% Z)
        }
        x <- x * y
        if (log10(min(x)) < -300) {
          stop("Warning, algorithm diverges")
          result <- t(t(x[,1] * A) * x[,1])
          return(result)
        }
        v <- x * (A %*% x)
        rk <- 1 - v
        rho_km1 <- drop(t(rk) %*% rk)
        rout <- rho_km1
        MVP <- MVP + k + 1
        # Update inner iteration stopping criterion.
        rat <- rout / rold
        rold <- rout
        res_norm <- sqrt(rout)
        eta_o <- eta
        eta <- g * rat
        if (g * eta_o^2 > 0.1) {
            eta <- max(c(eta, g * eta_o^2))
        }
        eta <- max(c(min(c(eta, etamax)), stop_tol / res_norm))
    }
    
    result <- t(t(x[,1] * A) * x[,1])
    return(result)
}

normalizeKR <- function(object, sizes = NULL) {
    if (is(object, "msscafExp")) {
        chromosome <- FALSE
    }
    else if (is(object, "msscafRefExp")) {
        chromosome <- TRUE
    }
    else {
        stop(paste0("Function 'normalizeKR' does not know what to do with parameter of type '", is(object), "'."))
    }
    #message("    KR normalization.")
    data <- object@interactionMatrix
    if (nrow(data) == 0) {
        # Matrix is empty, do not normalize it
        return(object)
    }
    if (chromosome) {
        mat <- makeFullMatrix(data)
    }
    else {
        mat <- makeFullMatrixGenome(data, sizes)
    }
    n   <- nrow(mat)
    nullRows <- which((colSums(mat) == 0) | (rowSums(mat) == 0))
    if (length(nullRows) > 0) {
        # message(paste0("      ", length(nullRows), " rows/columns are empty."))
        diag(mat)[nullRows] <- 1
    }
    mat <- tryCatch(
        expr = {
             KR(mat)
        },
        error = function(e) {
            # message("KR did not converge, resorting to ICE normalization.")
            mat <- tryCatch(
                expr = {
                     ICE(mat)
                },
                error = function(e) {
                    message("Neither KR nor ICE converged.")
                    return(object)
                })
        })
    if (length(nullRows) > 0) {
        mat[nullRows, ] <- 0
        mat[, nullRows] <- 0
    }
    if (chromosome) {
        object@interactionMatrix <- makeSparseMatrix(mat)
    }
    else {
        object@interactionMatrix <- makeSparseMatrixGenome(mat, sizes)
    }
    return(object)
}

normalizeMD <- function(object) {
    #message("    MD normalization.")
    data <- object@interactionMatrix
    if (nrow(data) == 0) {
        # Matrix is empty, do not normalize it
        return(object)
    }
    nCellsNotDiag <- data %>%
        filter(bin1 != bin2) %>%
        nrow()
    if (nCellsNotDiag <= 2 * object@size) {
        # Matrix is too sparse, everything is on diagonal
        data %<>% filter(bin1 != bin2)
        medianValue <- median(data$count)
        data %<>% 
            mutate(count = log2((count + 0.0001) / (medianValue + 0.0001))) %>%
            filter(count != 0)
        object@interactionMatrix <- data
        return(object)
    }
    data %<>%
        #makeFullTibble() %>%
        mutate(distance = abs(bin1 - bin2))
    if (max(data$distance) == 0) {
        # Everything is on the diagonal, skip
        return(object)
    }
    sampled <- data %>%
        sample_n(size = min(object@parameters@sampleSize, nrow(data))) %>%
        rename(sampledDistance = distance) %>%
        dplyr::select(c(sampledDistance, count)) %>%
        arrange(sampledDistance)
    
    optimizeSpan <- function(model, spans = c(0.01, 0.9)) {
        getLoessCriterion <- function(x) {
            traceL <- x$trace.hat
            sigma2 <- sum(x$residuals^2)/(x$n - 1)
            gcv    <- x$n * sigma2/(x$n - traceL)^2
            return(gcv)
        }
        
        optimizationFunction <- function(model = model, span) {
            updatedModel <- update(model, span = span)
            getLoessCriterion(updatedModel)
        }
        
        result <- optimize(optimizationFunction, model = model, spans)
        return(result$minimum)
    }
    
    l <- loess(count ~ sampledDistance, data = sampled)
    span <- optimizeSpan(l)
    l <- loess(count ~ sampledDistance, span = span, data = sampled)
    sampled %<>%
        dplyr::mutate(loess = predict(l)) %>%
        dplyr::mutate(loess = pmax(loess, 0)) %>%
        dplyr::select(-count) %>%
        unique()
    sampledDistances <- unique(sort(sampled$sampledDistance))
    uniqueDistances <- unique(sort(data$distance))
    valueMap <- tibble::tibble(distance = uniqueDistances,
                       sampledDistance =
                           sapply(uniqueDistances, function(x) {
                               sampledDistances[which.min(abs(x - sampledDistances))]
                           })) %>%
        purrr::left_join(sampled, by = "sampledDistance") %>%
        dplyr::select(-sampledDistance)
    data %<>% purrr::left_join(valueMap, by = "distance") %>%
        dplyr::mutate(count = log2((count + 0.0001) / (loess + 0.0001))) %>%
        dplyr::select(-c(distance, loess)) %>%
        dplyr::filter(count != 0)
    object@interactionMatrix <- data
    return(object)
}

.normalizeHighCountRows <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    message(paste0("\tDataset ", object@name))
    normalizeHighCountRowsCpp(object@interactionMatrix, sizes)
}

normalizeHighCountRows <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Trimming high count lines.")
    purrr::walk(object@data, .normalizeHighCountRows, sizes = object@sizes)
}


normalize <- function(object) {
    object <- normalizeKR(object)
    #plot.10XRef(object)
    object <- normalizeMD(object)
    #plot.10XRef(object)
    return(object)
}
