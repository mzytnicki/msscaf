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

normalizeKR <- function (data) {
    mat <- makeFullMatrix(data)
    n   <- nrow(mat)
    nullRows <- which((colSums(mat) == 0) | (rowSums(mat) == 0))
    if (length(nullRows) > 0) {
        message(paste0("    ", length(nullRows), " rows/columns are empty."))
    }
    diag(mat)[nullRows] <- 1
    matKR <- KR(mat)
    matKR[nullRows, ] <- 0
    matKR[, nullRows] <- 0
    vecKR <- as.vector(t(matKR))
    vecKR[is.na(vecKR)] <- 0
    return(makeTibbleFromList(vecKR, n))
}

normalizeMD <- function (data) {
    data %<>%
        makeFullTibble() %>%
        mutate(distance = abs(bin1 - bin2))
    sampled <- data %>%
        sample_n(size = min(object@sampleSize, nrow(data))) %>%
        rename(sampledDistance = distance) %>%
        select(c(sampledDistance, count)) %>%
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
        mutate(loess = predict(l)) %>%
        mutate(loess = pmax(loess, 0)) %>%
        select(-count) %>%
        unique()
    sampledDistances <- unique(sort(sampled$sampledDistance))
    uniqueDistances <- unique(sort(data$distance))
    valueMap <- tibble(distance = uniqueDistances,
                       sampledDistance =
                           sapply(uniqueDistances, function(x) {
                               sampledDistances[which.min(abs(x - sampledDistances))]
                           })) %>%
        left_join(sampled, by = "sampledDistance") %>%
        select(-sampledDistance)
    data %<>% left_join(valueMap, by = "distance") %>%
        mutate(count = log2((count + 0.0001) / (loess + 0.0001))) %>%
        select(-c(distance, loess)) %>%
        filter(bin1 <= bin2) %>%
        filter(count != 0)
    return(data)
}

normalize <- function(object) {
    message(paste0("  Normalizing ", object@chromosome, "."))
    object@interactionMatrix <- normalizeKR(object@interactionMatrix)
    #plot.10XRef(object)
    object@interactionMatrix <- normalizeMD(object@interactionMatrix)
    #plot.10XRef(object)
    return(object)
}

randomizeData <- function (data) {
    data %>% mutate(count = sample(count))
}

computeMeanRectangle <- function(value, data = data) {
    data %>%
        filter(bin1 <= value) %>%
        filter(bin2 >= value) %>%
        summarise(nCells = n(), meanCount = mean(count))
}

computeMeanRectangles <- function(data) {
    n <- max(data$bin1, data$bin2)
    map_df(seq(n), computeMeanRectangle, data = data) %>%
        rowid_to_column(var = "bin")
}

checkBreak <- function(object, meanCountThreshold, nCellsThreshold) {
    message(paste0("  Checking ", object@chromosome, "."))
    data       <- object@interactionMatrix
    dataRandom <- randomizeData(data)
    rectangles <- computeMeanRectangles(data) %>%
        add_column(randMeanCount = computeMeanRectangles(dataRandom)$meanCount) %>%
        mutate(fcMeanCount = meanCount - randMeanCount)
    rectangles %>%
        gather(key = "type", value = "value", fcMeanCount, nCells) %>%
        ggplot(aes(x = bin, y = value)) +
            geom_line() +
            facet_grid(rows = vars(type), scales = "free_y")
    breakPoint <- rectangles %>%
        arrange(fcMeanCount, desc(nCells)) %>%
        slice(1) %>%
        select(bin, fcMeanCount, nCells) %>%
        as.list()
    if (breakPoint$fcMeanCount <= meanCountThreshold &
        breakPoint$nCells >= nCellsThreshold) {
        return(list(ref = object@chromosome, bin = breakPoint$bin))
    }
    return(list(ref = object@chromosome, bin = -1))
}

normalizeAndBreak <- function(data, threshold = threshold, nCells = nCells) {
    data <- normalize(data)
    data <- checkBreak(data, threshold, nCells)
}

checkBreaks <- function (object) {
    message("Splitting matrix.")
    data <- splitByRef(object@interactionMatrix)
    output <- lapply(data,
                     normalizeAndBreak,
                     threshold = object@breakThreshold,
                     nCells    = object@breakNCells)
}
