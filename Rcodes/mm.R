anova.asreml <-
function (x, ...) 
{
    UseMethod("wald")
}
asreml <-
function (fixed = y ~ 1, random, sparse, rcov, G.param, R.param, 
    predict = predict.asreml(), constraints = asreml.constraints(), 
    data = sys.parent(), subset, family = asreml.gaussian(), 
    weights = NULL, offset = NULL, na.method.Y = "include", na.method.X = "fail", 
    keep.order = FALSE, ran.order = "user", fixgammas = FALSE, 
    as.multivariate = NULL, model.frame = FALSE, start.values = FALSE, 
    update.Gcon = TRUE, update.Rcon = TRUE, dump.model = FALSE, 
    model = FALSE, debug = FALSE, control = asreml.control(...), 
    ...) 
{
    if (!asreml.Rsys) 
        ObjectSize <- options()$object.size
    dpmv <- -1e-37
    if (mode(model) == "logical" && model == TRUE) 
        stop("\nMust supply an object to argument model\n")
    if (mode(model) != "logical") {
        ws <- max(model$control$workspace, model$control$pworkspace)
        if (!asreml.Rsys) 
            if (object.size() < ws * 8) 
                options(object.size = ws * 8)
        asr.fit <- asreml.call(model$data, model$inter, model$struc, 
            model$glm, model$asr.predict$predict, model$asr.predict, 
            model$weights, model$offset, model$multivariate, 
            model$ntrt, model$control, model$keep, model$constraints, 
            dpmv, model$debug)
        asr.fit$call <- model$call
        asr.fit$distribution <- model$glm$id
        asr.fit$link <- model$glm$link
        asr.fit$G.param <- asreml.gupdt(model$G.param, asr.fit$gammas, 
            names(asr.fit$gammas.con), VCov = TRUE)
        asr.fit$R.param <- asreml.rupdt(model$R.param, asr.fit$gammas, 
            names(asr.fit$gammas.con), VCov = TRUE)
        asr.fit$factor.names <- model$asr.inter$facnam[model$asr.inter$neword]
        asr.fit$fixed.formula <- model$fixed
        asr.fit$random.formula <- model$random.Inter
        asr.fit$sparse.formula <- model$sparse
        if (model.frame) 
            asr.fit$model.frame <- model$data
        oldClass(asr.fit) <- "asreml"
        if (!asreml.Rsys) 
            options(object.size = max(ObjectSize, object.size(asr.fit)))
        return(asr.fit)
    }
    ll <- c(as.list(control$call)[-1], list(...))
    control <- do.call("asreml.control", ll[!duplicated(names(ll))])
    if (debug) {
        cat("   Checking args... ")
        ptime <- proc.time()[1]
    }
    call <- match.call()
    if (!inherits(fixed, "formula")) 
        stop("\nfixed must be a formula")
    if (length(fixed) != 3) 
        stop("\nFixed model formula must be of the form \"resp ~ pred\"")
    if (missing(random)) {
        random <- ~NULL
    }
    else {
        if (!inherits(random, "formula")) 
            stop("\nrandom must be a formula")
        if (length(random) != 2) 
            stop("\nRandom model formula must be of form \" ~ pred\"")
    }
    if (missing(sparse)) {
        sparse <- ~NULL
    }
    else {
        if (!inherits(sparse, "formula")) 
            stop("\nsparse must be a formula")
        if (length(sparse) != 2) 
            stop("\nSparse model formula must be of form \" ~ pred\"")
    }
    if (missingRcov <- missing(rcov)) {
        rcov <- ~NULL
    }
    else {
        if (!inherits(rcov, "formula")) 
            stop("\nrcov must be a formula")
        if (length(rcov) != 2) 
            stop("\nRcov model formula must be of the form \" ~ pred\"")
    }
    if (missing(data)) 
        data <- environment(fixed)
    varNames <- names(data)
    groups <- control$group
    if (length(groups) > 0) {
        for (i in names(groups)) {
            if (is.numeric(groups[[i]])) 
                groups[[i]] <- varNames[groups[[i]]]
            else {
                if (any(is.na(match(groups[[i]], varNames)))) 
                  stop(paste("Object in group", names(groups[1]), 
                    "not in data."))
            }
        }
    }
    mbfList <- list()
    if (length(control$mbf) > 0) {
        mbfList <- asreml.getMbf(control$mbf, data)
    }
    if (!missing(weights)) 
        weights <- as.character(substitute(weights))
    if (!missing(offset)) 
        offset <- as.character(substitute(offset))
    if (!missing(as.multivariate)) 
        as.multivariate <- as.character(substitute(as.multivariate))
    form <- asreml.ModelFormula(fixed, random, sparse, rcov, 
        weights, offset, groups, mbfList$mbfAttr, ignore = c("trait", 
            "units"), varNames)
    model.y <- attr(form, "model.y")
    noDatFr <- (missing(data) || asreml.ie(length(model.y) > 
        1, is.null(data[model.y]), is.null(data[[model.y]])))
    if (length(model.y) == 1 && noDatFr && is.matrix(eval(as.name(model.y)))) 
        model.y <- dimnames(eval(as.name(model.y)))[[2]]
    multivariate <- (length(model.y) > 1 || !is.null(as.multivariate))
    if (missing(subset)) 
        mf <- as.call(list(as.name("asreml.mf.default"), formula = form, 
            data = data, drop.unused.levels = control$drop.unused.levels, 
            na.action = "na.pass"))
    else mf <- as.call(list(as.name("asreml.mf.default"), formula = form, 
        data = call$data, subset = call$subset, drop.unused.levels = control$drop.unused.levels, 
        na.action = "na.pass"))
    mf[[1]] <- as.name("asreml.mf.default")
    attr(mf$formula, ".Environment") <- sys.frame(sys.parent())
    data <- eval(mf)
    if (missingRcov) 
        rcov <- asreml.formula(parse(text = "~ units"))
    family <- asreml.getFamily(family)
    asr.glm <- asreml.glm(family)
    asrAttrib <- list(Control = control, GROUP = groups, MBF = mbfList)
    for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    data <- asreml.data(data, model.y, fixed, random, sparse, 
        rcov, asr.glm, asrAttrib, na.method.X, na.method.Y, as.multivariate, 
        weights, offset, ran.order)
    thrlevels <- attr(data, "thrlevels")
    S2 <- attr(data, "S2")
    keep <- attr(data, "keep")
    ntrt <- attr(data, "ntrt")
    model.y <- attr(data, "model.y")
    fixed <- attr(data, "fixed")
    sparse <- attr(data, "sparse")
    random <- attr(data, "random")
    random.Inter <- attr(data, "random.Inter")
    asr.glm$thrlevels <- thrlevels
    asr.glm$kwtcol <- ifelse(asr.glm$id == 9, 2, 1)
    if (ntrt == 1 | asr.glm$id == 9) 
        S2 <- asr.glm$dispersion
    if (is.null(weights)) 
        wts <- NULL
    else wts <- data[, weights]
    asr.ainv <- asreml.ainv(wts, asr.glm, attr(data, "nRow"), 
        control$ginverse)
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Building interaction table... ")
        ptime <- proc.time()[1]
    }
    if ((length(control$Csparse[[2]]) > 0) & !inherits(control$Csparse, 
        "formula")) 
        stop("Csparse must be a formula")
    asr.inter <- asreml.inter(data, fixed, random.Inter, sparse, 
        control$Csparse, model.y, keep.order, asr.ainv$locgiv, 
        asr.glm, multivariate)
    control$Ftest <- asreml.Fown(control$Ftest, asr.inter)
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Generating RCOV list... ")
        ptime <- proc.time()[1]
    }
    if (missing(R.param)) {
        R.param <- asreml.rdflt(y = model.y, rcov = rcov, data = data, 
            dispersion = S2)
    }
    else {
        Rdflt <- asreml.rdflt(y = model.y, rcov = rcov, data = data, 
            dispersion = S2)
        R.param <- asreml.setR(Rdflt, R.param, update.Rcon)
        if (length(asr.glm$dispersion) > 1) 
            stop("\nDispersion parameter must be a scalar\n")
        R.param <- lapply(R.param, function(x, a) {
            nam <- names(x$variance$s2)
            if (!is.na(a)) {
                x$variance$s2 <- a
                x$variance$con <- "F"
            }
            names(x$variance$s2) <- nam
            x
        }, asr.glm$dispersion)
    }
    if (length(R.param) > 1 & length(model.y) == 1) 
        scale <- var(data[[model.y]][!is.na(data[[model.y]])])/2
    else scale <- 1
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Generating GCOV list... ")
        ptime <- proc.time()[1]
    }
    if (missing(G.param)) 
        G.param <- asreml.gdflt(random, data, scale = scale, 
            control = control)
    else {
        Gdflt <- asreml.gdflt(random, data, scale = scale, control = control)
        G.param <- asreml.setG(Gdflt, G.param, update.Gcon)
    }
    if (!missing(start.values)) {
        if ((is.logical(start.values) & start.values == TRUE) | 
            is.character(start.values)) {
            sv.list <- list(G.param = G.param, R.param = R.param, 
                gammas.table = asreml.gammas(list(G.param = G.param, 
                  R.param = R.param)))
            if (!asreml.Rsys) 
                options(object.size = max(ObjectSize, object.size(sv.list)))
        }
        if (is.character(start.values)) {
            if (asreml.Rsys) 
                write.table(sv.list$gammas.table, start.values, 
                  row.names = FALSE, sep = ",")
            else write.table(sv.list$gammas.table, start.values, 
                dimnames.write = "col", sep = ",", quote.strings = TRUE)
        }
        if (debug) 
            cat((proc.time()[1] - ptime), "seconds\n")
        return(sv.list)
    }
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Building variance structure table... ")
        ptime <- proc.time()[1]
    }
    asr.struc <- asreml.struc(data, asr.inter, random, G.param, 
        rcov, R.param, fixgammas, asr.ainv)
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Label solutions, set up call...\n")
        ptime <- proc.time()[1]
    }
    asr.predict <- asreml.prdStruc(predict, data, asr.inter, 
        asr.struc, asr.glm)
    if (multivariate) 
        ntrt <- 1
    if (dump.model) {
        model.list <- list(data = data, inter = asr.inter, struc = asr.struc, 
            glm = asr.glm, weights = weights, offset = offset, 
            multivariate = multivariate, ntrt = ntrt, control = control, 
            G.param = G.param, R.param = R.param, call = call, 
            keep = keep, constraints = constraints, predict = predict, 
            asr.predict = asr.predict, fixed = fixed, random.Inter = random.Inter, 
            sparse = sparse, debug = debug)
        if (!asreml.Rsys) 
            options(object.size = max(ObjectSize, object.size(model.list)))
        return(model.list)
    }
    if (!asreml.Rsys) 
        if (object.size() < ws * 8) 
            options(object.size = ws * 8)
    asr.fit <- asreml.call(data, asr.inter, asr.struc, asr.glm, 
        asr.predict$predict, asr.predict, weights, offset, multivariate, 
        ntrt, control, keep, constraints, dpmv, debug)
    if (debug) {
        cat((proc.time()[1] - ptime), "seconds\n")
        cat("   Done!\n")
    }
    asr.fit$call <- call
    asr.fit$distribution <- asr.glm$id
    asr.fit$link <- asr.glm$link
    asr.fit$family <- family
    if (!is.null(weights)) {
        asr.fit$prior.weights <- rep(0, length(keep))
        asr.fit$prior.weights[keep] <- data[, weights]
    }
    asr.fit$fitted.values <- family$inverse(asr.fit$linear.predictors)
    asr.fit$control <- control
    if (debug) {
        asr.fit$G.param <- G.param
        asr.fit$R.param <- R.param
    }
    else {
        asr.fit$G.param <- asreml.gupdt(G.param, asr.fit$gammas, 
            names(asr.fit$gammas.con), VCov = TRUE)
        asr.fit$R.param <- asreml.rupdt(R.param, asr.fit$gammas, 
            names(asr.fit$gammas.con), VCov = TRUE)
    }
    asr.fit$factor.names <- asr.inter$facnam[asr.inter$neword]
    asr.fit$fixed.formula <- fixed
    asr.fit$random.formula <- random.Inter
    asr.fit$sparse.formula <- sparse
    if (model.frame) 
        asr.fit$model.frame <- data
    if (debug) 
        return(list(control = control, asr.inter = asr.inter, 
            asr.struc = asr.struc, asr.glm = asr.glm, asr.predict = asreml.ie(is.null(predict), 
                NULL, asr.predict), data = data, asr.fit = asr.fit, 
            fixed.formula = fixed, random.formula = random.Inter, 
            sparse.formula = sparse))
    class(asr.fit) <- "asreml"
    if (!asreml.Rsys) 
        options(object.size = max(ObjectSize, object.size(asr.fit)))
    asr.fit
}
asreml.About <-
function () 
{
    where <- paste(asreml.findMe(), collapse = " ")
    pp <- packageDescription("asreml")
    package <- unclass(unlist(as.POSIXlt(strptime(pp$Packaged, 
        "%Y-%m-%d"))))
    msg <- paste(rep(" ", 1024), collapse = "")
    arch <- R.Version()$platform
    what <- .C("about", where, msg = msg, arch, pp$Version, as.integer(package["mday"]), 
        as.integer(package["mon"]), as.integer(package["year"] + 
            1900))
    cat(what$msg)
    invisible(what$msg)
}
asreml.aexp <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "aexp"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ai"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.aexph <-
function (x, y, init = NA, data, ...) 
{
    fun <- "aexph"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ah"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.aexpv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "aexpv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "av"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.agau <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "agau"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ai"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.agauh <-
function (x, y, init = NA, data, ...) 
{
    fun <- "agauh"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ah"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.agauv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "agauv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "av"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.ainv <-
function (weight, glm, nrow, ginverse) 
{
    base <- 1
    nainv <- 0
    locwts <- 0
    locgiv <- vector(length = 0)
    lcainv <- matrix(0, nrow = 5, ncol = length(ginverse$ainverse) + 
        as.numeric(!is.null(weight)) + as.numeric(glm$id > 1))
    ainv <- vector(length = 0)
    llainv <- vector(length = 0, mode = "numeric")
    aipt <- 0
    llpt <- base
    dpmv <- -1e-37
    if (!is.null(weight)) {
        nainv <- 1
        aipt <- aipt + 2
        if (glm$id == 9) {
            ainv <- c(ainv, 0, as.vector(rbind(rep(0, length(weight)), 
                weight))[-1])
            lcainv[, nainv] <- c(nrow, 3, aipt, glm$kwtcol, glm$thrlevels)
        }
        else {
            ainv <- c(ainv, 0, weight)
            lcainv[, nainv] <- c(nrow, 8, aipt, llpt, 1)
        }
        locwts <- aipt
        aipt <- length(ainv)
    }
    else if (glm$id > 1) {
        nainv <- 1
        if (glm$id == 9) {
            ainv <- c(ainv, 0, as.vector(rbind(rep(0, nrow), 
                rep(1, nrow)))[-1])
            aipt <- aipt + 2
            lcainv[, nainv] <- c(nrow, 3, aipt, glm$kwtcol, glm$thrlevels)
        }
        else {
            ainv <- c(ainv, 0, rep(1, nrow))
            aipt <- aipt + 2
            lcainv[, nainv] <- c(nrow, 8, aipt, llpt, 1)
        }
        locwts <- aipt
        aipt <- length(ainv)
    }
    if (length(ginverse) > 0) {
        who <- NULL
        where <- NULL
        for (a in 1:length(ginverse$factor)) {
            aa <- names(ginverse$factor)[a]
            if (a == 1 || (a > 1 && !is.element(ginverse$gname(aa, 
                ginverse), who))) {
                nainv <- nainv + 1
                grp <- ginverse$groups(aa, ginverse)
                grpOfst <- ginverse$groupoffset(aa, ginverse)
                locgiv <- c(locgiv, nainv)
                urow <- unique(ginverse$row(aa, ginverse))
                order <- length(urow)
                tt <- table(ginverse$row(aa, ginverse))
                ainv <- c(ainv, 0, ginverse$ai(aa, ginverse))
                aipt <- aipt + 2
                llainv <- c(llainv, cumsum(c(0, tt[-length(tt)])), 
                  ginverse$col(aa, ginverse))
                lcainv[, nainv] <- c(order, 7, aipt, llpt, grp)
                llpt <- llpt + order + length(ginverse$col(aa, 
                  ginverse))
                aipt <- length(ainv)
                who <- c(who, ginverse$gname(aa, ginverse))
                where <- c(where, nainv)
                ainv[lcainv[3, nainv] - 1] <- dpmv
            }
            else locgiv <- c(locgiv, where[match(ginverse$gname(aa, 
                ginverse), who)])
        }
        names(locgiv) <- names(ginverse$factor)
    }
    ainv[is.na(ainv)] <- dpmv
    return(list(nainv = nainv, lcainv = lcainv, ainv = ainv, 
        llainv = llainv, locwts = locwts, locgiv = locgiv))
}
asreml.Ainverse <-
function (pedigree, fgen = list(character(0), 0.01), gender = character(0), 
    groups = 0, groupOffset = 0, method = 0, selfing = NA, inBreed = NA, 
    mgs = FALSE, mv = c("NA", "0", "*"), psort = FALSE) 
{
    newped <- checkPedigree(pedigree, fgen, gender, mv)
    if (psort) 
        return(newped)
    nan <- nrow(newped)
    if ((!is.na(selfing)) & (!is.na(inBreed))) 
        stop("Cannot specify both partial selfing and inbreeding coefficient\n")
    if (is.na(selfing)) 
        selfing <- 0
    ibl <- 1
    if (is.na(inBreed)) {
        inBreed <- 0
        ibl <- 0
    }
    fgsx <- rep(0, nan)
    fmode <- FALSE
    xlink <- FALSE
    if (ncol(newped) == 4) {
        fgsx <- newped[, 4]
        if (attr(newped, "col4Type") == "fgen") 
            fmode <- TRUE
        else if (attr(newped, "col4Type") == "gender") 
            xlink <- TRUE
    }
    lvls <- attr(newped, "rowNames")
    idx <- vector(mode = "numeric", length = nan)
    pmat <- matrix(0, nrow = nan, ncol = 2)
    vped <- match(c(newped[, 1], newped[, 2], newped[, 3]), lvls, 
        nomatch = 0)
    storage.mode(nan) <- "integer"
    storage.mode(vped) <- "integer"
    storage.mode(pmat) <- "integer"
    storage.mode(idx) <- "integer"
    error <- paste(rep(" ", 256), collapse = "")
    pmat.out <- .C("pmat", vped, nan, pmat = pmat, idx = idx, 
        nan, error = error)
    if (nchar(pmat.out$error) > 0) 
        stop(pmat.out$error)
    pmat <- t(pmat.out$pmat)
    lenai <- nan * 5
    llainv <- matrix(0, nrow = 2, ncol = lenai)
    ginv <- vector(mode = "numeric", length = lenai * 3)
    ainv <- vector(mode = "numeric", length = lenai)
    inbreeding <- vector(mode = "numeric", length = nan + 1)
    if (mgs) 
        mgs <- 1
    else mgs <- 0
    linv <- 0
    debug <- 0
    ifault <- 0
    det <- 0
    fmode <- as.numeric(fmode)
    xlink <- as.numeric(xlink)
    storage.mode(fgsx) <- "double"
    storage.mode(fmode) <- "integer"
    storage.mode(xlink) <- "integer"
    storage.mode(llainv) <- "integer"
    storage.mode(ainv) <- "double"
    storage.mode(ginv) <- "double"
    storage.mode(pmat) <- "integer"
    storage.mode(nan) <- "integer"
    storage.mode(lenai) <- "integer"
    storage.mode(linv) <- "integer"
    storage.mode(groups) <- "integer"
    storage.mode(groupOffset) <- "double"
    storage.mode(method) <- "integer"
    storage.mode(selfing) <- "double"
    storage.mode(ibl) <- "integer"
    storage.mode(inBreed) <- "double"
    storage.mode(debug) <- "integer"
    storage.mode(ifault) <- "integer"
    storage.mode(det) <- "double"
    storage.mode(mgs) <- "integer"
    where <- asreml.findMe()
    pdg <- .C("pedigree", llainv = llainv, ainv = ainv, ginv = ginv, 
        pmat, fgsx, inbreeding = inbreeding, nan, lenai, linv = linv, 
        groups, groupOffset, fmode, xlink, method, selfing, ibl, 
        inBreed, mgs, det = det, debug, ifault = ifault, where, 
        error)
    if (pdg$ifault < 0) 
        stop(error)
    ginv = data.frame(matrix(pdg$ginv[1:(3 * pdg$linv)], ncol = 3, 
        byrow = TRUE))
    names(ginv) <- c("Row", "Column", "Ainverse")
    inbreeding = pdg$inbreeding[2:(nan + 1)]
    names(inbreeding) <- lvls[pmat.out$idx]
    attr(ginv, "rowNames") <- lvls[pmat.out$idx]
    attr(ginv, "geneticGroups") <- c(groups, groupOffset)
    list(pedigree = newped, ginv = ginv, inbreeding = inbreeding, 
        det = pdg$det, ifault = pdg$ifault)
}
asreml.ainvPow <-
function (Ainv, Power, Coords, Pstruc) 
{
    sum <- as.vector(matrix(1, nrow = 1, ncol = nrow(Coords)) %*% 
        Coords)
    where <- seq(1, ncol(Coords))[!is.na(sum)]
    newCoords <- matrix(Coords[, where], nrow = nrow(Coords), 
        ncol = length(where), byrow = FALSE)
    temp <- vector(mode = "character", length = ncol(newCoords))
    for (i in 1:ncol(newCoords)) temp[i] <- paste(newCoords[, 
        i], collapse = ",")
    where <- seq(1, ncol(newCoords))[!duplicated(temp)]
    newCoords <- matrix(Coords[, where], nrow = nrow(Coords), 
        ncol = length(where), byrow = FALSE)
    buffer <- vector(mode = "numeric", length = 5)
    buffer[1] <- ncol(newCoords)
    buffer[2] <- 6000 + Pstruc
    buffer[3] <- ifelse(Ainv$nainv == 0, 1, length(Ainv$ainv) + 
        1)
    buffer[4] <- ifelse(Ainv$nainv == 0, 1, length(Ainv$llainv) + 
        1)
    buffer[5] <- 1
    if (Ainv$nainv == 0) {
        Ainv$ainv <- as.vector(newCoords)
        Ainv$llainv <- c(Power, nrow(newCoords), 0, 0, 0)
        Ainv$lcainv <- matrix(buffer, ncol = 1)
    }
    else {
        Ainv$ainv <- c(Ainv$ainv, as.vector(newCoords))
        Ainv$llainv <- c(Ainv$llainv, Power, nrow(newCoords), 
            0, 0, 0)
        Ainv$lcainv <- cbind(Ainv$lcainv, buffer)
    }
    Ainv$nainv <- Ainv$nainv + 1
    return(Ainv)
}
asreml.and <-
function (obj, times = 1, data, ...) 
{
    wgt <- function(times) {
        inter34 <- vector(length = 2)
        x <- -1
        rem <- 1
        while (rem != 0) {
            x <- x + 1
            rem <- times * 10^x - trunc(times * 10^x)
            inter34 <- c(times * 10^x, x)
        }
        inter34
    }
    out <- vector(mode = "list", length = 13)
    sobj <- substitute(obj)
    if (mode(sobj) == "call" && sobj[[1]] != ":" && inherits(obj, 
        "asreml.special")) {
        out[[1]] <- "and"
        out[[2]] <- paste("and(", obj[[2]], asreml.ie(times == 
            1, NULL, paste(",", times)), ")", sep = "")
        if (is.element(obj$Fun, c("ped", "giv", "ide"))) {
            out[[3]] <- obj$Call
            out[[12]] <- obj$Call
        }
        else {
            out[[3]] <- obj[[3]]
            out[[12]] <- obj[[12]]
        }
        out[[4]] <- obj[[4]]
        out[[5]] <- obj[[5]]
        out[[6]] <- obj[[6]]
        out[[7]] <- ""
        out[[8]] <- 2
        out[[8]] <- obj[[8]]
        out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(-12, 0, wgt(times))
        out[[11]] <- obj[[11]]
        out[[13]] <- NA
    }
    else {
        if (mode(sobj) == "call") {
            if (sobj[[1]] == ":") 
                obj <- deparse(sobj, width.cutoff = 500)
            else stop("Bad argument to and()\n")
        }
        else obj <- as.character(substitute(obj))
        tt <- unlist(lapply(strsplit(obj, ":", fixed = TRUE)[[1]], 
            function(x, data) {
                ee <- eval(asreml.Eval(x, data))
                if (inherits(ee, "asreml.special")) 
                  return(ee$Call)
                else return(x)
            }, data))
        out[[1]] <- "and"
        out[[2]] <- paste("and(", paste(tt, collapse = ":"), 
            asreml.ie(times == 1, NULL, paste(",", times)), ")", 
            sep = "")
        out[[3]] <- obj
        out[[4]] <- out[[2]]
        out[[5]] <- NA
        out[[6]] <- ""
        out[[7]] <- ""
        out[[8]] <- 2
        out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(-12, 0, wgt(times))
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- strsplit(obj, ":", fixed = TRUE)[[1]]
        out[[13]] <- NA
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.ante <-
function (obj, k = 1, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ante", "var", numeric(0), "matrix(0.1,nrow=n,ncol=n)+diag(0.05,nrow=n)", 
        k, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 12, 0, 0, k + 1, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0)
    out
}
asreml.ar1 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar1", "cor", 0.1, numeric(0), NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar1h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar1h", "cor", 0.1, "rep(0.1,n)", NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar1v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar1v", "cor", 0.1, 0.1, NA, obj, init, 
        data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
        0)
    out
}
asreml.ar2 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar2", "cor", "rep(0.1,2)", numeric(0), 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar2h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar2h", "cor", "rep(0.1,2)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar2v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar2v", "cor", "rep(0.1,2)", 0.1, NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 
        0)
    out
}
asreml.ar3 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar3", "cor", "rep(0.1,3)", numeric(0), 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar3h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar3h", "cor", "rep(0.1,3)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ar3v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ar3v", "cor", "rep(0.1,3)", 0.1, NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0, 
        0)
    out
}
asreml.arma <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("arma", "cor", "rep(0.1,2)", numeric(0), 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.armah <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("armah", "cor", "rep(0.1,2)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.armav <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("armav", "cor", "rep(0.1,2)", 0.1, NA, 
        obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 
        0)
    out
}
asreml.asUnivariate <-
function (stack, data, response = "y", traitName = "Trait", traitsWithinUnits = TRUE) 
{
    if (is.character(stack)) {
        stack <- match(stack, names(data))
        if (any(is.na(stack))) 
            stop("names in stack not matched")
    }
    if (max(stack) > ncol(data)) 
        stop("Column in stack out of range")
    if (length(colClass <- unique(sapply(data[, stack, drop = FALSE], 
        data.class))) != 1) 
        stop("All columns to be stacked must be of the same type.")
    if (colClass != "numeric") 
        stop("Stack columns must be numeric")
    trait <- names(data)[stack]
    noStack <- seq(1, ncol(data))[-stack]
    nRows <- sapply(data[, stack], length)
    if (traitsWithinUnits) {
        if (length(unique(nRows)) != 1) 
            stop("Unequal length columns")
        n <- nRows[1]
        vec <- data.frame(V1 = as.vector(t(as.matrix(data[, stack]))), 
            V2 = rep(trait, n))
        rep.vec <- rep(1:n, rep(length(nRows), n))
    }
    else {
        vec <- as.data.frame(cbind(as.vector(as.matrix(data[, 
            stack])), rep(trait, nRows)))
        rep.vec <- unlist(lapply(nRows, function(x) {
            seq(1, x)
        }))
    }
    repData <- as.data.frame(cbind(vec, data[rep.vec, noStack, 
        drop = FALSE]))
    names(repData) <- c(response, traitName, names(data)[noStack])
    return(repData)
}
asreml.at <-
function (obj, lvls, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to at() must be a simple object (factor)\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "at"
    out[[2]] <- asreml.spCall(sys.call())
    obj <- as.character(substitute(obj))
    out[[3]] <- obj
    data <- data[[obj]]
    if (!is.factor(data)) 
        stop(paste("Argument", obj, "to at() must be a factor\n"))
    if (missing(lvls)) 
        at.lev <- levels(data)
    else {
        at.lev <- eval(lvls)
    }
    if (is.numeric(at.lev)) 
        at.lev <- levels(data)[at.lev]
    if (any(is.na(match(at.lev, levels(data))))) 
        stop(paste("Not all levels specified in at() match those in", 
            obj, "\n"))
    out[[4]] <- at.lev
    out[[5]] <- NA
    out[[6]] <- ""
    out[[7]] <- ""
    out[[8]] <- 2
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(-15, 0, NA, NA)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- obj
    out[[13]] <- NA
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.binomial <-
function (link = "logit", dispersion = 1) 
{
    link <- as.character(substitute(link))
    misnames <- c("logit", "probit", "cloglog", "Logit", "Probit", 
        "clog-log", "Cloglog", "Clog-log")
    corresp <- c(1, 2, 3, 1, 2, 3, 3, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Binomial links are \"logit\", \"probit\", or \"cloglog\"\n")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("binomial", link = link)
    fam$dispersion <- dispersion
    fam
}
asreml.call <-
function (dataFrame, asr.inter, asr.struc, asr.glm, predict, 
    asr.predict, weight, offset, multivariate, ntrt, control, 
    keep, constraints, dpmv, debug) 
{
    nRow <- attr(dataFrame, "nRow")
    mbfList <- attr(dataFrame, "MBF")
    for (i in names(control)) {
        if (i %in% c("maxiter", "knots", "nsppoints", "Ftest", 
            "aom", "trace")) 
            storage.mode(control[[i]]) <- "integer"
        else if (i %in% c("workspace", "pworkspace", "splstp", 
            "splinescale", "stepsize")) 
            storage.mode(control[[i]]) <- "double"
    }
    for (i in names(asr.inter)) {
        if (i %in% c("inter", "flev", "nlev", "IBV", "JBV", "neword", 
            "submds", "submod", "neq", "neqd")) 
            storage.mode(asr.inter[[i]]) <- "integer"
        else if (i %in% c("points", "faconst")) 
            storage.mode(asr.inter[[i]]) <- "double"
    }
    gkrige = 1
    asr.predict$lpvals <- length(asr.predict$pvals)
    asr.predict$pvals <- c(asr.predict$pvals, rep(0, sum(asr.inter$flev)))
    asr.predict$lnpvec <- length(asr.predict$pvals)
    pvrow <- asr.predict$lenPrdvals
    pvcol <- 3
    if (asr.glm$link > 1) 
        pvcol <- 5
    asr.predict$prdvals <- matrix(0, nrow = pvrow, ncol = pvcol)
    asr.predict$vpvmat <- rep(0, max(asr.predict$lenVpvals, 1))
    asr.predict$avsed <- rep(0, 10)
    names(asr.predict$avsed) <- c("overall", "min", "mean", "max", 
        "", "", "", "", "", "")
    asr.predict$pvtext <- raw(32768)
    for (i in names(asr.predict)) {
        if (i %in% c("npred", "nfinmd", "lpvals", "lpwts", "nptabl", 
            "npfact", "ivals", "lnpvec")) 
            storage.mode(asr.predict[[i]]) <- "integer"
        else if (i %in% c("pvals", "prdvals", "vpvmat", "avsed")) 
            storage.mode(asr.predict[[i]]) <- "double"
    }
    ddf <- control$denDF
    typeiii <- control$ssType - 1
    if ((asr.inter$neq < 10000) & (asr.inter$neqd < 1000)) {
        if (ddf == 0) 
            ddf <- 6
    }
    else if (ddf >= 0) 
        ddf <- 1
    aodev <- 0
    ny <- 1
    nycol <- 1
    if (is.null(weight)) 
        wtcol <- 0
    else wtcol <- match(weight, names(dataFrame))
    if (is.null(offset)) 
        offsetcol <- 0
    else offsetcol <- match(offset, names(dataFrame))
    storage.mode(asr.glm$id) <- "integer"
    storage.mode(asr.glm$link) <- "integer"
    asr.glm$glmphi <- ifelse(is.null(asr.glm$phi), 1, asr.glm$phi)
    storage.mode(asr.glm$glmphi) <- "double"
    storage.mode(asr.glm$thrlevels) <- "integer"
    ngamma <- length(asr.struc$gammas)
    gammas <- asr.struc$gammas
    asr.struc$printGammas <- rep(1, ngamma)
    asr.struc$printGammas[asreml.grep("<NotEstimated>$", names(gammas))] <- 0
    asr.struc$gammas <- c(gammas, 0)
    for (i in names(asr.struc)) {
        if (i %in% c("struc", "nspat", "nrsect", "tgamma", "printGammas", 
            "tyfact", "gmcstr", "nfrstg")) 
            storage.mode(asr.struc[[i]]) <- "integer"
        else if (i %in% c("gammas")) 
            storage.mode(asr.struc[[i]]) <- "double"
    }
    constr <- asreml.gmcon(gammas, constraints)
    storage.mode(constr$gmcon) <- "double"
    storage.mode(constr$gequal) <- "integer"
    storage.mode(constr$ngmcon) <- "integer"
    if (debug) {
        cat("NGMCON ", constr$ngmcon, "\n")
        print(cbind(seq(1, length(gammas)), gammas, constr$gequal, 
            constr$gmcon, asr.struc$gmcstr))
    }
    spline <- c(-3, -5, -7, -13, -24)
    mxSpline <- max(1, sum(as.numeric(is.element(asr.inter$inter[1, 
        ], spline))))
    naiopt <- 60
    aiopt <- rep(0, naiopt)
    win95 <- 0
    if (version$os == "Microsoft Windows") 
        win95 <- 1
    ifault <- 0
    nan <- 0
    if (length(control$ginverse) > 0) 
        nan <- max(sapply(control$ginverse$ainverse, function(x) length(attr(x, 
            "rowNames"))))
    iarg <- c(ncol(dataFrame), nrow(dataFrame), ny, ntrt, nycol, 
        wtcol, offsetcol, ny, ddf, typeiii, aodev, gkrige, mxSpline, 
        nan, win95, as.numeric(debug))
    storage.mode(iarg) <- "integer"
    minus12 <- (asr.inter$inter[1, asr.inter$neword] != -12)
    whichYssq <- ((asr.inter$submds <= asr.inter$submod) & minus12)
    nfactdsub <- sum(asr.inter$fixedFlag[whichYssq])
    nfactd <- sum(asr.inter$fixedFlag)
    mywv <- 1
    if (ddf > 0) {
        aiopt[42] <- length(asr.struc$ainv$ainv) + 1
        lssp <- (asr.inter$neqd * (asr.inter$neqd + 1))/2
        ii <- ngamma - nfactd + 2
        asr.struc$ainv$ainv <- c(asr.struc$ainv$ainv, rep(0, 
            ii * lssp + lssp + 2 * nfactd))
        mywv <- sum(as.numeric((asr.struc$tyfact >= 1) & (asr.struc$tyfact <= 
            9999))) + (ngamma - asr.struc$nrsect[2, 1] + 1)
    }
    if (length(asr.struc$ainv$llainv) == 0) 
        asr.struc$ainv$llainv <- c(0, 0)
    aiopt[45] <- length(asr.struc$ainv$llainv) + 1
    asr.struc$ainv$llainv <- c(asr.struc$ainv$llainv, rep(-9, 
        nfactd * nfactd + 2))
    asr.struc$ainv$llainv[aiopt[45]] <- -113
    storage.mode(asr.struc$ainv$nainv) <- "integer"
    storage.mode(asr.struc$ainv$lcainv) <- "integer"
    storage.mode(asr.struc$ainv$llainv) <- "integer"
    storage.mode(asr.struc$ainv$ainv) <- "double"
    if (length(mbfList) > 0) {
        mbfList$mbfAddr <- rbind(sapply(mbfList$mbfAttr, function(x) {
            x$nux
        }), sapply(mbfList$mbfAttr, function(x) {
            x$ncz
        }))
    }
    else mbfList$mbfAddr <- 0
    storage.mode(mbfList$mbfAddr) <- "integer"
    storage.mode(mbfList$mbfX) <- "double"
    asr <- list()
    nfactp <- length(asr.inter$neword)
    ndense <- asr.inter$neqd - ny
    if (control$Cfixed) {
        asr$Cfixed <- vector(mode = "double", length = ndense * 
            (ndense + 1)/2)
        asr$Cfixed[1] <- 1
    }
    else asr$Cfixed <- -1
    asr$CfacNum <- 0
    asr$nCfac <- 0
    if (length(control$Csparse[[2]]) > 0) {
        Clen <- sum(asr.inter$Coptions$Cflev * (asr.inter$Coptions$Cflev + 
            1)/2)
        asr$Csparse <- vector(mode = "double", length = Clen * 
            4)
        asr$CfacNum <- match(asr.inter$Coptions$Cfacnam, asr.inter$facnam)
        asr$nCfac <- length(asr$CfacNum)
    }
    else asr$Csparse <- -1
    if (control$aom) {
        aiopt[50] <- 4
        aiopt[51] <- 3
    }
    asr$aovTbl <- rep(0, nfactdsub * 7)
    asr$apstvar <- rep(-1e+37, mywv * (mywv + 3))
    storage.mode(asr$aovTbl) <- "double"
    storage.mode(asr$apstvar) <- "double"
    asr$loglik <- vector(mode = "double", length = 1)
    asr$soln <- matrix(0, nrow = asr.inter$neq, ncol = 2 * as.numeric(control$aom) + 
        1)
    storage.mode(asr$soln) <- "double"
    asr$vsoln <- as.double(rep(0, asr.inter$neq))
    asr$resid <- matrix(0, nrow = nRow, ncol = 3 * as.numeric(control$aom) + 
        1)
    storage.mode(asr$resid) <- "double"
    asr$linearPredictors <- as.double(rep(0, nRow))
    asr$hat <- as.double(rep(0, nRow))
    asr$sigma2 <- vector(mode = "double", length = 1)
    asr$deviance <- vector(mode = "double", length = 1)
    asr$nedf <- vector(mode = "integer", length = 1)
    asr$ai <- vector(mode = "double", length = ngamma * ngamma)
    asr$nwv <- vector(mode = "integer", length = 1)
    asr$noeff <- vector(mode = "integer", length = nfactp)
    asr$nsing <- vector(mode = "integer", length = 1)
    asr$yssqu <- vector(mode = "double", length = nfactp)
    storage.mode(asr$Cfixed) <- "double"
    storage.mode(asr$Csparse) <- "double"
    storage.mode(asr$CfacNum) <- "integer"
    storage.mode(asr$nCfac) <- "integer"
    asr$prvgam <- as.double(rep(0, ngamma))
    asr$score <- as.double(rep(0, ngamma))
    asr$trArray <- matrix(NA, nrow = ngamma + 3, ncol = control$maxiter + 
        3)
    storage.mode(asr$trArray) <- "double"
    if (asr.struc$ainv$locwts > 0 && !multivariate) 
        aiopt[8] <- asr.struc$ainv$lcainv[3, 1]
    aiopt[19] <- 99
    aiopt[2] <- control$eqorder
    aiopt[38] <- as.integer(0.75 * control$maxiter)
    aiopt[39] <- 5
    aiopt[40] <- control$tol[1]
    aiopt[41] <- control$tol[2]
    aiopt[57] <- control$threads
    aiopt[58] <- 0
    storage.mode(aiopt) <- "integer"
    if (length(control$ginverse) > 0) {
        grpFac <- asr.inter$facnam[asr.inter$nlev > 1]
        for (i in names(control$ginverse$factor)) {
            if ((length(grpFac) > 0 && is.na(match(i, grpFac))) | 
                (length(grpFac) == 0)) {
                dataFrame[[i]] <- factor(as.character(dataFrame[[i]]), 
                  levels = control$ginverse$rownames(i, control$ginverse))
            }
        }
    }
    asr.vec <- as.vector(t(sapply(dataFrame, function(x, dpmv) {
        if (is.factor(x)) {
            which <- (is.na(x) | (as.character(x) == "NA"))
            x <- as.numeric(x)
        }
        else which <- is.na(x)
        x[which] <- dpmv
        x
    }, dpmv)))
    storage.mode(asr.vec) <- "double"
    asr$errtxt <- raw(80)
    asr$converge <- 0
    storage.mode(asr$converge) <- "integer"
    asr$where <- asreml.findMe()
    asr$license <- raw(1024)
    storage.mode(dpmv) <- "double"
    if (debug) {
        print(cbind(names(asr.predict), sapply(asr.predict, function(x) storage.mode(x))))
        print(cbind(names(control), sapply(control, function(x) storage.mode(x))))
        print(cbind(names(asr), sapply(asr, function(x) storage.mode(x))))
    }
    ifault <- .Call("casr", asr.vec, iarg, dpmv, asr.inter, asr.struc, 
        asr.struc$ainv, asr.glm, asr.predict, constr, mbfList, 
        control, aiopt, asr)
    ltxt <- (seq(along = asr$license)[asr$license == as.raw(0)])[1] - 
        1
    license <- rawToChar(asr$license[1:ltxt])
    ltxt <- (seq(along = asr$errtxt)[asr$errtxt == as.raw(0)])[1] - 
        1
    if (ltxt == 0) 
        errchar <- character(0)
    else errchar <- rawToChar(asr$errtxt[1:ltxt])
    if (ifault != 0 | aiopt[36] != 0) {
        xtra <- {
            if (ifault == -2) 
                ""
            else "Results may be erroneous"
        }
        if (asreml.Rsys) 
            warning(paste("Abnormal termination\n", errchar, 
                "\n", xtra, sep = ""), call. = FALSE)
        else {
            warning(paste("Abnormal termination\n", errchar, 
                "\n", xtra, sep = ""))
            stop()
        }
    }
    else if (control$trace || !asr$converge) 
        cat(errchar, "\n")
    asr.struc$gammas <- asr.struc$gammas[-length(asr.struc$gammas)]
    gammas <- asr.struc$gammas[asr.struc$printGammas == 1]
    prvgam <- asr$prvgam[asr.struc$printGammas == 1]
    if (asr$converge) {
        asr$prvgam[abs(asr$prvgam) < .Machine$double.eps] <- NA
        pc <- abs((gammas - prvgam)/prvgam)
        if (any(pc[!is.na(pc)] > 0.01)) {
            which <- seq(1, length(pc))[!is.na(pc)]
            which <- which[pc[which] > 0.01]
            pc <- matrix(pc[which], ncol = 1)
            dimnames(pc) <- list(names(gammas)[which], "Change(%)")
            warning("At least one parameter changed by more than 1% on the last iteration\n")
            print(round(pc * 100, 2))
        }
    }
    dimnames(asr$soln) <- list(asreml.labsoln(ny, asr.inter$equations, 
        asr.inter$neq, asr.inter$facnam, asr.inter$varLevels), 
        c("bu", "ginvU", "stdCondBlup")[seq(1, 2 * as.numeric(control$aom) + 
            1)])
    dimnames(asr$resid) <- list(NULL, c("e", "rinvE", "stdCondRes", 
        "")[seq(1, 3 * as.numeric(control$aom) + 1)])
    if (!is.na(match("mv", names(dataFrame)))) 
        asr$resid[!is.na(dataFrame$mv), ] <- NA
    resid <- linearPredictors <- hat <- rep(NA, length(keep))
    resid[keep] <- asr$resid[seq(1, nrow(dataFrame)), 1]
    linearPredictors[keep] <- asr$linearPredictors[seq(1, nrow(dataFrame))]
    hat[keep] <- asr$hat[seq(1, nrow(dataFrame))]
    if (control$aom) {
        asr$resid[, 3][asr$resid[, 3] <= 0] <- NA
        asr$resid[, 3] <- asr$resid[, 2]/sqrt(asr$resid[, 3])
        asr$soln[, 3][asr$soln[, 3] <= 0] <- NA
        asr$soln[, 3] <- asr$soln[, 2]/sqrt(asr$soln[, 3])
        aom <- list(R = asr$resid[, c(2, 3)], G = asr$soln[asr.inter$eqr, 
            c(2, 3)])
    }
    else aom <- list()
    if (debug) 
        kol <- ncol(asr$trArray)
    else kol <- max(1, sum(as.numeric(apply(asr$trArray, 2, function(x) {
        !all(is.na(x))
    }))))
    monitor <- as.data.frame(asr$trArray[seq(1, sum(asr.struc$printGammas) + 
        3), 1:kol])
    if (!debug) {
        row.names(monitor) <- c("loglik", "S2", "df", names(gammas))
        names(monitor) <- as.character(seq(1, kol))
        monitor$constraint <- c(NA, NA, NA, asreml.guzpfx(asr.struc$gmcstr[asr.struc$printGammas == 
            1]))
    }
    coefficients <- vector(mode = "list", length = 3)
    if (length(asr.inter$eqf) > 0) {
        coefficients[[1]] <- asr$soln[asr.inter$eqf, 1]
        names(coefficients[[1]]) <- dimnames(asr$soln)[[1]][asr.inter$eqf]
    }
    if (length(asr.inter$eqr) > 0) {
        coefficients[[2]] <- asr$soln[asr.inter$eqr, 1]
        names(coefficients[[2]]) <- dimnames(asr$soln)[[1]][asr.inter$eqr]
    }
    if (length(asr.inter$eqs) > 0) {
        coefficients[[3]] <- asr$soln[asr.inter$eqs, 1]
        names(coefficients[[3]]) <- dimnames(asr$soln)[[1]][asr.inter$eqs]
    }
    names(coefficients) <- c("fixed", "random", "sparse")
    attr(coefficients, "Terms") <- c(asr.inter$facnam[seq(1, 
        ny)], rep(asr.inter$facnam[asr.inter$equations], asr.inter$flev[asr.inter$equations]))[c(asr.inter$eqf, 
        asr.inter$eqr, asr.inter$eqs)]
    vcoeff <- vector(mode = "list", length = 3)
    vcoeff[[1]] <- asr$vsoln[asr.inter$eqf]
    vcoeff[[2]] <- asr$vsoln[asr.inter$eqr]
    vcoeff[[3]] <- asr$vsoln[asr.inter$eqs]
    names(vcoeff) <- c("fixed", "random", "sparse")
    if (control$Cfixed) {
        eqf1 <- asr.inter$eqf - 1
        xy <- expand.grid(eqf1, eqf1)
        xy <- xy[xy$Var1 <= xy$Var2, ]
        whichCfixed <- (xy$Var2 * (xy$Var2 - 1)/2) + xy$Var1
        Cfixed <- asr$Cfixed[whichCfixed]
    }
    else Cfixed <- NULL
    if (length(control$Csparse[[2]]) > 0) {
        asr$Csparse <- as.data.frame(matrix(asr$Csparse, ncol = 4, 
            byrow = TRUE))
        names(asr$Csparse) <- c("term", "row(i)", "column(j)", 
            "Cij")
        asr$Csparse <- asr$Csparse[asr$Csparse$Cij != 0, ]
        asr$Csparse$term <- asr.inter$facnam[asr.inter$neword[asr$Csparse[, 
            1]]]
    }
    else asr$Csparse <- NULL
    nsect <- ncol(asr.struc$nrsect)
    whichScore <- seq(1, sum(asr.struc$printGammas))
    if (nsect > 1) 
        whichScore <- seq(2, sum(asr.struc$printGammas) + 1)
    yssqu <- asr$yssqu[whichYssq]
    noeff <- asr$noeff[whichYssq]
    nolev <- asr.inter$flev[asr.inter$neword][whichYssq]
    names(yssqu) <- asr.inter$facnam[asr.inter$neword][whichYssq]
    names(noeff) <- asr.inter$facnam[asr.inter$neword][whichYssq]
    names(nolev) <- asr.inter$facnam[asr.inter$neword][whichYssq]
    yssqu <- yssqu[asr.inter$fixedFlag[whichYssq]]
    nzdf <- (noeff[names(yssqu)] > 0)
    aovTbl <- matrix(asr$aovTbl, ncol = 7, byrow = TRUE)
    if (length(nzdf) > 0) {
        aovTbl <- aovTbl[seq(1, sum(as.numeric(nzdf))), , drop = FALSE]
        dimnames(aovTbl) <- list(rev(names(yssqu)[nzdf]), c("id", 
            "df", "denDF", "Finc", "Fcon", "M", "Fprob"))
    }
    else aovTbl <- NULL
    ai <- asr$ai[seq(1, asr$nwv * (asr$nwv + 1)/2)]
    which <- asr.struc$printGammas[seq(nfactp + 1, ngamma)]
    newAI <- matrix(as.logical(which %o% which), nrow = length(which))
    newAI <- ai[newAI[!lower.tri(newAI)]]
    if (all(asr$apstvar < -1e+36)) 
        sv <- NULL
    else {
        asr$apstvar <- round(asr$apstvar, 10)
        awv <- (sqrt(9 + 4 * sum(as.numeric(asr$apstvar > -1e+37))) - 
            3)/2
        asr$apstvar <- asr$apstvar[asr$apstvar > -1e+37]
        sv <- matrix(asr$apstvar[1:(awv * (awv + 3))], nrow = awv, 
            byrow = FALSE)
        fn <- sv[, 1]
        fn[asr.struc$tyfact[fn] == 0] <- 0
        nn <- asr.inter$facnam[asr.inter$neword[fn[fn > 0]]]
        svnam <- c(nn[rank(match(nn, names(asr.struc$gammas)))], 
            names(asr.struc$gammas)[asr.struc$nrsect[2, 1]:ngamma])[1:awv]
        sv <- matrix(sv[, 2:(awv + 3)], nrow = nrow(sv))
        dimnames(sv) <- list(svnam, c("df", "Variance", svnam))
    }
    if (asr.predict$npred > 0) {
        ltxt <- (seq(along = asr.predict$pvtext)[asr.predict$pvtext == 
            as.raw(0)])[1] - 1
        pvtext <- rawToChar(asr.predict$pvtext[1:ltxt])
        asr.predict$vpvmat[asr.predict$vpvmat > 9.9e+36] <- NA
    }
    predictions <- asreml.prdList(predict, asr.predict, asr.inter, 
        asr.predict$prdvals, asr.predict$vpvmat, asr.predict$avsed, 
        pvtext)
    gammas.con = asr.struc$gmcstr[asr.struc$printGammas == 1]
    names(gammas.con) <- asreml.guzpfx(gammas.con)
    if (debug) 
        return(list(monitor = monitor, loglik = asr$loglik, gammas = gammas, 
            gammas.type = asr.struc$tgamma[asr.struc$printGammas == 
                1], gammas.con = gammas.con, stratumVariances = sv, 
            score = asr$score[whichScore], coefficients = coefficients, 
            vcoeff = vcoeff, predictions = predictions, residuals = resid, 
            linear.predictors = linearPredictors - resid, hat = hat, 
            sigma2 = asr$sigma2, deviance = -asr$deviance, nedf = asr$nedf, 
            ai = newAI, nwv = asr$nwv, noy = ntrt, nolev = nolev, 
            noeff = noeff, nsing = asr$nsing, yssqu = yssqu, 
            Cfixed = Cfixed, Csparse = asr$Csparse, aovTbl = aovTbl, 
            asr.vec = matrix(asr.vec, nrow = nrow(dataFrame), 
                byrow = TRUE), aom = aom))
    list(monitor = monitor, loglik = asr$loglik, gammas = gammas, 
        gammas.type = asr.struc$tgamma[asr.struc$printGammas == 
            1], gammas.con = gammas.con, stratumVariances = sv, 
        score = asr$score[whichScore], coefficients = coefficients, 
        vcoeff = vcoeff, predictions = predictions, residuals = resid, 
        linear.predictors = linearPredictors - resid, hat = hat, 
        aom = aom, sigma2 = asr$sigma2, deviance = -asr$deviance, 
        nedf = asr$nedf, ai = newAI, nwv = asr$nwv, noy = ntrt, 
        nolev = nolev, noeff = noeff, nsing = asr$nsing, yssqu = yssqu, 
        Cfixed = Cfixed, Csparse = asr$Csparse, aovTbl = aovTbl, 
        converge = as.logical(asr$converge), last.message = errchar, 
        license = license, ifault = ifault)
}
asreml.CallEval <-
function (expr, data, Rcov, AND) 
{
    if (is.character(expr)) 
        expr <- parse(text = expr)
    switch(mode(expr[[1]]), call = {
        if (asreml.Rsys) expr11 <- as.character(expr[[1]][[1]]) else expr11 <- expr[[1]][[1]]
        if (!is.na(match(expr11, asreml.Spcls$Fun))) {
            if (!is.null(data)) {
                ll <- length(expr[[1]]) + 1
                nn <- names(expr[[1]])
                if (is.null(nn)) nn <- rep("", (ll - 1))
                expr[[1]][[ll]] <- as.name(data)
                names(expr[[1]]) <- c(nn, "data")
            }
            ll <- length(expr[[1]]) + 1
            nn <- names(expr[[1]])
            if (is.null(nn)) nn <- rep("", (ll - 1))
            expr[[1]][[ll]] <- Rcov
            names(expr[[1]]) <- c(nn, "Rcov")
            expr[[1]][[1]] <- as.name(paste("asreml.", expr[[1]][[1]], 
                sep = ""))
        } else expr[1] <- expr[1]
        if (length(expr[[1]]) > 1) {
            for (i in 2:length(expr[[1]])) {
                expr[[1]][i] <- asreml.CallEval(expr[[1]][i], 
                  data, Rcov)
            }
        }
    }, numeric = {
        expr[1] <- expr[1]
    }, name = {
        expr[1] <- expr[1]
    }, character = {
        expr[1] <- expr[1]
    }, logical = {
        expr[1] <- expr[1]
    }, stop(paste("Dunno what to do with", expr, "\n")))
    expr
}
asreml.chol <-
function (obj, k = 1, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("chol", "var", numeric(0), "matrix(0.1,nrow=n,ncol=n)+diag(0.05,nrow=n)", 
        k, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 11, 0, 0, k + 1, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0)
    out
}
asreml.cholc <-
function (obj, k = 1, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("cholc", "var", numeric(0), "matrix(0.1,nrow=n,ncol=n)+diag(0.05,nrow=n)", 
        k, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 11, 0, 0, -(k + 1), 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0)
    out
}
asreml.cir <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "cir"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ii"
    struc <- c(0, 0, 6, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.cirh <-
function (x, y, init = NA, data, ...) 
{
    fun <- "cirh"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ih"
    struc <- c(0, 0, 6, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.cirv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "cirv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "iv"
    struc <- c(0, 0, 6, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.con <-
function (obj, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to con() must be a factor\n")
    out <- asreml.spc("con", "mdl", 0, 0, NA, substitute(obj), 
        NA, data, match.call()$Rcov)
    out[[4]] <- out[[4]][-length(out[[4]])]
    out[[10]] <- c(length(out[[4]]), 0, NA, -1)
    out
}
asreml.constraints <-
function (form, gammas, drop.unused.levels = TRUE, intercept = FALSE, 
    na.action = na.include) 
{
    if (missing(form)) 
        return(NULL)
    if (missing(gammas)) 
        stop("\nmissing gammas argument")
    if (is.null(gammas$Gamma)) 
        stop("\ngammas argument must be a data frame with a component named 'Gamma'\n")
    if (!inherits(form, "formula")) 
        stop("\nArgument must be a formula")
    tt <- terms(form)
    if (attr(tt, "intercept") == 1 && !intercept) 
        form <- asreml.formula(parse(text = paste("~", paste(c("-1", 
            as.character(form[2])), collapse = "+"), sep = " ")))
    Mf <- model.frame.default(form, data = gammas, drop.unused.levels = drop.unused.levels, 
        na.action = na.action)
    Mm <- model.matrix(form, data = Mf)
    dimnames(Mm) <- list(as.character(gammas$Gamma), dimnames(Mm)[[2]])
    return(Mm)
}
asreml.control <-
function (knots = 50, nsppoints = 21, splinepoints = list(), 
    predictpoints = list(), grid = TRUE, splinescale = -1, spline.step = list(dev = 10000, 
        pol = 10000), pwrpoints = list(), ginverse = list(), 
    mbf = list(), group = list(), denDF = -1, ssType = 1, Ftest = formula("~NULL"), 
    drop.unused.levels = FALSE, Cfixed = FALSE, Csparse = formula("~NULL"), 
    aom = FALSE, equate.levels = character(0), trace = TRUE, 
    maxiter = 13, stepsize = 0.1, tol = c(0, 0), emflag = TRUE, 
    workspace = 1.6e+07, pworkspace = 1.6e+07, eqorder = 3, threads = 1, 
    ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!missing(Ftest) && length(Ftest[[2]]) == 0) 
        call$Ftest <- Ftest
    if (!missing(ginverse)) 
        call$ginverse <- substitute(ginverse)
    splstp <- c(0, 0, 10000, 0, 10000, 0, 10000, 0, 0, 10000)
    names(splstp) <- c("", "", "spl", "", "dev", "", "pol", "", 
        "", "")
    if (any(is.na(j <- match(names(spline.step), names(splstp))))) 
        stop("Unrecognzed component of spline.step\n")
    splstp[sort(j)] <- unlist(spline.step)[match(names(splstp), 
        names(spline.step), nomatch = 0)]
    if (length(predictpoints) > 0) {
        if (is.list(predictpoints)) {
            if (length(grid) == 1) 
                grid <- rep(grid, length(predictpoints))
            for (i in 1:length(grid)) {
                if (is.list(predictpoints[[i]])) {
                  if (length(predictpoints[[i]]) > 1) {
                    if (grid[i]) {
                      grd <- expand.grid(predictpoints[[i]][[2]], 
                        predictpoints[[i]][[1]])
                      predictpoints[[i]][[1]] <- as.numeric(as.character(grd[, 
                        2]))
                      predictpoints[[i]][[2]] <- as.numeric(as.character(grd[, 
                        1]))
                    }
                    else {
                      if (length(predictpoints[[i]][[1]]) != 
                        length(predictpoints[[i]][[2]])) 
                        stop("Predict points components 'x' and 'y' must be equal length\n")
                    }
                  }
                }
                else {
                  temp <- predictpoints[[i]]
                  predictpoints[[i]] <- vector(mode = "list")
                  predictpoints[[i]]$x <- temp
                }
            }
        }
        else {
            temp <- predictpoints[[i]]
            predictpoints[[i]] <- vector(mode = "list")
            predictpoints[[i]]$x <- temp
        }
    }
    if (length(ginverse) > 0) {
        giv <- vector(mode = "list", length = 9)
        names(giv) <- c("factor", "ainverse", "gname", "rownames", 
            "groups", "groupoffset", "row", "col", "ai")
        giv$factor <- rep(0, length(ginverse))
        names(giv$factor) <- names(ginverse)
        cc <- as.call(parse(text = deparse(call$ginverse)))
        tgt <- unlist(lapply(cc[[1]][-1], function(x) paste(as.character(x), 
            collapse = "|")))
        names(tgt) <- names(ginverse)
        giv$ainverse <- vector(mode = "list", length = length(unique(tgt)))
        names(giv$ainverse) <- unique(tgt)
        giv$gname <- function(fac, where) {
            names(where$ainverse)[where$factor[fac]]
        }
        giv$rownames <- function(fac, where) {
            attr(where$ainverse[[where$factor[fac]]], "rowNames")
        }
        giv$groups <- function(fac, where) {
            attr(where$ainverse[[where$factor[fac]]], "geneticGroups")[1]
        }
        giv$groupoffset <- function(fac, where) {
            attr(where$ainverse[[where$factor[fac]]], "geneticGroups")[2]
        }
        giv$row <- function(fac, where) {
            (where$ainverse[[where$factor[fac]]])[, 1]
        }
        giv$col <- function(fac, where) {
            (where$ainverse[[where$factor[fac]]])[, 2]
        }
        giv$ai <- function(fac, where) {
            (where$ainverse[[where$factor[fac]]])[, 3]
        }
        for (i in 1:length(ginverse)) {
            who <- names(ginverse)[i]
            if (is.data.frame(ginverse[[i]])) {
                if (ncol(ginverse[[i]]) < 3) 
                  stop("Less than 3 columns in ginverse data frame\n")
                if (ncol(ginverse[[i]]) > 3) 
                  warning("Extra columns in ginverse ignored\n")
                which <- match(c("row", "column"), casefold(names(ginverse[[i]]))[1:3])
                if (any(is.na(which))) 
                  stop("Cannot match 'row' or 'column' in ginverse\n")
                for (w in which) {
                  if (is.factor(ginverse[[i]][, w])) {
                    ginverse[[i]][, w] <- as.numeric(ginverse[[i]][, 
                      w])
                    warning("Factor in ginverse coerced to numeric")
                  }
                }
                which <- c(which, (1:3)[-which])
                if (any(diff(order(ginverse[[i]][, which[1]], 
                  ginverse[[i]][, which[2]])) < 0)) 
                  stop("Ginverse must be sorted as columns within rows\n")
                if (is.null(attr(ginverse[[i]], "rowNames"))) {
                  warning("Missing rowNames attribute for ginverse\n")
                  rowNames <- seq(1, length(unique(ginverse[[i]][, 
                    which[1]])))
                }
                else rowNames <- attr(ginverse[[i]], "rowNames")
                if (is.null(attr(ginverse[[i]], "geneticGroups"))) 
                  grp <- 0
                else grp <- attr(ginverse[[i]], "geneticGroups")
                ginverse[[i]] <- as.matrix(ginverse[[i]][, which])
                attr(ginverse[[i]], "rowNames") <- rowNames
                attr(ginverse[[i]], "geneticGroups") <- grp
            }
            else if (is.matrix(ginverse[[i]])) {
                if (nrow(ginverse[[i]]) == ncol(ginverse[[i]])) {
                  if (ncol(ginverse[[i]]) <= 3) 
                    stop("Convert ginverse matrix to a data frame\n")
                  if (length(dimnames(ginverse[[i]])[[1]]) == 
                    0) {
                    warning("Missing row component of dimnames for ginverse matrix\n")
                    rowNames <- seq(1, nrow(ginverse[[i]]))
                  }
                  else rowNames <- dimnames(ginverse[[i]])[[1]]
                  if (is.null(attr(ginverse[[i]], "geneticGroups"))) {
                    grp <- 0
                  }
                  else grp <- attr(ginverse[[i]], "geneticGroups")
                  mat <- matrix((!is.na(as.vector(ginverse[[i]]))), 
                    nrow = nrow(ginverse[[i]]), ncol = ncol(ginverse[[i]]))
                  flag <- as.vector(t(lower.tri(ginverse[[i]], 
                    diag = TRUE)))
                  flag <- flag & t(mat)
                  row <- as.vector(col(mat))[flag]
                  col <- as.vector(row(mat))[flag]
                  ginverse[[i]] <- cbind(row, col, as.vector(t(ginverse[[i]]))[flag])
                  attr(ginverse[[i]], "rowNames") <- rowNames
                  attr(ginverse[[i]], "geneticGroups") <- grp
                }
                else if (ncol(ginverse[[i]]) == 3) {
                  which <- match(c("row", "column"), dimnames(ginverse[[i]])[[2]])
                  if (any(is.na(which))) 
                    stop("Cannot match 'row' or 'column' in dimnames(ginverse)[[2]]\n")
                  which <- c(which, (1:3)[-which])
                  if (any(diff(order(ginverse[[i]][, which[1]], 
                    ginverse[[i]][, which[2]]))) < 0) 
                    stop("Ginverse must be sorted as columns within rows\n")
                  if (is.null(attr(ginverse[[i]], "rowNames"))) {
                    warning("Missing rowNames attribute for ginverse\n")
                    rowNames <- seq(1, length(unique(ginverse[[i]][, 
                      which[1]])))
                  }
                  else rowNames <- attr(ginverse[[i]], "rowNames")
                  if (is.null(attr(ginverse[[i]], "geneticGroups"))) 
                    grp <- 0
                  else grp <- attr(ginverse[[i]], "geneticGroups")
                  ginverse[[i]] <- ginverse[[i]][, which]
                  attr(ginverse[[i]], "rowNames") <- rowNames
                  attr(ginverse[[i]], "geneticGroups") <- grp
                }
            }
            else if (class(ginverse[[i]]) != "list" && is.vector(ginverse[[i]])) {
                ginverse[[i]][ginverse[[i]] == 0] <- NA
                n <- (sqrt(1 + 8 * length(ginverse[[i]])) - 1)/2
                mat <- matrix(nrow = n, ncol = n)
                flag <- as.vector(t(lower.tri(mat, diag = TRUE)))
                row <- as.vector(col(mat))[flag]
                col <- as.vector(row(mat))[flag]
                flag <- !is.na(ginverse[[i]])
                if (is.null(attr(ginverse[[i]], "rowNames"))) {
                  warning("Missing rowNames attribute for ginverse\n")
                  rowNames <- seq(1, length(unique(row[flag])))
                }
                else rowNames <- attr(ginverse[[i]], "rowNames")
                if (is.null(attr(ginverse[[i]], "geneticGroups"))) {
                  grp <- 0
                }
                else grp <- attr(ginverse[[i]], "geneticGroups")
                ginverse[[i]] <- cbind(row[flag], col[flag], 
                  ginverse[[i]][flag])
                attr(ginverse[[i]], "rowNames") <- rowNames
                attr(ginverse[[i]], "geneticGroups") <- grp
            }
            else stop("G-inverse must be a dataframe, matrix or vector of full lower triangle row X row")
            which <- match(tgt[who], unique(tgt))
            giv$factor[who] <- which
            if (is.null(giv$ainverse[[which]])) 
                giv$ainverse[[which]] <- ginverse[[i]]
        }
    }
    else giv <- NULL
    if (length(mbf) > 0) {
        if (!is.list(mbf)) 
            stop("mbf must be a list\n")
        lapply(mbf, function(x) {
            nn <- names(x)
            if (is.na(match("dataFrame", nn)) | is.na(match("key", 
                nn))) 
                stop("Elements of mbf must have components 'dataFrame' and 'key'\n")
        })
    }
    if (length(group) > 0) {
        if (!is.list(group)) 
            stop("group must be a list object\n")
    }
    if (!is.logical(trace)) 
        stop("\n trace must be one of TRUE or FALSE.\n")
    splXtras <- vector(mode = "list", length = 9)
    names(splXtras) <- c("knots", "nsppoints", "splstp", "splinescale", 
        "IBV", "JBV", "points", "sPtr", "dimTbl")
    if (length(splinepoints) > 0 || length(predictpoints) > 0) 
        splXtras$points <- vector(mode = "double", length = 0)
    if (length(splinepoints) > 0) {
        splXtras$sPtr <- matrix(0, nrow = length(splinepoints), 
            ncol = 2)
        dimnames(splXtras$sPtr) <- list(names(splinepoints), 
            c("loc", "len"))
        for (i in 1:length(splinepoints)) {
            if (length(splinepoints[[i]]) < 3) 
                stop("Insufficient knots for spline (must have at least 3)\n")
            splXtras$sPtr[i, 1] <- length(splXtras$points) + 
                1
            splXtras$sPtr[i, 2] <- length(splinepoints[[i]])
            splXtras$points <- c(splXtras$points, splinepoints[[i]])
        }
    }
    if (length(predictpoints) > 0) {
        splXtras$IBV <- matrix(0, ncol = length(predictpoints), 
            nrow = 4)
        np <- names(predictpoints)
        dimnames(splXtras$IBV) <- list(c("flag", "loc", "data", 
            "len"), np)
        ndim2 <- (sapply(predictpoints, function(x) {
            length(x)
        }) > 1)
        if (sum(as.numeric(ndim2)) > 0) {
            splXtras$JBV <- matrix(0, ncol = sum(as.numeric(ndim2)), 
                nrow = 4)
            dimnames(splXtras$JBV) <- list(c("flag", "loc", "data", 
                "len"), np)
            splXtras$dimTbl <- vector(mode = "list", length = 2)
            names(splXtras$dimTbl) <- c("term", "axis")
            splXtras$dimTbl$term <- vector(mode = "character")
            splXtras$dimTbl$axis <- vector(mode = "character")
        }
        for (i in 1:length(predictpoints)) {
            splXtras$IBV[2, i] <- length(splXtras$points) + 1
            splXtras$IBV[4, i] <- length(predictpoints[[i]][[1]])
            splXtras$IBV[1, i] <- 1
            splXtras$points <- c(splXtras$points, predictpoints[[i]][[1]])
            nax <- names(predictpoints[[i]])
            for (j in 1:length(predictpoints[[i]])) {
                splXtras$dimTbl$term <- c(splXtras$dimTbl$term, 
                  np[i])
                splXtras$dimTbl$axis <- c(splXtras$dimTbl$axis, 
                  nax[j])
            }
            if (ndim2[i]) {
                splXtras$JBV[2, i] <- length(splXtras$points) + 
                  1
                splXtras$JBV[4, i] <- length(predictpoints[[i]][[2]])
                splXtras$JBV[1, i] <- 1
                splXtras$points <- c(splXtras$points, predictpoints[[i]][[2]])
            }
        }
    }
    splXtras$dimTbl <- as.data.frame(splXtras$dimTbl)
    splXtras$knots <- knots
    splXtras$nsppoints <- nsppoints
    splXtras$splstp <- splstp
    splXtras$splinescale <- splinescale
    return(list(splXtras = splXtras, knots = knots, nsppoints = nsppoints, 
        splinepoints = splinepoints, predictpoints = predictpoints, 
        splinescale = splinescale, splstp = splstp, pwrpoints = pwrpoints, 
        ginverse = giv, mbf = mbf, group = group, denDF = denDF, 
        ssType = ssType, Ftest = Ftest, drop.unused.levels = drop.unused.levels, 
        Cfixed = Cfixed, Csparse = Csparse, aom = aom, equate.levels = equate.levels, 
        trace = trace, maxiter = maxiter, stepsize = stepsize, 
        tol = tol, emflag = emflag, workspace = workspace, pworkspace = pworkspace, 
        eqorder = eqorder, threads = threads, call = call))
}
asreml.cor <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("cor", "cor", 0.1, numeric(0), NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.corb <-
function (obj, k = 1, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- obj[[1]]
        out[[2]] <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k))
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k) 
                stop("corb(): Wrong number of initial values\n")
            init <- c(init, rep(0, w - k))
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)))
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k))), NULL))
        out[[8]] <- rep(3, w)
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- obj[[10]]
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
        out[[13]] <- FALSE
    }
    else {
        out[[1]] <- "corb"
        obj <- as.character(substitute(obj))
        out[[2]] <- obj
        out[[3]] <- obj
        drop.unused.levels <- attr(data, "Control")$drop.unused.levels
        data <- data[[obj]]
        if (is.factor(data)) 
            out[[4]] <- asreml.levels(data, as.logical(match.call()$Rcov))
        else stop(paste(obj, "must be a factor\n"))
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k))
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k) 
                stop("corb(): Wrong number of initial values\n")
            init <- c(init, rep(0, w - k))
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)))
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k))), NULL))
        out[[8]] <- rep(3, w)
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- FALSE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.corbh <-
function (obj, k = 1, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- obj[[1]]
        out[[2]] <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k), rep(0.1, n))
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k + n) 
                stop("corbh(): Wrong number of initial values\n")
            tmp <- names(init)
            init <- c(init[1:k], rep(0, w - k), init[(k + 1):(k + 
                n)])
            if (length(tmp) > 0) 
                names(init) <- c(tmp[1:k], rep("", w - k), tmp[(k + 
                  1):(k + n)])
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)), rep("P", 
                n))
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k)), sep = ""), NULL), paste("!", out[[4]], 
                ".var", sep = ""))
        out[[8]] <- c(rep(3, w), rep(2, n))
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- obj[[10]]
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
        out[[13]] <- TRUE
    }
    else {
        out[[1]] <- "corbh"
        obj <- as.character(substitute(obj))
        out[[2]] <- obj
        out[[3]] <- obj
        drop.unused.levels <- attr(data, "Control")$drop.unused.levels
        data <- data[[obj]]
        if (is.factor(data)) 
            out[[4]] <- asreml.levels(data, as.logical(match.call()$Rcov))
        else stop(paste(obj, "must be a factor\n"))
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k), rep(0.1, n))
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k + n) 
                stop("corbh(): Wrong number of initial values\n")
            tmp <- names(init)
            init <- c(init[1:k], rep(0, w - k), init[(k + 1):(k + 
                n)])
            if (length(tmp) > 0) 
                names(init) <- c(tmp[1:k], rep("", w - k), tmp[(k + 
                  1):(k + n)])
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)), rep("P", 
                n))
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k)), sep = ""), NULL), paste("!", out[[4]], 
                ".var", sep = ""))
        out[[8]] <- c(rep(3, w), rep(2, n))
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.corbv <-
function (obj, k = 1, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- obj[[1]]
        out[[2]] <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k), 0.1)
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k + 1) 
                stop("corbv(): Wrong number of initial values\n")
            tmp <- names(init)
            init <- c(init[1:k], rep(0, w - k), init[k + 1])
            if (length(tmp) > 0) 
                names(init) <- c(tmp[1:k], rep("", w - k), tmp[k + 
                  1])
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)), "P")
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k)), sep = ""), NULL), "!var")
        out[[8]] <- c(rep(3, w), 2)
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- obj[[10]]
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
        out[[13]] <- TRUE
    }
    else {
        out[[1]] <- "corbv"
        obj <- as.character(substitute(obj))
        out[[2]] <- obj
        out[[3]] <- obj
        drop.unused.levels <- attr(data, "Control")$drop.unused.levels
        data <- data[[obj]]
        if (is.factor(data)) 
            out[[4]] <- asreml.levels(data, as.logical(match.call()$Rcov))
        else stop(paste(obj, "must be a factor\n"))
        n <- length(out[[4]])
        w <- n - 1
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, k), rep(0, w - k), 0.1)
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != k + 1) 
                stop("corbv(): Wrong number of initial values\n")
            tmp <- names(init)
            init <- c(init[1:k], rep(0, w - k), init[k + 1])
            if (length(tmp) > 0) 
                names(init) <- c(tmp[1:k], rep("", w - k), tmp[k + 
                  1])
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", k), rep("F", (w - k)), "P")
        else out[[6]] <- con
        out[[7]] <- c(paste("!", "cor", seq(1, k), sep = ""), 
            asreml.ie(k < w, paste("!", rep("<NotEstimated>", 
                (w - k)), sep = ""), NULL), "!var")
        out[[8]] <- c(rep(3, w), 2)
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 5, 0, 0, k, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.corg <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("corg", "cor", "matrix(0.1,nrow=n,ncol=n)", 
        numeric(0), NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.corgh <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("corgh", "cor", "matrix(0.1,nrow=n,ncol=n)", 
        "rep(0.1,n)", NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.corgv <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("corgv", "cor", "matrix(0.1,nrow=n,ncol=n)", 
        0.1, NA, obj, init, data, match.call()$Rcov)
    N <- length(out[[4]]) * (length(out[[4]]) - 1)/2
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, N, 0, 0, 
        0)
    out
}
asreml.corh <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("corh", "cor", 0.1, "rep(0.1,n)", NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.corv <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("corv", "cor", 0.1, 0.1, NA, obj, init, 
        data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
        0)
    out
}
asreml.data <-
function (data, model.y, fixed, random, sparse, rcov, glm, asrAttrib, 
    na.method.X, na.method.Y, as.multivariate, weights, offset, 
    ran.order) 
{
    control <- attr(data, "Control")
    mvtrait <- character(0)
    if (length(model.y) > 1) {
        nn <- names(data)
        nn1 <- nn[-1]
        data <- cbind(as.data.frame(eval(data[, 1])), data[, 
            -1])
        names(data) <- c(model.y, nn1)
        for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    }
    else names(data) <- c(model.y, names(data)[-1])
    nRow <- nrow(data)
    data$units <- factor(seq(1, nRow))
    ntrt <- 1
    S2 <- NA
    thrlevels <- 1
    if (length(model.y) > 1 && is.factor(data[, model.y])) {
        if (glm$id == 9) {
            lvls <- unique(as.character(data[, model.y]))
            data[, model.y] <- match(as.character(data[, model.y]), 
                lvls)
            data$trait <- factor(data[, model.y])
            levels(data$trait) <- c(lvls[-length(lvls)], NA)
            ntrt <- thrlevels <- length(levels(data$trait))
            S2 <- asr.glm$dispersion
            nRow <- nrow(data) * (length(lvls) - 1)
        }
        else stop("Response is a factor\n")
    }
    else if (glm$id == 9) 
        stop("Response must be a factor for multinomial models\n")
    if ((ntrt <- length(model.y)) > 1) {
        traits <- model.y
        S2 <- -1
        y <- data.frame(list(y = as.vector(t(data[, model.y])), 
            trait = factor(rep(traits, nrow(data)), levels = traits)))
        fixed <- asreml.formula(parse(text = paste("y ~", as.character(fixed[3]))))
        which <- is.na(match(names(data), traits))
        data <- cbind(y, lapply(data[which], function(x, n) {
            rep(x, rep(n, length(x)))
        }, ntrt))
        if (any(is.na(y)) && (na.method.Y == "omit" | na.method.X == 
            "omit")) 
            stop("Cannot use na.method=omit in multivariate analysis\n")
        nRow <- nrow(data)
        model.y <- "y"
        mvtrait <- "trait"
        for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    }
    if (!is.null(as.multivariate)) {
        mvtrait <- as.multivariate
        ntrt <- length(unique(data[[mvtrait]]))
        S2 <- -1
        nRow <- nrow(data)
        if ((nRow%%ntrt) != 0) 
            stop("Number of observations must be an integer multiple of the number of traits")
        nu <- nRow/ntrt
        if (all(!duplicated(data[[mvtrait]][ntrt]))) 
            data$units <- factor(rep(seq(1, nu), rep(ntrt, nu)))
        else data$units <- factor(rep(seq(1, nu), ntrt))
    }
    MV <- FALSE
    if (na.method.Y == "include") 
        MV <- TRUE
    form <- unique(c(model.y, attr(terms(asreml.ModelFormula(fixed, 
        random, sparse, rcov, weights, offset, attr(data, "GROUP"), 
        attr(data, "MBF")$mbfAttr, ignore = character(0), names(data))), 
        "term.labels"), "units"))
    rn <- row.names(data)
    data <- asreml.modelFrame(form, model.y, data = data, na.method.Y = na.method.Y, 
        na.method.X = na.method.X, control$drop.unused.levels)
    MV <- MV & any(is.na(data[, model.y]))
    if (MV) {
        data$mv <- asreml.mvfact(data[, model.y])
        if (length(sparse[[2]]) == 0) 
            form <- asreml.formula(parse(text = "~ mv"))
        else form <- asreml.formula(parse(text = paste("~", paste(c("mv", 
            as.character(sparse[2])), collapse = "+"), sep = " ")))
        sparse <- form
    }
    for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    if (length(eql <- control$equate.levels) > 1) {
        if (any(is.na(match(eql, names(data))))) 
            stop("Elements of 'equate.levels' not in data")
        if (!all(sapply(data[eql], function(x) is.factor(x)))) 
            stop("Not all elements of 'equate.levels' are factors")
        ulev <- unique(unlist(lapply(data[eql], function(x) levels(x))))
        for (i in eql) data[[i]] <- factor(as.character(data[[i]]), 
            levels = ulev)
    }
    keep <- !is.na(match(rn, row.names(data)))
    if (glm$id == 9) 
        keep <- rep(keep, length(levels(data$trait)))
    data <- asreml.powerFac(random, data, Rcov = 0)
    data <- asreml.powerFac(rcov, data, Rcov = 1)
    random <- asreml.Gorder(random, data, ran.order)
    temp <- asreml.randomStr(random, data)
    random.Inter <- temp$form
    if ((ndum <- length(temp$dummyTerms)) > 0) {
        nad <- names(temp$dummyTerms)
        nn <- names(data)
        data <- cbind(data, matrix(0, nrow = nrow(data), ncol = ndum))
        names(data) <- c(nn, nad)
        for (i in nad) attr(data[[i]], "levels") <- seq(1, temp$dummyTerms[[i]])
    }
    for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    attr(data, "nRow") <- nRow
    attr(data, "S2") <- S2
    attr(data, "thrlevels") <- thrlevels
    attr(data, "keep") <- keep
    attr(data, "ntrt") <- ntrt
    attr(data, "model.y") <- model.y
    attr(data, "sparse") <- sparse
    attr(data, "fixed") <- fixed
    attr(data, "random") <- random
    attr(data, "random.Inter") <- random.Inter
    attr(data, "mvtrait") <- mvtrait
    data
}
asreml.deparse <-
function (x) 
{
    if (asreml.Rsys) 
        paste(deparse(x), collapse = "")
    else deparse(x)
}
asreml.dev <-
function (obj, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to dev() must be a simple variate\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "dev"
    out[[2]] <- asreml.spCall(sys.call())
    obj <- as.character(substitute(obj))
    out[[3]] <- obj
    if (is.factor(data[[obj]])) 
        stop(paste(obj, "is already a factor.\n"))
    splXtras <- attr(data, "Control")$splXtras
    knots <- splXtras$knots
    kl <- 0
    kptr <- 0
    inter <- c(-5, kptr, 0, kl)
    ibv <- inter
    jbv <- inter
    if (!is.null(splXtras$IBV)) {
        which <- match(obj, dimnames(splXtras$IBV)[[2]])
        if (!is.na(which)) {
            ibv <- splXtras$IBV[, which]
            ibv[1] <- 1
            ibv[3] <- 1
        }
    }
    N <- length(unique(data[[obj]])) + knots + max(splXtras$nsppoints, 
        ibv[4])
    lx <- length(data[[obj]])
    one <- 1
    knotpoints <- vector(mode = "double", length = N)
    nux <- vector(mode = "integer", length = 1)
    ncz <- vector(mode = "integer", length = 1)
    fail <- 0
    temp <- .C("splinek", as.double(data[[obj]]), as.integer(one), 
        as.integer(lx), as.integer(inter), as.integer(ibv), as.integer(jbv), 
        as.double(splXtras$splstp), as.integer(knots), as.integer(splXtras$nsppoints), 
        as.double(splXtras$splinescale), as.double(splXtras$points), 
        as.integer(N), ncz = as.integer(ncz), nux = as.integer(nux), 
        knotpoints = as.double(knotpoints), fail = as.integer(fail), 
        NAOK = TRUE)
    if (temp$fail > 0) 
        stop
    out[[4]] <- temp$knotpoints[seq(1, temp$nux)]
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- 0.1
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != 1) 
            stop("dev(): Wrong number of initial values\n")
    }
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) 
        out[[6]] <- "P"
    else out[[6]] <- con
    out[[7]] <- ""
    out[[8]] <- 2
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(-5, 0, match(obj, names(data)), 0)
    out[[11]] <- list(faconst = NULL, nuxPoints = temp$knotpoints[seq(1, 
        temp$nux)], coords = NULL)
    out[[12]] <- obj
    out[[13]] <- ifelse(IniFlag, TRUE, FALSE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.diag <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("diag", "var", numeric(0), "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.Eval <-
function (expr, data = NULL, Rcov = 0) 
{
    z <- as.character(substitute(data))
    expr <- asreml.CallEval(expr, z, Rcov)
    if (length(expr[[1]]) == 1 & mode(expr[[1]]) == "name") 
        expr <- as.character(expr[1])
    expr
}
asreml.exp <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
        x <- substitute(x)
    fun <- "exp"
    type <- "i"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.expandGrid <-
function (flist, plist) 
{
    asreml.expand.grid <- function(...) {
        nargs <- length(args <- list(...))
        if (!nargs) 
            return(as.data.frame(list()))
        if (nargs == 1 && is.list(a1 <- args[[1]])) 
            nargs <- length(args <- a1)
        if (nargs == 0) 
            return(as.data.frame(list()))
        cargs <- args
        nmc <- paste("Var", 1:nargs, sep = "")
        nm <- names(args)
        if (is.null(nm)) 
            nm <- nmc
        else if (any(ng0 <- nchar(nm) > 0)) 
            nmc[ng0] <- nm[ng0]
        names(cargs) <- nmc
        rep.fac <- 1
        d <- sapply(args, length)
        orep <- prod(d)
        for (i in seq(nargs, 1, -1)) {
            x <- args[[i]]
            nx <- length(x)
            orep <- orep/nx
            x <- x[rep.int(rep.int(seq(1, nx), rep.int(rep.fac, 
                nx)), orep)]
            if (!is.factor(x) && is.character(x)) 
                x <- factor(x, levels = unique(x))
            cargs[[i]] <- x
            rep.fac <- rep.fac * nx
        }
        rn <- as.integer(seq(1, prod(d)))
        structure(cargs, class = "data.frame", row.names = rn)
    }
    result <- NULL
    if ((nfl <- length(flist)) != 0) {
        result <- asreml.expand.grid(flist)
    }
    if (length(plist) == 0) 
        return(result)
    if (length(plist) == 1 && is.list(plist[[1]])) 
        plist <- plist[[1]]
    dim <- sapply(plist, length)
    if (all(dim == 0)) {
        result2 <- do.call("data.frame", plist)
    }
    if (!is.null(result)) 
        len <- nrow(result)
    else len <- max(dim)
    result2 <- vector("list", ndim <- length(dim))
    out.names <- names(plist)
    if (length(out.names) == 0) 
        out.names <- paste("Var", seq(length = ndim), sep = "")
    names(result2) <- out.names
    for (i in seq(length = ndim)) {
        this <- plist[[i]]
        if (isnum <- is.numeric(this)) 
            labels <- format(this)
        else {
            flabels <- levels(this)
            labels <- as.character(this)
            if (length(flabels) == 0) {
                flabels <- unique(labels)
                this <- match(this, flabels)
            }
            else this <- levelsIndex(this)
        }
        vari <- rep(this, length = len)
        result2[[i]] <- if (isnum) 
            vari
        else factor(vari, levels = seq(flabels), labels = flabels)
    }
    result2 <- data.frame(result2)
    if (!is.null(result)) 
        return(cbind(result, result2))
    else return(result2)
}
asreml.exph <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
        x <- substitute(x)
    fun <- "exph"
    type <- "h"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.expv <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
        x <- substitute(x)
    fun <- "expv"
    type <- "v"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.fa <-
function (obj, k = 1, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- "fa"
        out[[2]] <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- c(obj[[4]], paste("Comp", seq(1, k), sep = ""))
        n <- length(out[[4]]) - k
        N <- k * n + n
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, n), rep(0.1, n * k))
        }
        if (!missing(init)) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("fa(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("P", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", temp, ".var", sep = ""), paste(paste(rep(paste("!", 
            temp, sep = ""), k), "fa", sep = "."), rep(seq(1, 
            k), rep(n, k)), sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        temp[seq(1, length(con))[con == "P"]] <- 1
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n + k, 0, NA, 0)
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
        out[[13]] <- TRUE
    }
    else if (is.numeric(substitute(obj))) {
        out[[1]] <- "fa"
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        n <- obj
        N <- k * n + n
        lvls <- seq(1, obj + k)
        out[[4]] <- lvls
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, n), rep(0.1, n * k))
        }
        if (IniFlag) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("fa(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("P", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", temp, ".var", sep = ""), paste(paste(rep(paste("!", 
            temp, sep = ""), k), "fa", sep = "."), rep(seq(1, 
            k), rep(n, k)), sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(obj, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 
            0, 0, 0, 0)
        out[[10]] <- c(obj, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    else {
        out[[1]] <- "fa"
        obj <- as.character(substitute(obj))
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        drop.unused.levels <- attr(data, "Control")$drop.unused.levels
        data <- data[[obj]]
        if (is.factor(data)) 
            out[[4]] <- c(asreml.levels(data, drop.unused.levels), 
                paste("Comp", seq(1, k), sep = ""))
        else stop(paste(obj, "must be a factor\n"))
        n <- length(out[[4]]) - k
        N <- k * n + n
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0.1, n), rep(0.1, n * k))
        }
        if (IniFlag) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("fa(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("P", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", obj, ".", temp, ".var", sep = ""), 
            paste(paste(rep(paste("!", obj, ".", temp, sep = ""), 
                k), "fa", sep = "."), rep(seq(1, k), rep(n, k)), 
                sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n + k, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.facv <-
function (obj, k = 1, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to facv() must be a factor\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "facv"
    obj <- as.character(substitute(obj))
    out[[2]] <- obj
    out[[3]] <- obj
    drop.unused.levels <- attr(data, "Control")$drop.unused.levels
    data <- data[[obj]]
    if (is.factor(data)) 
        out[[4]] <- asreml.levels(data, drop.unused.levels)
    else stop(paste(obj, "must be a factor\n"))
    n <- length(out[[4]])
    N <- k * n + n
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- c(rep(0.1, n * k), rep(0.1, n))
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != N) 
            stop("facv(): Wrong number of initial values\n")
    }
    if (k > 1) 
        init[c(!lower.tri(matrix(nrow = n, ncol = k), diag = TRUE), 
            rep(FALSE, n))] <- 0
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) {
        temp <- rep("U", N)
        if (k > 1) 
            temp[c(!lower.tri(matrix(nrow = n, ncol = k), diag = TRUE), 
                rep(FALSE, n))] <- "F"
        out[[6]] <- temp
    }
    else out[[6]] <- con
    out[[7]] <- c(paste(paste(rep(paste("!", obj, ".", out[[4]], 
        sep = ""), k), "fa", sep = "."), rep(seq(1, k), rep(n, 
        k)), sep = ""), paste("!", obj, ".", out[[4]], ".var", 
        sep = ""))
    out[[8]] <- c(rep(2, n * k), rep(2, n))
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 13, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out[[10]] <- c(n, 0, NA, 0)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- obj
    out[[13]] <- TRUE
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.fak <-
function (lev, lab) 
{
    A <- length(lev)
    B <- length(lab)
    return(((A - 1) - sqrt((1 - A)^2 - 4 * (B - A)))/2)
}
asreml.findMe <-
function () 
{
    db <- find("asreml")
    where <- vector(mode = "character", length = length(db))
    for (i in seq(1, length(db))) {
        if (substring(db[i], 1, 5) == "file:") 
            where[i] <- substring(db[i], 6, nchar(db[i]))
        else if (substring(db[i], 1, 7) == "package") 
            where[i] <- system.file(package = "asreml")
        else if (db[i] == ".GlobalEnv") 
            where[i] <- getwd()
        else where[i] <- database.path(db[i])
        pthLen <- nchar(where[i])
        if (!(substring(where[i], pthLen, pthLen) == "/" | substring(where[i], 
            pthLen, pthLen) == "\\")) {
            if (substring(casefold(where[i]), (pthLen - 4), pthLen) == 
                casefold("_Data") | substring(casefold(where[i]), 
                (pthLen - 4), pthLen) == casefold(".Data")) 
                where[i] <- paste(substring(where[i], 1, pthLen - 
                  6), "/", sep = "")
            else if (substring(casefold(where[i]), (pthLen - 
                5), pthLen) == casefold(".Rdata")) 
                where[i] <- paste(substring(where[i], 1, pthLen - 
                  7), "/", sep = "")
            else where[i] <- paste(where[i], "/", sep = "")
        }
    }
    where
}
asreml.formula <-
function (object) 
{
    if (asreml.Rsys) 
        structure(object[[1]], class = "formula")
    else formula(object)
}
asreml.Fown <-
function (ftest, inter) 
{
    if (length(as.formula(ftest)) != 2) 
        stop("Expression for F test must be a formula (with no response variate)")
    if (length(ftest[[2]]) == 0) 
        return(0)
    ff <- asreml.getOp(ftest[[2]], "|")
    if (is.null(ff)) {
        ffL <- attr(terms(ftest), "term.labels")
        ffR <- character(0)
    }
    else {
        ffL <- attr(terms(formula(paste("~", deparse(ff[[2]])))), 
            "term.labels")
        tt <- terms(formula(paste("~", deparse(ff[[3]]))))
        ffR <- asreml.ie(attr(tt, "intercept"), c("(Intercept)", 
            attr(tt, "term.labels")), attr(tt, "term.labels"))
    }
    if (length(ffL) != 1) 
        stop("Only a single term can be tested in 'Ftest'")
    ownvec <- c(0, 0, match(ffL, inter$facnam), -match(ffR, inter$facnam))
    if (any(is.na(ownvec))) 
        stop("term in 'Ftest' does not appear in the model")
    ownvec[1] <- length(ownvec)
    ownvec[2] <- ownvec[1] - 2
    return(ownvec)
}
asreml.Gamma <-
function (link = "inverse", dispersion = 1, phi = 1) 
{
    link <- as.character(substitute(link))
    misnames <- c("inverse", "log", "identity", "reciprocal", 
        "1/mu", "Inverse", "Reciprocal", "Log", "Identity")
    corresp <- c(1, 2, 3, 1, 1, 1, 1, 2, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Gamma links are \"inverse\",  \"log\" or \"identity\"")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("Gamma", link = link, phi = phi)
    fam$dispersion <- dispersion
    fam$phi <- phi
    fam
}
asreml.gammas <-
function (object) 
{
    if (mode(object) != "list") 
        stop("\n object must be of mode list\n")
    if (is.na(match("G.param", names(object)))) 
        stop("\n object must have a component named G.param\n")
    if (is.na(match("R.param", names(object)))) 
        stop("\n object must have a component named R.param\n")
    gammas <- vector(mode = "numeric")
    gmcstr <- vector(mode = "character")
    gcov.list <- object[["G.param"]]
    rcov.list <- object[["R.param"]]
    if (length(gcov.list) > 0) {
        nterm <- length(gcov.list)
        for (n in seq(1, nterm)) {
            nfact <- length(gcov.list[[n]])
            for (j in 1:nfact) {
                gammas <- c(gammas, gcov.list[[n]][[j]]$initial)
                gmcstr <- c(gmcstr, gcov.list[[n]][[j]]$con)
            }
        }
    }
    if (length(rcov.list) > 0) {
        nsect <- length(rcov.list)
        nspat <- max(sapply(rcov.list, length)) - 1
        for (n in seq(1, nsect)) {
            ndim <- length(rcov.list[[n]]) - 1
            gammas <- c(gammas, rcov.list[[n]]$variance$s2)
            gmcstr <- c(gmcstr, rcov.list[[n]]$variance$con)
            for (j in 1:ndim) {
                gammas <- c(gammas, rcov.list[[n]][[j]]$initial)
                gmcstr <- c(gmcstr, rcov.list[[n]][[j]]$con)
            }
        }
    }
    if (length(gammas) == 0) 
        stop("\n No variance components found - missing G & R lists?\n")
    not.id <- !is.na(gammas)
    gammas <- gammas[not.id]
    gmcstr <- gmcstr[not.id]
    x <- data.frame(Gamma = names(gammas), Value = gammas, Constraint = gmcstr, 
        row.names = NULL, stringsAsFactors = FALSE)
    return(x)
}
asreml.gau <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
        x <- substitute(x)
    fun <- "gau"
    type <- "i"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.gauh <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(obj, "asreml.special"))) 
        x <- substitute(x)
    fun <- "gauh"
    type <- "h"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.gaussian <-
function (link = "identity", dispersion = NA) 
{
    link <- as.character(substitute(link))
    misnames <- c("inverse", "log", "identity", "reciprocal", 
        "1/mu", "Inverse", "Reciprocal", "Log", "Identity")
    corresp <- c(1, 2, 3, 1, 1, 1, 1, 2, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Gaussian links are \"log\", \"inverse\" or \"identity\"\n")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("gaussian", link = link)
    fam$dispersion <- dispersion
    fam
}
asreml.gauv <-
function (x, init = NA, dist = NA, data, ...) 
{
    if (!(mode(substitute(x)) == "call" && inherits(obj, "asreml.special"))) 
        x <- substitute(x)
    fun <- "gauv"
    type <- "v"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    out <- asreml.specialsOned(fun, x, type, data, struc, dist, 
        init)
    out
}
asreml.gdflt <-
function (gcov, data, scale = 1, control = asreml.control(...), 
    ...) 
{
    asreml.gdflt2 <- function(Y, Var, term, scale) {
        ndim <- length(Var)
        level2 <- vector(mode = "list", length = ndim)
        tmp <- NULL
        for (i in seq(1, ndim)) {
            if (length(Y) == ndim && length(unique(Var)) < length(Y)) 
                X <- Y[[i]]
            else X <- Y[[Var[i]]]
            X$Initial[X$Tgamma == 2] <- X$Initial[X$Tgamma == 
                2] * scale
            if (any(is.na(X$Initial))) 
                level2[[i]]$levels <- character(0)
            else level2[[i]]$levels <- X$Lvls
            level2[[i]]$initial <- X$Initial
            level2[[i]]$con <- X$Con
            level2[[i]]$model <- X$Fun
            names(level2[[i]]$initial) <- paste(term, X$Lab, 
                sep = "")
            tmp <- c(tmp, X$Call)
        }
        names(level2) <- tmp
        level2
    }
    Terms <- function(form) {
        tt <- terms(form, specials = asreml.Spcls$Fun, keep.order = TRUE)
        terms.fac <- attr(tt, "factors")
        which.and <- attr(tt, "specials")$and
        terms.var <- dimnames(terms.fac)[[1]]
        terms.lab <- dimnames(terms.fac)[[2]]
        terms.seq <- seq(along = terms.lab)
        if (length(which.and) > 0) {
            terms.var <- terms.var[-which.and]
            is.and <- apply((terms.fac[which.and, , drop = FALSE] > 
                0), 2, any)
            terms.seq <- seq(along = terms.lab)[!is.and]
            terms.lab <- terms.lab[-seq(along = terms.lab)[is.and]]
        }
        list(var = terms.var, lab = terms.lab, seq = terms.seq, 
            tt = tt)
    }
    strDflt <- function(G) {
        S <- vector(mode = "list", length = length(G$Lvls))
        for (s in 1:length(S)) {
            S[[s]] <- list()
            S[[s]]$Lvls <- G$Lvls[[s]]
            S[[s]]$Initial <- G$Initial[[s]]
            S[[s]]$Lab <- G$Lab[[s]]
            S[[s]]$Tgamma <- G$Tgamma[[s]]
            S[[s]]$Con <- G$Con[[s]]
            S[[s]]$Fun <- names(G$Initial)[s]
            S[[s]]$Call <- G$Call[s]
            S[[s]]$isVariance <- G$isVariance[s]
        }
        names(S) <- sapply(S, function(x) x$Fun)
        S
    }
    whichFUN <- function(Y, fname, Var) {
        sapply(Var, function(v, y, f) {
            ifelse(!is.null(y[[v]]$Fun) && (is.element(y[[v]]$Fun, 
                f) || is.element(v, f)), TRUE, FALSE)
        }, Y, fname)
    }
    fixPed <- function(Y, Var) {
        which.ped <- whichFUN(Y, c("ped", "ide", "giv", "spl", 
            "dev"), Var)
        if (sum(as.numeric(which.ped)) == 0) 
            return(Y)
        at <- whichFUN(Y, "at", Var)
        others <- sapply(Var, function(v, y) {
            ifelse(is.na(y[[v]]$isVariance) || !y[[v]]$isVariance, 
                TRUE, FALSE)
        }, Y)
        which <- (which.ped & others)
        if (all(which[!at])) {
            Y[[Var[which.ped][1]]]$isVariance <- TRUE
            if (length(which.ped) > 1) {
                for (w in seq(along = which.ped)[which.ped][-1]) {
                  Y[[Var[w]]]$Initial <- NA
                  Y[[Var[w]]]$Con <- ""
                }
            }
        }
        else if (sum(as.numeric(which))) {
            for (w in seq(along = which)[which]) {
                Y[[Var[w]]]$Initial <- NA
                Y[[Var[w]]]$Con <- ""
                Y[[Var[w]]]$isVariance <- FALSE
            }
        }
        Y
    }
    if (!inherits(gcov, "formula")) 
        stop("\nGcov must be a formula")
    if (length(gcov[[2]]) == 0) 
        return(vector(mode = "list", length = 0))
    else if (length(gcov) > 2) 
        stop("\nGcov model formula must be of the form \" ~ pred\"")
    splstp <- control$splstp
    tt <- Terms(gcov)
    terms.seq <- tt$seq
    terms.lab <- tt$lab
    terms.var <- tt$var
    nlab <- length(terms.lab)
    nvar <- length(terms.var)
    G <- lapply(terms.var, function(x, data) {
        g <- eval(asreml.Eval(x, data))
        if (is.na(match("Fun", names(g)))) {
            g <- do.call("asreml.id", list(obj = x, data = data, 
                Rcov = 0))
        }
        return(g)
    }, data)
    names(G) <- terms.var
    gdflt <- vector(mode = "list")
    Gnames <- vector(mode = "character")
    tt <- terms.order(gcov, TRUE)
    k <- 0
    for (z in terms.seq) {
        Y <- G
        Var <- dimnames(attr(tt[[z]], "factors"))[[1]]
        Term <- dimnames(attr(tt[[z]], "factors"))[[2]]
        cal <- unlist(lapply(Var, function(v, y) {
            if (inherits(y[[v]], "asreml.special") && substring(y[[v]]$Fun[1], 
                1, 2) == "id" && mode(y[[v]]$Obj) == "numeric") 
                return(NULL)
            if (inherits(y[[v]], "asreml.special") && y[[v]]$Fun[1] == 
                "grp") 
                return(y[[v]]$Obj)
            ifelse(is.null(y[[v]]$Call), v, y[[v]]$Call)
        }, Y))
        which.str <- whichFUN(Y, "str", Var)
        if (sum(as.numeric(which.str)) > 1) 
            stop("Direct product of multiple structures not allowed")
        which <- whichFUN(Y, "at", Var)
        if (sum(which) > 1) 
            stop("Ilegal use of at()\n")
        Y <- fixPed(Y, Var)
        which.id <- whichFUN(Y, c("id", "grp", "lin", "mbf"), 
            Var)
        which.cor <- sapply(Var, function(v, y) {
            ifelse(is.na(y[[v]]$isVariance) || !y[[v]]$isVariance, 
                TRUE, FALSE)
        }, Y)
        which.id[which] <- FALSE
        which.cor[which] <- FALSE
        which.idv <- which.id & which.cor
        v <- Var
        if (all(which.cor[!which])) {
            if (sum(as.numeric(which.idv))) {
                vv <- (seq(along = Var)[which.idv])[1]
                tmp <- Y[[Var[vv]]]$Obj
                Y[[Var[vv]]] <- do.call("asreml.idv", list(obj = tmp, 
                  data = data, Rcov = 0))
            }
            else warning(paste("Term", Term, "is a correlation model.\n"))
        }
        if (any(which.str)) {
            ws <- seq(along = Var)[which.str]
            S <- strDflt(Y[[Var[ws]]])
            S <- fixPed(S, names(S))
            str.tt <- Terms(formula(paste("~", Y[[Var[ws]]]$Obj)))
            str.var <- str.tt$var
            str.lab <- str.tt$lab
            str.seq <- str.tt$seq
            at.form <- character(0)
            for (zz in str.seq) {
                str.var <- dimnames(attr(str.tt$tt[zz], "factors"))[[1]]
                str.term <- dimnames(attr(str.tt$tt[zz], "factors"))[[2]]
                sG <- lapply(str.var, function(x, data) {
                  sg <- eval(asreml.Eval(x, data))
                  if (is.na(match("Fun", names(sg)))) {
                    sg <- do.call("asreml.id", list(obj = x, 
                      data = data, Rcov = 0))
                  }
                  return(sg)
                }, data)
                names(sG) <- sapply(sG, function(x) {
                  if (x$Fun == "grp") 
                    x$Obj
                  else x$Call
                })
                str.var <- names(sG)
                str.lab <- paste(str.var, collapse = ":")
                str.at <- whichFUN(sG, "at", str.var)
                if (any(str.at)) {
                  at.form <- unique(c(at.form, paste(unlist(lapply(sG[str.var[str.at]], 
                    function(x) {
                      paste("at(", x$Obj, ", ", x$Lvls, ")", 
                        sep = "")
                    })), str.var[!str.at], sep = ":")))
                }
                else at.form <- unique(c(at.form, str.lab))
            }
            at.form <- paste(at.form, collapse = "+")
            k <- k + 1
            gdflt[[k]] <- asreml.gdflt2(S, names(S), at.form, 
                scale)
            Gnames[k] <- at.form
        }
        else if (sum(which)) {
            X <- Y[[Var[which]]]
            if (sum(as.numeric(!(which | which.id))) > 0) {
                for (what in Var[!(which | which.id)]) {
                  if (!inherits(Y[[what]], "asreml.special")) {
                    if (any(tapply(data[[Y[[what]]$Obj]], data[[X$Obj]], 
                      function(x) {
                        length(unique(x))
                      }) != length(unique(data[[Y[[what]]$Obj]])))) 
                      warning(paste("Not all levels of", Y[[what]]$Obj, 
                        "appear in each level of", X$Obj, "\n"))
                  }
                }
            }
            for (lvl in X$Lvls) {
                cal[which] <- paste("at(", X$Obj, ", ", lvl, 
                  ")", sep = "")
                term <- paste(cal, collapse = ":")
                k <- k + 1
                gdflt[[k]] <- asreml.gdflt2(Y, Var[!which], term, 
                  scale)
                Gnames[k] <- term
            }
        }
        else {
            term <- paste(cal, collapse = ":")
            k <- k + 1
            gdflt[[k]] <- asreml.gdflt2(Y, Var, term, scale)
            Gnames[k] <- term
        }
    }
    names(gdflt) <- Gnames
    return(gdflt)
}
asreml.generate <-
function (lvls) 
{
    N <- length(lvls)
    x <- vector(length = (N + 2))
    x[1] <- x[N + 2] <- 1
    x[2:(N + 1)] <- sapply(lvls, function(x) {
        length(x)
    })
    indx <- matrix("", nrow = prod(x[2:(N + 1)]), ncol = N)
    for (i in seq(2, (N + 1))) indx[, i - 1] <- rep(rep(lvls[[i - 
        1]], rep(prod(x[seq(i + 1, N + 2)]), x[i])), prod(x[seq(1, 
        i - 1)]))
    indx
}
asreml.getData <-
function (d = "data") 
{
    data <- NULL
    for (i in 1:length(sys.frames())) {
        if (exists(d, envir = sys.frame(-i), inherits = FALSE)) {
            data <- get(d, envir = sys.frame(-i), inherits = FALSE)
            break
        }
    }
    data
}
asreml.getFamily <-
function (thing) 
{
    if (is.null(thing)) 
        family <- asreml.gaussian()
    else if (is.character(thing)) 
        family <- eval(parse(text = thing, sep = ""))
    else if (is.name(thing)) 
        family <- eval(parse(text = paste(as.character(thing), 
            "()", sep = "")))
    else if (is.function(thing)) 
        family <- thing()
    else family <- eval(thing)
    family
}
asreml.getInterExp <-
function (expr) 
{
    if (is.name(expr) || !is.language(expr)) 
        return(NULL)
    if (expr[[1]] == as.name("(")) 
        return(asreml.getInterExp(expr[[2]]))
    if (!is.call(expr)) 
        stop("expr must be of class call")
    if (expr[[1]] == as.name(":")) 
        return(expr)
    if (length(expr) == 2) 
        return(asreml.getInterExp(expr[[2]]))
    c(asreml.getInterExp(expr[[2]]), asreml.getInterExp(expr[[3]]))
}
asreml.getMbf <-
function (mbfr, data) 
{
    nn <- names(mbfr)
    mbfAttr <- list()
    mbFrames <- list()
    mbfX <- vector(mode = "numeric")
    for (i in nn) {
        mdf <- eval(as.name(mbfr[[i]]$dataFrame))
        minus <- match(mbfr[[i]]$key[[2]], names(mdf))
        colNames <- names(mdf)[-minus]
        colNames <- c(mbfr[[i]]$key[1], colNames)
        k1 <- mbfr[[i]]$key[1]
        k2 <- mbfr[[i]]$key[2]
        ky1 <- unique(as.character(data[, k1]))
        idx <- match(ky1, as.character(mdf[, k2]))
        if (any(which <- is.na(idx))) {
            mdf <- as.list(mdf)
            mdf[[1]] <- c(mdf[[1]], ky1[which])
            for (i in seq(2, length(mdf))) mdf[[i]] <- rep(NA, 
                length(which[which]))
            mdf <- data.frame(mdf)
        }
        what <- match(as.character(mdf[, k2]), as.character(data[, 
            k1]))
        if (is.factor(data[, k1])) 
            mdf[, k2] <- as.numeric(data[what, k1])
        else if (is.numeric(data[, k1])) 
            mdf[, k2] <- data[what, k1]
        else mdf[, k2] <- as.numeric(factor(data[, k1]))
        srt <- order(mdf[, k2])
        mbfAttr[[i]]$colNames <- colNames
        mbfAttr[[i]]$key <- mbfr[[i]]$key
        mbfAttr[[i]]$ncz <- length(colNames) - 1
        mbfAttr[[i]]$nux <- nrow(mdf)
        kk <- seq(1, ncol(mdf))
        kk[seq(1, minus)] <- kk[seq(1, minus)] - 1
        kk[1] <- minus
        mbFrames[[i]] <- mdf[srt, kk]
        mbfX <- c(mbfX, as.vector(as.matrix(mdf[srt, kk])))
    }
    mbfX[is.na(mbfX)] <- 0
    names(mbFrames) <- nn
    list(mbfAttr = mbfAttr, mbFrames = mbFrames, mbfX = mbfX)
}
asreml.getOp <-
function (expr, op) 
{
    if (is.name(expr) || !is.language(expr)) 
        return(NULL)
    if (expr[[1]] == as.name("(")) 
        return(asreml.getOp(expr[[2]], op))
    if (!is.call(expr)) 
        stop("expr must be of class call")
    if (expr[[1]] == as.name(op)) 
        return(expr)
    if (length(expr) == 2) 
        return(asreml.getOp(expr[[2]], op))
    c(asreml.getOp(expr[[2]], op), asreml.getOp(expr[[3]], op))
}
asreml.getPdim <-
function (dataFrame, obj) 
{
    dims <- attr(dataFrame, "Control")$splXtras$dimTbl
    which <- match(obj, dims$axis)
    if (is.na(which)) 
        return(0)
    term <- dims$term[which]
    ndim <- sum(as.numeric(as.character(dims$term) == term))
    return(ndim)
}
asreml.getPpoints <-
function (dataFrame, obj) 
{
    splXtras <- attr(dataFrame, "Control")$splXtras
    dims <- splXtras$dimTbl
    pp <- vector(mode = "double", length = 0)
    which <- match(obj, dims$axis)
    if (is.na(which)) 
        return(pp)
    term <- dims$term[which]
    what <- (as.character(dims$term) == term)
    N <- match(obj, dims$axis[what])
    if (N == 1) {
        pp <- splXtras$points[seq(splXtras$IBV[2, term], splXtras$IBV[2, 
            term] + splXtras$IBV[4, term] - 1)]
    }
    else if (N == 2) {
        pp <- splXtras$points[seq(splXtras$JBV[2, term], splXtras$JBV[2, 
            term] + splXtras$JBV[4, term] - 1)]
    }
    return(pp)
}
asreml.ginv <-
function (x, tol = 1e-08) 
{
    S <- svd(x)
    D <- S$d
    delcol <- abs(D) > tol
    D <- D[delcol]
    V <- S$v[, delcol]
    U <- S$u[, delcol]
    t(U %*% (t(V)/D))
}
asreml.giv <-
function (obj, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- "giv"
        out[[2]] <- obj[[3]]
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        IniFlag <- FALSE
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- 0.1
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != 1) 
                stop("giv(): Wrong number of initial values\n")
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- "P"
        else out[[6]] <- con
        out[[7]] <- paste("!", obj[[2]], ".giv", sep = "")
        out[[8]] <- 2
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, -7, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        names(out[[9]]) <- c("", "", out[[3]], rep("", 13))
        out[[10]] <- obj[[10]]
        attr(out[[10]], "Hack") <- obj[[1]]
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
    }
    else {
        out[[1]] <- "giv"
        obj <- as.character(substitute(obj))
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        v <- match(obj, names(data))
        if (!is.factor(data[[obj]])) 
            stop("Argument to giv() must be a factor)\n")
        gv <- attr(data, "Control")$ginverse
        if (is.null(gv$rownames)) 
            stop("Deprecated 'ginverse' component of the fitted object, rerun with a newer version of asreml()")
        glvls <- gv$rownames(obj, gv)
        if (length(flvls <- levels(data[[obj]])) == 0) 
            flvls <- unique(data[[obj]])
        if (any(is.na(match(flvls, glvls)))) 
            stop(paste(obj, "has levels in data missing in ginverse\n"))
        out[[4]] <- glvls
        IniFlag <- FALSE
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- 0.1
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != 1) 
                stop("giv(): Wrong number of initial values\n")
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- "P"
        else out[[6]] <- con
        out[[7]] <- ".giv"
        out[[8]] <- 2
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, -7, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        names(out[[9]]) <- c("", "", obj, rep("", 13))
        out[[10]] <- c(-2, 0, v, NA)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
    }
    out[[13]] <- ifelse(IniFlag, TRUE, FALSE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.glist <-
function (gcov, data) 
{
    asreml.glist2 <- function(Y, Var, is.str = FALSE) {
        ndim <- length(Var)
        level2 <- vector(mode = "list", length = ndim)
        tmp <- NULL
        for (i in seq(1, ndim)) {
            if (length(Y) == ndim && length(unique(Var)) < length(Y)) 
                X <- Y[[i]]
            else X <- Y[[Var[i]]]
            level2[[i]]$order <- length(X$Lvls)
            level2[[i]]$struc <- X$Struc
            level2[[i]]$tgamma <- abs(X$Tgamma)
            level2[[i]]$coords <- X$Coords$coords
            level2[[i]]$obj <- X$Obj
            level2[[i]]$is.str <- is.str
            tmp <- c(tmp, X$Call)
        }
        names(level2) <- tmp
        return(level2)
    }
    Terms <- function(form) {
        tt <- terms(form, specials = asreml.Spcls$Fun, keep.order = TRUE)
        terms.fac <- attr(tt, "factors")
        which.and <- attr(tt, "specials")$and
        terms.var <- dimnames(terms.fac)[[1]]
        terms.lab <- dimnames(terms.fac)[[2]]
        terms.seq <- seq(along = terms.lab)
        if (length(which.and) > 0) {
            terms.var <- terms.var[-which.and]
            is.and <- apply((terms.fac[which.and, , drop = FALSE] > 
                0), 2, any)
            terms.seq <- seq(along = terms.lab)[!is.and]
            terms.lab <- terms.lab[-seq(along = terms.lab)[is.and]]
        }
        list(var = terms.var, lab = terms.lab, seq = terms.seq, 
            tt = tt)
    }
    strList <- function(G) {
        S <- vector(mode = "list", length = length(G$Lvls))
        for (s in 1:length(S)) {
            S[[s]] <- list()
            S[[s]]$Lvls <- G$Lvls[[s]]
            S[[s]]$Initial <- G$Initial[[s]]
            S[[s]]$Lab <- G$Lab[[s]]
            S[[s]]$Tgamma <- G$Tgamma[[s]]
            S[[s]]$Con <- G$Con[[s]]
            S[[s]]$Fun <- names(G$Initial)[s]
            S[[s]]$Call <- G$Call[s]
            S[[s]]$Struc <- G$Struc[[s]]
            S[[s]]$Coords <- G$Coords
            S[[s]]$Obj <- G$Obj[1]
            S[[s]]$isVariance <- G$isVariance[s]
        }
        names(S) <- sapply(S, function(x) x$Fun)
        S
    }
    whichFUN <- function(Y, fname, Var) {
        sapply(Var, function(v, y, f) {
            ifelse(!is.null(y[[v]]$Fun) && (is.element(y[[v]]$Fun, 
                f) || is.element(v, f)), TRUE, FALSE)
        }, Y, fname)
    }
    fixPed <- function(Y, Var) {
        which.ped <- whichFUN(Y, c("ped", "ide", "giv", "spl", 
            "dev"), Var)
        if (sum(as.numeric(which.ped)) == 0) 
            return(Y)
        at <- whichFUN(Y, "at", Var)
        others <- sapply(Var, function(v, y) {
            ifelse(is.na(y[[v]]$isVariance) || !y[[v]]$isVariance, 
                TRUE, FALSE)
        }, Y)
        which <- (which.ped & others)
        if (all(which[!at])) {
            Y[[Var[which.ped][1]]]$isVariance <- TRUE
            if (length(which.ped) > 1) {
                for (w in seq(along = which.ped)[which.ped][-1]) {
                  Y[[Var[w]]]$Initial <- NA
                  Y[[Var[w]]]$Con <- ""
                }
            }
        }
        else if (sum(as.numeric(which))) {
            for (w in seq(along = which)[which]) {
                Y[[Var[w]]]$Initial <- NA
                Y[[Var[w]]]$Con <- ""
                Y[[Var[w]]]$isVariance <- FALSE
            }
        }
        Y
    }
    if (length(gcov[[2]]) == 0) 
        return(vector(mode = "list", length = 0))
    tt <- Terms(gcov)
    terms.seq <- tt$seq
    terms.lab <- tt$lab
    terms.var <- tt$var
    nlab <- length(terms.lab)
    nvar <- length(terms.var)
    G <- lapply(terms.var, function(x, data) {
        g <- eval(asreml.Eval(x, data))
        if (is.na(match("Fun", names(g)))) {
            g <- do.call("asreml.id", list(obj = x, data = data, 
                Rcov = 0))
        }
        return(g)
    }, data)
    names(G) <- terms.var
    glist <- vector(mode = "list")
    Gnames <- vector(mode = "character")
    tt <- terms.order(gcov, TRUE)
    k <- 0
    for (z in terms.seq) {
        Y <- G
        Var <- dimnames(attr(tt[[z]], "factors"))[[1]]
        Term <- dimnames(attr(tt[[z]], "factors"))[[2]]
        cal <- unlist(lapply(Var, function(v, y) {
            if (inherits(y[[v]], "asreml.special") && substring(y[[v]]$Fun[1], 
                1, 2) == "id" && mode(y[[v]]$Obj) == "numeric") 
                return(NULL)
            if (inherits(y[[v]], "asreml.special") && y[[v]]$Fun[1] == 
                "grp") 
                return(y[[v]]$Obj)
            ifelse(is.null(y[[v]]$Call), v, y[[v]]$Call)
        }, Y))
        which.str <- whichFUN(Y, "str", Var)
        if (sum(as.numeric(which.str)) > 1) 
            stop("Direct product of multiple structures not allowed")
        which <- whichFUN(Y, "at", Var)
        if (sum(which) > 1) 
            stop("Ilegal use of at()\n")
        Y <- fixPed(Y, Var)
        which.id <- whichFUN(Y, c("id", "grp", "lin", "mbf"), 
            Var)
        which.cor <- sapply(Var, function(v, y) {
            ifelse(is.na(y[[v]]$isVariance) || !y[[v]]$isVariance, 
                TRUE, FALSE)
        }, Y)
        which.id[which] <- FALSE
        which.cor[which] <- FALSE
        which.idv <- which.id & which.cor
        v <- Var
        if (all(which.cor[!which])) {
            if (sum(as.numeric(which.idv))) {
                vv <- (seq(along = Var)[which.idv])[1]
                tmp <- Y[[Var[vv]]]$Obj
                Y[[Var[vv]]] <- do.call("asreml.idv", list(obj = tmp, 
                  data = data, Rcov = 0))
            }
            else warning(paste("Term", Term, "is a correlation structure.\n"))
        }
        if (any(which.str)) {
            ws <- seq(along = Var)[which.str]
            S <- strList(Y[[Var[ws]]])
            S <- fixPed(S, names(S))
            str.tt <- Terms(formula(paste("~", Y[[Var[ws]]]$Obj)))
            str.var <- str.tt$var
            str.lab <- str.tt$lab
            str.seq <- str.tt$seq
            at.form <- character(0)
            for (zz in str.seq) {
                str.var <- dimnames(attr(str.tt$tt[zz], "factors"))[[1]]
                str.term <- dimnames(attr(str.tt$tt[zz], "factors"))[[2]]
                sG <- lapply(str.var, function(x, data) {
                  sg <- eval(asreml.Eval(x, data))
                  if (is.na(match("Fun", names(sg)))) {
                    sg <- do.call("asreml.id", list(obj = x, 
                      data = data, Rcov = 0))
                  }
                  return(sg)
                }, data)
                names(sG) <- sapply(sG, function(x) {
                  if (x$Fun == "grp") 
                    x$Obj
                  else x$Call
                })
                str.var <- names(sG)
                str.lab <- paste(str.var, collapse = ":")
                str.at <- whichFUN(sG, "at", str.var)
                if (any(str.at)) {
                  at.form <- unique(c(at.form, paste(unlist(lapply(sG[str.var[str.at]], 
                    function(x) {
                      paste("at(", x$Obj, ", ", x$Lvls, ")", 
                        sep = "")
                    })), str.var[!str.at], sep = ":")))
                }
                else at.form <- unique(c(at.form, str.lab))
            }
            at.form <- paste(at.form, collapse = "+")
            k <- k + 1
            glist[[k]] <- asreml.glist2(S, names(S), is.str = TRUE)
            Gnames[k] <- at.form
        }
        else if (sum(which)) {
            X <- Y[[Var[which]]]
            for (lvl in X$Lvls) {
                cal[which] <- paste("at(", X$Obj, ", ", lvl, 
                  ")", sep = "")
                term <- paste(cal, collapse = ":")
                k <- k + 1
                glist[[k]] <- asreml.glist2(Y, Var[!which])
                Gnames[k] <- term
            }
        }
        else {
            term <- paste(cal, collapse = ":")
            k <- k + 1
            glist[[k]] <- asreml.glist2(Y, Var)
            Gnames[k] <- term
        }
    }
    names(glist) <- Gnames
    return(glist)
}
asreml.glm <-
function (fm) 
{
    asr.id <- c(1, 2, 3, 4, -1, -1, 5, 9, 1, 2, 3, 4, -1, -1, 
        5, 9)
    names(asr.id) <- c("Gaussian", "Binomial", "Poisson", "Gamma", 
        "Inverse Gaussian", "Log Normal", "Negative Binomial", 
        "Multinomial", "gaussian", "binomial", "poisson", "gamma", 
        "inverse.gaussian", "log.normal", "negative.binomial", 
        "multinomial")
    asr.link <- c(2, 3, 4, 1, 7, 5, -1, 6, 2, 3, 4, 1, 7, 5, 
        -1, 6)
    names(asr.link) <- c("Logit", "Probit", "Complementary Log", 
        "Identity", "Inverse", "Log", "Inverse Square", "Square Root", 
        "logit", "probit", "cloglog", "identity", "inverse", 
        "log", "1/mu^2", "sqrt")
    dist <- fm$family[1]
    lnk <- substring(fm$family[2], 1, match(":", substring(fm$family[2], 
        1:nchar(fm$family[2]), 1:nchar(fm$family[2]))) - 1)
    id <- asr.id[match(dist, names(asr.id))]
    link <- asr.link[match(lnk, names(asr.link))]
    if (is.null(fm$dispersion)) 
        dispersion <- ifelse(id == 1, NA, 1)
    else dispersion <- fm$dispersion
    names(dispersion) <- names(id)
    phi <- ifelse(is.null(fm$phi), 1, fm$phi)
    list(id = id, link = link, dispersion = dispersion, phi = phi)
}
asreml.glmDeriv <-
function (link) 
{
    switch(link, identity = function(mu) {
        1
    }, logit = function(mu) {
        d <- mu * (1 - mu)
        if (any(tiny <- (d < .Machine$double.eps))) {
            warning("Model unstable; fitted probabilities of 0 or 1")
            d[tiny] <- .Machine$double.eps
        }
        1/d
    }, cloglog = function(mu) {
        mu2 <- 1 - mu
        d <- -(mu2 * log(mu2))
        if (any(tiny <- (d < .Machine$double.eps))) {
            warning("Model unstable; fitted probabilities of 0 or 1")
            d[tiny] <- .Machine$double.eps
        }
        1/d
    }, probit = function(mu) {
        sqrt(2 * pi) * exp((qnorm(mu)^2)/2)
    }, log = function(mu) {
        1/mu
    }, inverse = function(mu) {
        -1/mu^2
    }, `1/mu^2` = function(mu) {
        -2/mu^3
    }, sqrt = function(mu) {
        1/(2 * sqrt(mu))
    })
}
asreml.glmDeviance <-
function (variance, residuals = FALSE) 
{
    switch(variance, constant = function(mu, y, w, residuals = substitute(residuals), 
        phi = NULL) {
        if (residuals) sqrt(w) * (y - mu) else sum(w * (y - mu)^2)
    }, `mu(1-mu)` = function(mu, y, w, residuals = substitute(residuals), 
        phi = NULL) {
        devy <- y
        nz <- y != 0
        devy[nz] <- y[nz] * log(y[nz])
        nz <- (1 - y) != 0
        devy[nz] <- devy[nz] + (1 - y[nz]) * log(1 - y[nz])
        devmu <- y * log(mu) + (1 - y) * log(1 - mu)
        if (any(small <- mu * (1 - mu) < .Machine$double.eps)) {
            warning("fitted values close to 0 or 1")
            smu <- mu[small]
            sy <- y[small]
            smu <- ifelse(smu < .Machine$double.eps, .Machine$double.eps, 
                smu)
            onemsmu <- ifelse((1 - smu) < .Machine$double.eps, 
                .Machine$double.eps, 1 - smu)
            devmu[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
        }
        devi <- 2 * (devy - devmu)
        if (residuals) sign(y - mu) * sqrt(abs(devi) * w) else sum(w * 
            devi)
    }, mu = function(mu, y, w, residuals = substitute(residuals), 
        phi = NULL) {
        nz <- y > 0
        devi <- -(y - mu)
        devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else 2 * 
            sum(w * devi)
    }, `mu^2` = function(mu, y, w, residuals = substitute(residuals), 
        phi = NULL) {
        nz <- y > 0
        devi <- (y - mu)/mu + log(mu)
        devi[nz] <- devi[nz] - log(y[nz])
        if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else 2 * 
            sum(w * devi)
    }, `mu^3` = function(mu, y, w, residuals = substitute(residuals), 
        phi = NULL) {
        devi <- ((y - mu)^2)/(mu^2 * y)
        if (residuals) sign(y - mu) * sqrt(w * abs(devi)) else sum(w * 
            devi)
    }, `mu+mu^2/phi` = function(mu, y, w, residuals = substitute(residuals), 
        phi = 1) {
        devi <- 2 * w * (y * log(pmax(1, y)/mu) - (y + phi) * 
            log((y + phi)/(mu + phi)))
        if (residuals) sign(y - mu) * sqrt(abs(devi)) else sum(devi)
    })
}
asreml.gmcon <-
function (gammas, constraints) 
{
    ngmcon <- 0
    ngamma <- length(gammas)
    gequal <- rep(0, ngamma)
    gmcon <- rep(0, ngamma)
    k <- 1
    if (!is.null(constraints)) {
        which <- match(dimnames(constraints)[[1]], names(gammas))
        if (sum(is.numeric(is.na(which))) > 0) 
            stop("\nCannot match rownames of constraint matrix with gamma vector\n")
        ngmcon <- ncol(constraints)
        gmcon <- matrix(0, nrow = ngamma, ncol = ngmcon)
        zeros <- rep(0, nrow(constraints))
        ones <- rep(1, nrow(constraints))
        k <- 0
        base <- min(which) - 1
        for (i in 1:ncol(constraints)) {
            mask1 <- constraints[, i] == ones
            mask0 <- constraints[, i] == zeros
            if (sum(as.numeric(mask1)) > 0 || sum(as.numeric(!mask1 & 
                !mask0)) > 0) {
                k <- k + 1
                gequal[which][!mask0] <- which[!mask0][1]
                gmcon[which, k] <- constraints[, i]
            }
        }
        gmcon <- rowSums(gmcon)
    }
    list(gmcon = gmcon, gequal = gequal, ngmcon = max(gequal))
}
asreml.Gorder <-
function (random, data, ran.order) 
{
    if (length(random[[2]]) == 0) 
        return(random)
    keep.order <- switch(ran.order, user = TRUE, R = FALSE, noeff = FALSE, 
        stop("ran.order must be one of 'user','R','noeff'"))
    if (keep.order) 
        return(random)
    if (ran.order == "R") 
        return(formula(paste("~", paste(attr(terms(random, keep.order = FALSE), 
            "term.labels"), collapse = "+"))))
    tt <- terms(random, specials = asreml.Spcls$Fun)
    tt.fac <- attr(tt, "factors")
    tt.var <- dimnames(tt.fac)[[1]]
    tt.lab <- dimnames(tt.fac)[[2]]
    which.spc <- unique(unlist(lapply(attr(tt, "specials"), function(x) x)))
    is.spc <- rep(FALSE, length(tt.var))
    is.spc[which.spc] <- TRUE
    names(is.spc) <- tt.var
    L <- unlist(sapply(tt.var, function(x, is.spc, data) {
        if (is.spc[x]) 
            length(eval(asreml.Eval(x, data))$Lvls)
        else if (is.factor(data[[x]])) 
            length(levels(data[[x]]))
        else 1
    }, is.spc, data))
    names(L) <- tt.var
    noeff <- vector(mode = "numeric", length = length(tt.lab))
    names(noeff) <- tt.lab
    for (i in 1:length(tt.lab)) noeff[i] <- prod(L[tt.fac[, i] > 
        0])
    return(formula(paste("~", paste(names(noeff)[order(noeff)], 
        collapse = "+"))))
}
asreml.grep <-
function (pattern, text) 
{
    seq(1, length(text))[regexpr(pattern, text) != -1]
}
asreml.grp <-
function (obj, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to grp() must be a numeric or character vector\n")
    out <- vector(mode = "list", length = 13)
    obj <- as.character(substitute(obj))
    out[[1]] <- "grp"
    out[[2]] <- asreml.spCall(sys.call())
    out[[3]] <- obj
    what <- match(obj, names(attr(data, "GROUP")))
    if (is.na(what)) 
        stop(paste("Object", obj, "is not a component of the group argument"))
    which <- attr(data, "GROUP")[[obj]]
    out[[4]] <- which
    out[[5]] <- NA
    out[[6]] <- ""
    out[[7]] <- ""
    out[[8]] <- 2
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(length(which), 0, NA, 0)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- which
    out[[13]] <- NA
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.gupdt <-
function (gcov.param, gammas, gmcstr = NULL, VCov = FALSE) 
{
    if (length(gcov.param) == 0) 
        return(NULL)
    gupdt <- lapply(gcov.param, function(x, gammas, gmcstr, VCov) {
        lapply(x, function(y, gammas, gmcstr, VCov) {
            tmp <- names(y$initial)
            i <- match(tmp, names(gammas), nomatch = NA)
            ii <- !is.na(i)
            i <- i[ii]
            if (length(i) > 0) {
                y$initial[ii] <- gammas[i]
                if (is.null(y$model)) 
                  stop("G.param 'model' component missing - rerun\n")
                if (VCov) {
                  y$initial <- switch(y$model, ante = asreml.udu(y$initial), 
                    chol = asreml.ldl(y$initial), cholc = asreml.ldl(y$initial), 
                    y$initial)
                  if (!is.null(gmcstr)) {
                    if (is.element(y$model, c("ante", "chol", 
                      "cholc"))) 
                      y$con <- y$con
                    else y$con[ii] <- casefold(substr(gmcstr[i], 
                      1, 1), upper = TRUE)
                  }
                }
                else if (!is.null(gmcstr)) 
                  y$con[ii] <- casefold(substr(gmcstr[i], 1, 
                    1), upper = TRUE)
                names(y$initial) <- tmp
            }
            y
        }, gammas, gmcstr, VCov)
    }, gammas, gmcstr, VCov)
    gupdt
}
asreml.guzpfx <-
function (rgmstr) 
{
    guzpfx <- c(" ", "Positive", "?", "Unconstrained", "Fixed", 
        "Constrained", "Singular", "Boundary")
    pc <- apply(matrix((rgmstr + 1), ncol = 1), 1, function(x) {
        if (is.na(x)) 
            return(3)
        y <- max(x, 1)
        if (y > 8 && y%%10 == 2) 
            y <- 3
        if (y > 8) 
            y <- 8
        y
    })
    guzpfx[pc]
}
asreml.id <-
function (obj, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("id", "cor", "numeric(0)", "numeric(0)", 
        NA, obj, NA, data, match.call()$Rcov)
    out
}
asreml.ide <-
function (obj, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- obj[[1]]
        out[[2]] <- paste("ide(", obj[[2]], ")", sep = "")
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        IniFlag <- FALSE
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- 0.1
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != 1) 
                stop("ide(): Wrong number of initial values\n")
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- "P"
        else out[[6]] <- con
        out[[7]] <- ""
        out[[8]] <- 2
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(-2, -1, NA, 0)
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
    }
    else {
        out[[1]] <- "ide"
        obj <- as.character(substitute(obj))
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        v <- match(obj, names(data))
        if (!is.factor(data[[obj]])) 
            stop("Argument to ped() must be a factor)\n")
        gv <- attr(data, "Control")$ginverse
        if (is.null(gv$rownames)) 
            stop("Deprecated 'ginverse' component of the fitted object, rerun with a newer version of asreml()")
        glvls <- gv$rownames(obj, gv)
        if (length(flvls <- levels(data[[obj]])) == 0) 
            flvls <- unique(data[[obj]])
        if (any(is.na(match(flvls, glvls)))) {
            print(flvls[is.na(match(flvls, glvls))])
            stop(paste(obj, "has levels in data missing in ginverse\n"))
        }
        out[[4]] <- glvls
        IniFlag <- FALSE
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- 0.1
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != 1) 
                stop("ide(): Wrong number of initial values\n")
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- "P"
        else out[[6]] <- con
        out[[7]] <- "!id"
        out[[8]] <- 2
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(-2, -1, v, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
    }
    out[[13]] <- ifelse(IniFlag, TRUE, FALSE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.idh <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("idh", "cor", "numeric(0)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.idv <-
function (obj, init = NA, data, Rcov, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("idv", "cor", "numeric(0)", 1, NA, obj, 
        init, data, match.call()$Rcov)
    out
}
asreml.ie <-
function (test, altT, altF) 
{
    if (test) 
        altT
    else altF
}
asreml.ieuc <-
function (x, y, init = NA, data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "ieuc"
    type <- "ii"
    struc <- c(0, 0, 6, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, dist = dist))
    out
}
asreml.ieuch <-
function (x, y, init = NA, data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "ieuch"
    type <- "ih"
    struc <- c(0, 0, 6, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, dist = dist))
    out
}
asreml.ieucv <-
function (x, y, init = NA, data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "ieucv"
    type <- "iv"
    struc <- c(0, 0, 6, 0, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, dist = dist))
    out
}
asreml.iexp <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "iexp"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ii"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.iexph <-
function (x, y, init = NA, data, ...) 
{
    fun <- "iexph"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ih"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.iexpv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "iexpv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "iv"
    struc <- c(0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.igau <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "igau"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ii"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.igauh <-
function (x, y, init = NA, data, ...) 
{
    fun <- "igauh"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ih"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.igauv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "igauv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "iv"
    struc <- c(0, 0, 6, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.iniCon <-
function (x) 
{
    n <- nchar(x)
    vec <- casefold(substring(x, 1:n, 1:n), upper = TRUE)
    where <- match(c("P", "U", "F"), vec)
    if (sum(!is.na(where)) == 0) {
        x <- paste(x, "F", sep = "")
        where <- c(NA, NA, n + 1)
    }
    if (sum(!is.na(where)) > 1) 
        stop("Only one of 'P','U','F' may be specified in 'mtrn'")
    which <- where[!is.na(where)]
    value <- as.numeric(substring(x, 1, which - 1))
    constraint <- casefold(substring(x, which, which), upper = TRUE)
    list(value = value, constraint = constraint)
}
asreml.inter <-
function (data, fixed, random, sparse, cinv, model.y, fix.ord, 
    locgiv, glm, multivariate) 
{
    rn <- function(str) {
        match(str, c("inter1", "inter2", "inter3", "inter4", 
            "flev", "nlev", "inform", "facord", "faconst"))
    }
    ntrt <- length(model.y)
    facnam <- vector(mode = "character")
    baseFac <- vector(mode = "list", length = 3)
    names(baseFac) <- c("facNum", "baseObj", "baseFun")
    inter <- matrix(0, nrow = 9, ncol = ntrt)
    varLevels <- vector(mode = "list")
    IBV <- attr(data, "Control")$splXtras$IBV
    JBV <- attr(data, "Control")$splXtras$JBV
    coords <- list()
    tt <- terms(fixed, specials = asreml.Spcls$Fun)
    terms.fac <- attr(tt, "factors")
    terms.var <- attr(tt, "variables")
    terms.rsp <- attr(tt, "response")
    terms.mu <- attr(tt, "intercept")
    inter[rn("nlev"), 1:ntrt] <- rep(1, ntrt)
    inter[rn("flev"), 1:ntrt] <- rep(1, ntrt)
    inter[rn("inter1"), 1:ntrt] <- rep(1, ntrt)
    inter[rn("inter3"), 1:ntrt] <- seq(1, ntrt)
    if (length(terms.fac) == 0) {
        facnam <- as.character(terms.var)
        if (asreml.Rsys) 
            facnam <- facnam[-1]
    }
    else facnam <- model.y
    inter[c(1, 2, 3, 4, 7, 8), 1] <- c(1, 0, 1, 0, 1, 1)
    if (ntrt == 1 & terms.mu & glm$id != 9 & !multivariate) {
        inter <- cbind(inter, c(-8, 0, 0, 0, 1, 0, 1, inter[rn("facord"), 
            ncol(inter)] + 1, 0))
        facnam <- c(facnam, "(Intercept)")
        baseFac$facNum <- c(baseFac$facNum, length(facnam))
        baseFac$baseFun <- c(baseFac$baseFun, "(Intercept)")
        baseFac$baseObj <- c(baseFac$baseObj, "(Intercept)")
        varLevels$"(Intercept)" <- "(Intercept)"
    }
    if (length(terms.fac) > 0) {
        fixed <- asreml.formula(parse(text = paste("~", as.character(fixed[3]))))
        asr.inter <- asreml.terms(data, fixed, facnam, varLevels, 
            baseFac, inter, IBV, JBV, coords, 1, fix.ord, locgiv)
        facnam <- asr.inter$facnam
        baseFac <- asr.inter$baseFac
        inter <- asr.inter$inter
        varLevels <- asr.inter$varLevels
        IBV <- asr.inter$IBV
        JBV <- asr.inter$JBV
        coords <- asr.inter$coords
    }
    if (length(random[[2]]) > 0) {
        keep.order <- TRUE
        asr.inter <- asreml.terms(data, random, facnam, varLevels, 
            baseFac, inter, IBV, JBV, coords, 2, keep.order, 
            locgiv)
        facnam <- asr.inter$facnam
        baseFac <- asr.inter$baseFac
        inter <- asr.inter$inter
        varLevels <- asr.inter$varLevels
        IBV <- asr.inter$IBV
        JBV <- asr.inter$JBV
        coords <- asr.inter$coords
    }
    if (length(sparse[[2]]) > 0) {
        keep.order <- FALSE
        asr.inter <- asreml.terms(data, sparse, facnam, varLevels, 
            baseFac, inter, IBV, JBV, coords, 3, keep.order, 
            locgiv)
        facnam <- asr.inter$facnam
        baseFac <- asr.inter$baseFac
        inter <- asr.inter$inter
        varLevels <- asr.inter$varLevels
        IBV <- asr.inter$IBV
        JBV <- asr.inter$JBV
        coords <- asr.inter$coords
    }
    if ((length(cinv[[2]]) > 0) & inherits(cinv, "formula")) {
        Cfacnam <- vector(mode = "character")
        Cfacnam <- "dummyIntercept"
        CbaseFac <- vector(mode = "list", length = 3)
        names(CbaseFac) <- c("facNum", "baseObj", "baseFun")
        Cinter <- matrix(0, nrow = 9, ncol = 1)
        keep.order <- T
        asr.inter <- asreml.terms(data, cinv, Cfacnam, varLevels, 
            CbaseFac, Cinter, IBV, JBV, coords, 4, keep.order, 
            locgiv)
        who <- (asr.inter$inter[rn("inform"), ] > 0)
        Cfacnam <- asr.inter$facnam[who]
        Cflev <- asr.inter$inter[rn("flev"), ][who]
    }
    else Cfacnam <- Cflev <- NULL
    flev <- inter[rn("flev"), ]
    nlev <- inter[rn("nlev"), ]
    inform <- inter[rn("inform"), ]
    facord <- inter[rn("facord"), ]
    faconst <- inter[rn("faconst"), ]
    inter <- inter[1:4, ]
    if (is.null(IBV)) 
        IBV <- matrix(c(0, 0, 0, 0), nrow = 1)
    if (is.null(JBV)) 
        JBV <- matrix(c(0, 0, 0, 0), nrow = 1)
    xsub <- rep(FALSE, length(inform))
    fn <- unique(c(match(baseFac$baseObj, facnam), match(baseFac$baseFun, 
        facnam)))
    fn <- fn[!is.na(fn)]
    xsub[fn] <- (inform < 0)[fn]
    inform[xsub] <- abs(inform[xsub])
    x <- (inform > 0) & (facord > 0)
    In.Form <- inform[x]
    Fac.Ord <- facord[x]
    Indx <- seq(1, length(inform))[x]
    fxd <- (Indx[In.Form == 1])[order(Fac.Ord[In.Form == 1])]
    rndm <- (Indx[In.Form == 2])[order(Fac.Ord[In.Form == 2])]
    sprs <- (Indx[In.Form == 3])[order(Fac.Ord[In.Form == 3])]
    fxd <- rev(fxd[-1])
    if (sum(which <- (inform == 1) & (xsub == 1) & (inter[1, 
        ] != -12))) 
        fxd <- c(fxd, seq(1, length(xsub))[which])
    if (sum(which <- (inform == 2) & (xsub == 1) & (inter[1, 
        ] != -12))) 
        rndm <- c(rndm, seq(1, length(xsub))[which])
    if (sum(which <- (inform == 3) & (xsub == 1) & (inter[1, 
        ] != -12))) 
        sprs <- c(sprs, seq(1, length(xsub))[which])
    i <- 1
    while (i < length(fxd)) {
        if (inter[1, fxd[i]] == -12) {
            tmp <- fxd[i]
            fxd[i] <- fxd[i + 1]
            fxd[i + 1] <- tmp
            i <- i + 1
        }
        i <- i + 1
    }
    neword <- c(fxd, rndm, sprs)
    fixedFlag <- c(rep(TRUE, length(fxd)), rep(FALSE, length(rndm) + 
        length(sprs)))
    sparseFlag <- c(rep(FALSE, length(fxd) + length(rndm)), rep(TRUE, 
        length(sprs)))
    equations <- neword[inter[1, neword] != -12]
    neq <- sum(flev[equations]) + 1
    neqd <- sum(flev[equations][inform[equations] == 1]) + 1
    neqf <- neqd - 1
    neqr <- sum(flev[equations][inform[equations] == 2])
    neqs <- sum(flev[equations][inform[equations] == 3])
    afxd <- NULL
    arndm <- NULL
    asprs <- NULL
    if (length(fxd) > 0) 
        afxd <- fxd[inter[1, fxd] != -12]
    if (length(rndm) > 0) 
        arndm <- rndm[inter[1, rndm] != -12]
    if (length(sprs) > 0) 
        asprs <- sprs[inter[1, sprs] != -12]
    nf <- length(afxd)
    nr <- length(arndm)
    ns <- length(asprs)
    sub <- !xsub[equations]
    eq <- seq(2, neq)
    eqf <- NULL
    if (nf > 0) 
        eqf <- seq(2, neqf + 1)[rep(sub[1:nf], flev[afxd])]
    eqr <- NULL
    if (nr > 0) 
        eqr <- seq(neqf + 2, neqf + neqr + 1)[rep(sub[seq(nf + 
            1, nf + nr)], flev[arndm])]
    eqs <- NULL
    if (ns > 0) 
        eqs <- seq(neqf + neqr + 2, neqf + neqr + neqs + 1)[rep(sub[seq(nf + 
            nr + 1, nf + nr + ns)], flev[asprs])]
    baseFac <- data.frame(baseFac)
    which <- seq(1, length(neword))[(inter[1, ][neword] == -12)]
    if (length(which) > 0) {
        conform <- vector(mode = "list", length = length(which))
        for (i in which) {
            k <- rank(which)[which == i]
            kk <- inter[2, ][neword[i]]
            lvl0 <- flev[kk]
            conform[[k]] <- facnam[kk]
            lvl1 <- 0
            for (j in seq((i - 1), 1, -1)) {
                lvl1 <- lvl1 + flev[neword[j]]
                conform[[k]] <- c(conform[[k]], facnam[neword[j]])
                if (lvl0 < lvl1) 
                  stop(paste("Factors '", paste(conform[[k]], 
                    collapse = ","), "' referenced by and() have levels that do not conform.\n"))
                if (lvl0 == lvl1) 
                  break
            }
        }
    }
    list(intercept = terms.mu, facnam = facnam, nlev = nlev, 
        flev = flev, inter = inter, inform = inform, submod = 0, 
        xsub = as.numeric(xsub), submds = as.numeric(xsub[neword]), 
        baseFac = baseFac, varLevels = varLevels, IBV = IBV, 
        JBV = JBV, points = attr(data, "Control")$splXtras$points, 
        coords = coords, neword = neword, neq = neq, neqd = neqd, 
        neqf = neqf, neqr = neqr, neqs = neqs, eqf = eqf, eqr = eqr, 
        eqs = eqs, fixedFlag = fixedFlag, sparseFlag = sparseFlag, 
        equations = equations, faconst = faconst, Coptions = list(Cfacnam = Cfacnam, 
            Cflev = Cflev))
}
asreml.inverse.gaussian <-
function (link = "1/mu^2", dispersion = NA) 
{
    link <- as.character(substitute(link))
    misnames <- c("1/mu^2")
    corresp <- c(1)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Inverse gaussiann link is \"1/mu^2\"")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("inverse.gaussian", link = link)
    fam$dispersion <- dispersion
    fam
}
asreml.isp <-
function (x, y, p = 1, init = NA, dist = NA, data, Rcov, ...) 
{
    fun <- "isp"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (is.character(p)) 
        p <- eval(parse(text = p))
    if (!is.name(substitute(y))) {
        type <- "i"
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
            x <- substitute(x)
        out <- asreml.specialsOned(fun, x, type, data, struc, 
            dist, init, p)
    }
    else {
        type <- "si"
        if (!is.na(dist)) 
            warning("dist argument ignored")
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
        out <- do.call("asreml.specialsTwod", list(fun = fun, 
            type = type, xx = xx, yy = yy, data = data, struc = struc, 
            init = init, Rcov = Rcov, p = p))
    }
    out
}
asreml.isph <-
function (x, y, p = 1, init = NA, dist = NA, data, Rcov, ...) 
{
    fun <- "isph"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (is.character(p)) 
        p <- eval(parse(text = p))
    if (!is.name(substitute(y))) {
        type <- "h"
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
            x <- substitute(x)
        out <- asreml.specialsOned(fun, x, type, data, struc, 
            dist, init, p)
    }
    else {
        type <- "sh"
        if (!is.na(dist)) 
            warning("dist argument ignored")
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
        out <- do.call("asreml.specialsTwod", list(fun = fun, 
            type = type, xx = xx, yy = yy, data = data, struc = struc, 
            init = init, Rcov = Rcov, p = p))
    }
    out
}
asreml.ispv <-
function (x, y, p = 1, init = NA, dist = NA, data, Rcov, ...) 
{
    fun <- "ispv"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (is.character(p)) 
        p <- eval(parse(text = p))
    if (!is.name(substitute(y))) {
        type <- "v"
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) 
            x <- substitute(x)
        out <- asreml.specialsOned(fun, x, type, data, struc, 
            dist, init, p)
    }
    else {
        type <- "sv"
        if (!is.na(dist)) 
            warning("dist argument ignored")
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
        out <- do.call("asreml.specialsTwod", list(fun = fun, 
            type = type, xx = xx, yy = yy, data = data, struc = struc, 
            init = init, Rcov = Rcov, p = p))
    }
    out
}
asreml.labsoln <-
function (ntrt, neword, neq, facnam, varLevels) 
{
    asreml.generate <- function(lvls) {
        N <- length(lvls)
        x <- vector(length = (N + 2))
        x[1] <- x[N + 2] <- 1
        x[2:(N + 1)] <- sapply(lvls, function(x) {
            length(x)
        })
        indx <- matrix("", nrow = prod(x[2:(N + 1)]), ncol = N)
        for (i in seq(2, (N + 1))) indx[, i - 1] <- rep(rep(lvls[[i - 
            1]], rep(prod(x[seq(i + 1, N + 2)]), x[i])), prod(x[seq(1, 
            i - 1)]))
        indx
    }
    labels <- vector(mode = "character", length = neq)
    labels[seq(1, ntrt)] <- facnam[seq(1, ntrt)]
    ip <- ntrt + 1
    for (n in neword) {
        i <- grep(":", substring(facnam[n], seq(1, nchar(facnam[n])), 
            seq(1, nchar(facnam[n]))))
        if (length(i) == 0) 
            var <- facnam[n]
        else var <- substring(facnam[n], c(1, i + 1), c(i - 1, 
            nchar(facnam[n])))
        lvls <- lapply(var, function(z, varLevels) {
            i <- match(z, names(varLevels))
            if (is.na(i)) 
                return(z)
            else if (length(varLevels[[i]]) == 1) {
                if (z == varLevels[[i]]) 
                  return(varLevels[[i]])
                else return(paste(z, varLevels[[i]], sep = "_"))
            }
            else return(paste(z, varLevels[[i]], sep = "_"))
        }, varLevels)
        z <- asreml.generate(lvls)
        lx <- dim(z)[1]
        fun <- sapply(seq(1, ncol(z)), function(x) {
            paste("z[,", x, "]")
        })
        fun <- paste("paste(", paste(fun, collapse = ","), ",sep=':')")
        x <- eval(parse(text = fun))
        labels[seq(ip, ip + lx - 1)] <- x
        ip <- ip + lx
    }
    labels
}
asreml.ldl <-
function (gammas) 
{
    N <- length(gammas)
    which <- rep(1, N)
    which[asreml.grep("<NotEstimated>$", names(gammas))] <- 0
    n <- (sqrt(8 * N + 1) - 1)/2
    q <- ((2 * n - 1) - sqrt((2 * n + 1)^2 - 8 * length(gammas[which == 
        1])))/2
    L <- D <- matrix(0, nrow = n, ncol = n)
    xx <- 0 <= (row(L) - col(L)) & (row(L) - col(L)) <= q
    L[xx] <- gammas[which == 1]
    diag(D) <- diag(L)
    diag(L) <- 1
    V <- L %*% D %*% t(L)
    V[!lower.tri(V)]
}
asreml.levels <-
function (x, Rcov = TRUE) 
{
    if (Rcov) {
        full <- levels(x)
        present <- unique(x)
        return(full[match(present, full)])
    }
    else return(levels(x))
}
asreml.lic <-
function (license = "asreml.lic", install = FALSE) 
{
    pv <- package_version(packageVersion("asreml"))
    ver <- c(pv$major, pv$minor)
    storage.mode(ver) <- "integer"
    stat <- as.integer(0)
    inst <- as.integer(ifelse(install, 1, 0))
    where <- paste(asreml.findMe(), collapse = " ")
    env <- names(Sys.getenv())
    exe <- as.integer(0)
    if (sum(length(grep("^RSTUDIO", env)), length(grep("^EMACS", 
        env)))) 
        exe <- as.integer(1)
    lic <- .C("ckLicense", where, ver, exe, lic = license, stat = stat, 
        inst = inst)
    return(invisible())
}
asreml.lin <-
function (obj, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to lin() must be a simple object (factor or variate)\n")
    out <- asreml.spc("lin", "mdl", 0, 0, NA, substitute(obj), 
        NA, data, match.call()$Rcov)
    out[[4]] <- out[[2]]
    out[[10]] <- c(1, 0, NA, 0)
    out
}
asreml.link <-
function (form, data, ...) 
{
    stop("link() is deprecated, please see str() instead")
}
asreml.ltri2mat <-
function (vec, diag = TRUE, cor = TRUE) 
{
    makeMat <- function(v, m, diagonal) {
        mat <- matrix(nrow = m, ncol = m)
        mat[t(lower.tri(mat, diagonal))] <- v
        mat <- t(mat)
        mat[!lower.tri(mat, TRUE)] <- t(mat)[!lower.tri(mat, 
            TRUE)]
        mat
    }
    if (is.logical(diag)) {
        n <- length(vec)
        dii <- numeric(0)
        if (diag) 
            m <- (sqrt(8 * n + 1) - 1)/2
        else {
            m <- (sqrt(8 * n + 1) + 1)/2
            dii <- rep(ifelse(cor, 1, 0), m)
        }
        mat <- makeMat(vec, m, diag)
    }
    else if (is.numeric(diag)) {
        nd <- length(diag)
        n <- length(vec) - nd
        dii <- {
            if (nd == 1) 
                rep(vec[diag], n)
            else vec[diag]
        }
        m <- (sqrt(8 * n + 1) + 1)/2
        mat <- makeMat(vec[1:n], m, FALSE)
    }
    if (length(dii)) 
        diag(mat) <- dii
    return(mat)
}
asreml.ma <-
function (obj, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to ma() must be a factor\n")
    out <- asreml.spc("ma", "mdl", 0, 0, NA, substitute(obj), 
        NA, data, match.call()$Rcov)
    out[[10]] <- c(-6, 0, NA, 0)
    out
}
asreml.ma1 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma1", "cor", 0.1, numeric(0), NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ma1h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma1h", "cor", 0.1, "rep(0.1,n)", NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ma1v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma1v", "cor", 0.1, 0.1, NA, obj, init, 
        data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
        0)
    out
}
asreml.ma2 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma2", "cor", "rep(0.1,2)", numeric(0), 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ma2h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma2h", "cor", "rep(0.1,2)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.ma2v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("ma2v", "cor", "rep(0.1,2)", 0.1, NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 
        0)
    out
}
asreml.makeFamily <-
function (dist, link, phi = NA) 
{
    variance.name <- switch(dist, binomial = "mu(1-mu)", gaussian = "constant", 
        Gamma = "mu^2", inverse.gaussian = "mu^3", poisson = "mu", 
        negative.binomial = "mu+mu^2/phi", multinomial = "mu(1-mu)")
    if (is.null(variance.name)) 
        stop(paste(dist, "not implemented\n"))
    if (asreml.Rsys) {
        fam.obj <- switch(dist, negative.binomial = do.call(dist, 
            list(link = link, phi = phi)), do.call(dist, list(link = link)))
    }
    else {
        if (dist == "gaussian") 
            fam.obj <- do.call(dist)
        else fam.obj <- do.call(dist, list(link = link))
    }
    if (asreml.Rsys) {
        fam <- list()
        fam$family <- c(fam.obj$family, paste(fam.obj$link, ":", 
            sep = ""))
        fam$link <- fam.obj$linkfun
        fam$inverse <- fam.obj$linkinv
        fam$deriv <- asreml.glmDeriv(link)
        fam$initialize <- fam.obj$initialize
        fam$variance <- fam.obj$variance
        fam$deviance <- asreml.glmDeviance(variance.name)
        fam$weight <- switch(dist, binomial = parse(text = "w*mu*(1.0-mu)"), 
            gaussian = parse(text = "w"), Gamma = parse(text = "w*mu^2"), 
            inverse.gaussian = parse(text = "w/((sqrt(family$variance(mu))*family$deriv(mu))^2.)"), 
            poisson = parse(text = "w*mu"))
    }
    else {
        fam <- fam.obj
        fam$family <- casefold(fam$family)
    }
    oldClass(fam) <- "asreml.family"
    fam
}
asreml.man <-
function (browser = "acroread") 
{
    unx <- FALSE
    where <- asreml.findMe()
    if (asreml.Rsys) 
        pdf <- "asreml-R.pdf"
    else pdf <- "asreml-S.pdf"
    if (version$os == "Microsoft Windows") 
        prog <- ""
    else if (asreml.Rsys) {
        where <- paste(where, "doc/", sep = "")
        if (length(grep("linux", version$os)) == 0) {
            prog <- ""
        }
        else {
            unx <- TRUE
            prog <- browser
        }
    }
    else {
        unx <- TRUE
        prog <- browser
    }
    if (unx) 
        prog <- paste(prog, paste(where, pdf, sep = ""), "&")
    else where <- paste(where, pdf, sep = "")
    what <- .C("Shelp", prog, where)
    invisible(NULL)
}
asreml.mat2df <-
function (ginverse, tol = 1e-08) 
{
    if (!is.matrix(ginverse)) 
        stop("Argument must be a matrix")
    if (nrow(ginverse) != ncol(ginverse)) 
        stop("Argument must be a square matrix")
    ginverse[abs(ginverse) < tol] <- NA
    if (length(dimnames(ginverse)[[1]]) == 0) {
        warning("Missing row component of dimnames for ginverse matrix\n")
        rowNames <- seq(1, nrow(ginverse))
    }
    else rowNames <- dimnames(ginverse)[[1]]
    mat <- matrix((!is.na(as.vector(ginverse))), nrow = nrow(ginverse), 
        ncol = ncol(ginverse))
    flag <- as.vector(t(lower.tri(ginverse, diag = TRUE)))
    flag <- flag & t(mat)
    row <- as.vector(col(mat))[flag]
    col <- as.vector(row(mat))[flag]
    ginverse <- data.frame(cbind(row, col, as.vector(t(ginverse))[flag]))
    names(ginverse) <- c("row", "column", "value")
    attr(ginverse, "rowNames") <- rowNames
    ginverse
}
asreml.matchCon <-
function (init) 
{
    what <- names(init)
    if (is.null(what)) 
        return(NULL)
    codes <- c("P", "U", "F")
    what <- casefold(what, upper = TRUE)
    which <- charmatch(what, codes)
    where <- is.na(which)
    if (all(where)) 
        return(NULL)
    else if (sum(as.numeric(where)) > 0) 
        what[where] <- "F"
    return(what)
}
asreml.matchInteraction <-
function (str, vec) 
{
    str <- strsplit(str, ":", fixed = TRUE)[[1]]
    ss <- unlist(lapply(strsplit(vec, ":", fixed = TRUE), function(x) paste(sort(x), 
        collapse = ":")))
    match(paste(sort(str), collapse = ":"), ss, nomatch = NA)
}
asreml.mbf <-
function (obj, data, ...) 
{
    obj <- as.character(substitute(obj))
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to mbf() must be a numeric or character vector\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "mbf"
    out[[2]] <- obj
    out[[3]] <- obj
    what <- match(obj, names(attr(data, "MBF")$mbfAttr))
    which <- attr(data, "MBF")$mbfAttr[[obj]]$colNames
    key <- attr(data, "MBF")$mbfAttr[[obj]]$key[1]
    out[[4]] <- which[-1]
    out[[5]] <- NA
    out[[6]] <- ""
    out[[7]] <- ""
    out[[8]] <- 2
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(-24, what, match(key, names(data)), length(which) - 
        1)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- which[1]
    out[[13]] <- NA
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.mf.default <-
function (formula, data, subset = NULL, na.action = na.fail, 
    drop.unused.levels = FALSE, ...) 
{
    formula <- as.formula(formula)
    if (!inherits(formula, "terms")) 
        formula <- terms(formula, data = data)
    env <- environment(formula)
    rownames <- .row_names_info(data, 0L)
    vars <- attr(formula, "variables")
    predvars <- attr(formula, "predvars")
    if (is.null(predvars)) 
        predvars <- vars
    varnames <- dimnames(attr(formula, "factors"))[[1]]
    variables <- eval(predvars, data, env)
    resp <- attr(formula, "response")
    if (is.null(rownames) && resp > 0L) {
        lhs <- variables[[resp]]
        rownames <- if (is.matrix(lhs)) 
            rownames(lhs)
        else names(lhs)
    }
    if (is.null(attr(formula, "predvars"))) {
        for (i in seq_along(varnames)) predvars[[i + 1L]] <- makepredictcall(variables[[i]], 
            vars[[i + 1L]])
        attr(formula, "predvars") <- predvars
    }
    extras <- substitute(list(...))
    extranames <- names(extras[-1L])
    extras <- eval(extras, data, env)
    subset <- eval(substitute(subset), data, env)
    data <- .External2(stats:::C_modelframe, formula, rownames, 
        variables, varnames, extras, extranames, subset, na.action, 
        PACKAGE = "stats")
    if (drop.unused.levels) {
        for (nm in names(data)) {
            x <- data[[nm]]
            if (is.factor(x) && length(unique(x[!is.na(x)])) < 
                length(levels(x))) 
                data[[nm]] <- data[[nm]][, drop = TRUE]
        }
    }
    attr(formula, "dataClasses") <- vapply(data, .MFclass, "")
    attr(data, "terms") <- formula
    data
}
asreml.ModelFormula <-
function (fixed, random, sparse, rcov, weights, offset, groups, 
    mbfKeys, ignore, namesInData) 
{
    myEval <- function(obj, namesInData) {
        expr <- obj$expr
        args <- obj$args
        if (is.character(expr)) 
            expr <- parse(text = expr)
        switch(mode(expr[[1]]), call = {
            expr11 <- as.character(expr[[1]][[1]])
            expr12 <- as.character(expr[[1]][[2]])
            if (length(expr[[1]]) > 1) {
                if (!is.na(which <- match(expr11, asreml.Spcls$Fun))) {
                  if (asreml.Spcls$ObjArgs[which] == -1) N <- length(expr[[1]]) else N <- asreml.Spcls$ObjArgs[which] + 
                    1
                } else N <- length(expr[[1]])
                if (N > 1) {
                  for (i in 2:N) {
                    zz <- myEval(list(expr = expr[[1]][i], args = args), 
                      namesInData)
                    expr[[1]][i] <- zz$expr
                    args <- zz$args
                  }
                }
            }
        }, numeric = {
        }, name = {
            if (length(namesInData) && is.element(as.character(expr[[1]]), 
                namesInData)) args <- c(args, as.character(expr[[1]])) else if (!length(namesInData)) args <- c(args, 
                as.character(expr[[1]]))
        }, character = {
        }, logical = {
        }, stop("Missing function argument(s) ??"))
        list(expr = expr, args = args)
    }
    rsp <- attr(terms(fixed), "response")
    y <- dimnames(attr(terms(fixed), "factors"))[[1]][rsp]
    if (is.null(y)) 
        y <- deparse(fixed[[2]])
    if (length(c(dimnames(attr(terms(fixed), "factors"))[[1]][-rsp], 
        dimnames(attr(terms(random), "factors"))[[1]], dimnames(attr(terms(sparse), 
            "factors"))[[1]], dimnames(attr(terms(rcov), "factors"))[[1]], 
        weights, offset)) == 0) {
        form <- formula(paste(y, "~ 1"))
        attr(form, "model.y") <- y
        return(form)
    }
    if (length(groups) > 0) 
        gnames <- unique(unlist(groups))
    else gnames <- NULL
    if (length(mbfKeys) > 0) 
        knames <- unlist(lapply(mbfKeys, function(x) x$key[1]))
    else knames <- NULL
    form <- formula(paste("~", paste(unique(c({
        if (!is.null(dim(attr(terms(fixed), "factors")))) dimnames(attr(terms(fixed), 
            "factors"))[[1]][-rsp]
    }, dimnames(attr(terms(random), "factors"))[[1]], dimnames(attr(terms(sparse), 
        "factors"))[[1]], dimnames(attr(terms(rcov), "factors"))[[1]], 
        gnames, knames, weights, offset)), collapse = "+")))
    tt <- terms(form, specials = asreml.Spcls$Fun)
    vars <- dimnames(attr(tt, "factors"))[[1]]
    is.special <- rep(FALSE, length(vars))
    is.special[unlist(sapply(attr(tt, "specials"), function(x) x))] <- TRUE
    args <- character(0)
    for (i in 1:length(vars)) {
        if (is.special[i]) 
            args <- myEval(list(expr = vars[i], args = args), 
                namesInData)$args
        else args <- c(args, vars[i])
    }
    args <- unique(args)
    model.y <- character(0)
    model.y <- myEval(list(expr = fixed[2], model.y), character(0))$args
    which <- is.na(match(args, ignore))
    if (sum(as.numeric(which)) > 0) 
        args <- args[which]
    else if (sum(as.numeric(which)) == 0) 
        args <- "1"
    form <- formula(paste(y, "~", paste(args, collapse = "+")))
    attr(form, "model.y") <- model.y
    form
}
asreml.modelFrame <-
function (formula, y, data, na.method.Y, na.method.X, drop.unused.levels) 
{
    formula <- formula[!is.na(match(formula, names(data)))]
    ntrt <- length(y)
    y <- formula[1:ntrt]
    x <- formula[-(1:ntrt)]
    x <- x[x != "mv" & x != "trait"]
    if (length(y) == 1 && is.factor(data[[y]])) 
        stop("Response is a factor\n")
    subset <- rep(TRUE, nrow(data))
    where <- logical(nrow(data))
    for (i in y) where <- (where | is.na(data[, i]))
    switch(na.method.Y, fail = {
        if (any(where)) stop("Missing values in the response, na.method.Y='fail'\n")
    }, omit = {
        subset <- !where
    })
    switch(na.method.X, fail = {
        for (i in x) {
            if (is.factor(data[, i])) {
                if (any(is.na(data[, i])) || any(as.character(data[, 
                  i]) == "NA") || any(is.na(data[, i] == "NA"))) stop(paste("Missing values in explanatory factor", 
                  "(", i, "):", "na.method.X='fail'\n"))
            } else {
                if (any(is.na(data[, i]))) stop(paste("Missing values in explanatory variable", 
                  "(", i, "):", "na.method.X='fail'\n"))
            }
        }
    }, omit = {
        for (i in x) {
            if (is.factor(data[, i])) where <- (is.na(data[[i]]) | 
                (as.character(data[, i]) == "NA") | (is.na(data[, 
                i] == "NA"))) else where <- is.na(data[[i]])
            subset <- !where & subset
        }
    })
    df <- data[subset, formula]
    rn <- row.names(df)
    names(df) <- NULL
    df <- lapply(df, function(x, drop.unused.levels) {
        if (inherits(x, "factor")) {
            if (drop.unused.levels) {
                return(factor(x))
            }
            else return(x)
        }
        else return(x)
    }, drop.unused.levels)
    df <- data.frame(df, row.names = rn)
    names(df) <- formula
    df
}
asreml.ModelFunctions <-
function (path = "./") 
{
    specials <- read.table(paste(path, "asreml.ModelFunctions.csv", 
        sep = ""), header = TRUE, sep = ",", row.names = NULL, 
        as.is = TRUE)
    invisible(specials)
}
asreml.monte <-
function (pedigree, fgen = list(character(0), 0), nsim = 2000, 
    rm = c("A", "D", "C", "Ct", "Dh", "Di", "E"), form = c("full", 
        "lower"), mv = c("NA", "0", "*"), psort = FALSE) 
{
    is.missing <- function(pcol, mv) {
        isna <- is.element(as.character(pcol), mv)
        isna <- isna | is.na(pcol)
        return(isna)
    }
    rm.vec <- c("A", "D", "C", "Ct", "Dh", "Di", "E")
    which.rm <- match(rm, rm.vec)
    if (any(is.na(which.rm))) 
        stop("'rm' must be one or more of 'A','D','C','Ct','Dh','Di','E'")
    if (any(dup <- duplicated(as.character(pedigree[, 1])))) {
        stop(cat("Duplicated individuals", "\n", as.character(pedigree[, 
            1])[dup], "\n"))
    }
    which <- seq(1, nrow(pedigree))[as.character(pedigree[, 1]) == 
        as.character(pedigree[, 2]) | as.character(pedigree[, 
        1]) == as.character(pedigree[, 3])]
    if (length(which <- which[!is.na(which)])) {
        cat(paste(which, pedigree[, 1][which]), sep = "\n")
        stop("Individuals appear as own parent")
    }
    form = match.arg(form)
    fmode <- FALSE
    col4 <- -1
    col4Name <- NULL
    fg <- as.double(fgen[[2]])
    if (length(fgen[[1]]) > 0) {
        if (length(col4 <- pedigree[, fgen[[1]]]) == 0) 
            stop("Argument 'fgen' specifies non-existent dataframe column.")
        else {
            fmode <- TRUE
            col4Name <- fgen[[1]]
            w <- is.na(pedigree[, col4Name])
            if (sum(as.numeric(w)) > 0) {
                if (sum(as.numeric(!(is.missing(pedigree[w, 2], 
                  mv) | is.missing(pedigree[w, 3], mv)))) > 0) 
                  stop("fgen is 'NA' for non-base individuals")
            }
        }
    }
    lped <- nrow(pedigree)
    charPed <- character(lped * 3)
    for (i in 1:3) charPed[seq((i - 1) * lped + 1, i * lped)] <- as.character(pedigree[, 
        i])
    charPed[is.na(charPed)] <- "NA"
    lvls <- unique(charPed)
    if (any(!is.na(which <- match(mv, lvls)))) {
        charPed[!is.na(match(charPed, mv))] <- "0"
        lvls <- lvls[-which[!is.na(which)]]
    }
    nan <- length(lvls)
    newped <- rep(paste(rep(" ", (max(nchar(lvls)) + 1)), collapse = ""), 
        3 * nan)
    if (nan > lped) {
        charPed <- c(charPed, rep(paste(rep(" ", (max(nchar(lvls)) + 
            1)), collapse = ""), 3 * (nan - lped)))
        if (fmode) 
            col4 <- c(col4, rep(0, (nan - lped)))
    }
    fgsx <- rep(0, nan)
    error <- paste(rep(" ", 256), collapse = "")
    storage.mode(nan) <- "integer"
    storage.mode(lped) <- "integer"
    storage.mode(fg) <- "double"
    storage.mode(col4) <- "double"
    storage.mode(fgsx) <- "double"
    out <- .C("pdsort", lped, nan, lvls, charPed, col4, fg, newped = newped, 
        fgsx = fgsx, error = error, NAOK = TRUE)
    if (nchar(out$error) > 0) {
        stop(out$error)
    }
    pds <- {
        if (length(col4Name)) 
            cbind(as.data.frame(matrix(out$newped, nrow = nan, 
                byrow = FALSE)), out$fgsx)
        else as.data.frame(matrix(out$newped, nrow = nan, byrow = FALSE))
    }
    names(pds) <- c(names(pedigree)[1:3], col4Name)
    if (psort) 
        return(pds)
    mped <- {
        if (length(col4Name)) 
            cbind(as.data.frame(matrix(match(out$newped, lvls, 
                nomatch = 0), nrow = nan, byrow = FALSE)), out$fgsx)
        else as.data.frame(matrix(match(out$newped, lvls, nomatch = 0), 
            nrow = nan, byrow = FALSE))
    }
    mped <- lapply(mped, function(x) {
        storage.mode(x) <- "double"
        x
    })
    storage.mode(nsim) <- "integer"
    x <- .Call("monte", mped, nsim)
    out <- vector(mode = "list", length = length(rm))
    names(out) <- rm
    if (form == "full") {
        for (i in which.rm) {
            out[[rm.vec[i]]] <- asreml.utri2mat(x$RM[, i])
            dimnames(out[[rm.vec[i]]]) <- list(lvls, lvls)
        }
    }
    else {
        for (i in which.rm) {
            out[[rm.vec[i]]] <- x$RM[, i]
            attr(out[[rm.vec[i]]], "rowNames") <- lvls
        }
    }
    out$Fi <- x$Fi
    out
}
asreml.mtrn <-
function (x, y, phi = NA, nu = 0.5, delta = 1, alpha = 0, lambda = 2, 
    data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "mtrn"
    type <- "mi"
    struc <- c(0, 0, 6, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    attr(phi, "missing") <- asreml.ie(missing(phi), TRUE, FALSE)
    attr(nu, "missing") <- asreml.ie(missing(nu), TRUE, FALSE)
    attr(delta, "missing") <- asreml.ie(missing(delta), TRUE, 
        FALSE)
    attr(alpha, "missing") <- asreml.ie(missing(alpha), TRUE, 
        FALSE)
    attr(lambda, "missing") <- asreml.ie(missing(lambda), TRUE, 
        FALSE)
    init <- NA
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, phi = phi, nu = nu, delta = delta, alpha = alpha, 
        lambda = lambda, dist = dist))
    out
}
asreml.mtrnh <-
function (x, y, phi = NA, nu = 0.5, delta = 1, alpha = 0, lambda = 2, 
    init = 0.1, data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "mtrnh"
    type <- "mh"
    struc <- c(0, 0, 6, 0, 0, 9, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0)
    attr(phi, "missing") <- asreml.ie(missing(phi), TRUE, FALSE)
    attr(nu, "missing") <- asreml.ie(missing(nu), TRUE, FALSE)
    attr(delta, "missing") <- asreml.ie(missing(delta), TRUE, 
        FALSE)
    attr(alpha, "missing") <- asreml.ie(missing(alpha), TRUE, 
        FALSE)
    attr(lambda, "missing") <- asreml.ie(missing(lambda), TRUE, 
        FALSE)
    attr(init, "missing") <- asreml.ie(missing(init), TRUE, FALSE)
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        if (!all(is.na(dist))) 
            warning("dist argument ignored")
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, phi = phi, nu = nu, delta = delta, alpha = alpha, 
        lambda = lambda, dist = dist))
    out
}
asreml.mtrnv <-
function (x, y, phi = NA, nu = 0.5, delta = 1, alpha = 0, lambda = 2, 
    init = 0.1, data, Rcov, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$dist)) 
        dist <- eval(call$dist)
    else dist <- list()
    fun <- "mtrnv"
    type <- "mv"
    struc <- c(0, 0, 6, 0, 0, 9, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0)
    attr(phi, "missing") <- asreml.ie(missing(phi), TRUE, FALSE)
    attr(nu, "missing") <- asreml.ie(missing(nu), TRUE, FALSE)
    attr(delta, "missing") <- asreml.ie(missing(delta), TRUE, 
        FALSE)
    attr(alpha, "missing") <- asreml.ie(missing(alpha), TRUE, 
        FALSE)
    attr(lambda, "missing") <- asreml.ie(missing(lambda), TRUE, 
        FALSE)
    attr(init, "missing") <- asreml.ie(missing(init), TRUE, FALSE)
    if (!is.name(substitute(y))) {
        if (!is.numeric(y)) 
            stop("y must be a numeric scalar or vector")
        y[is.na(y)] <- 0
        yy <- y
        if (!(mode(substitute(x)) == "call" && inherits(x, "asreml.special"))) {
            xx <- as.character(substitute(x))
        }
        else {
            xx <- x
        }
    }
    else {
        if (!all(is.na(dist))) 
            warning("dist argument ignored")
        xx <- as.character(substitute(x))
        yy <- as.character(substitute(y))
    }
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov, phi = phi, nu = nu, delta = delta, alpha = alpha, 
        lambda = lambda, dist = dist))
    out
}
asreml.multinomial <-
function (link = "logit", dispersion = 1) 
{
    link <- as.character(substitute(link))
    misnames <- c("logit", "probit", "cloglog", "Logit", "Probit", 
        "clog-log", "Cloglog", "Clog-log")
    corresp <- c(1, 2, 3, 1, 2, 3, 3, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Multinomial links are \"logit\", \"probit\", or \"cloglog\"\n")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("multinomial", link = link)
    fam$dispersion <- dispersion
    fam
}
asreml.mvfact <-
function (y) 
{
    if (is.matrix(y)) 
        stop("Missing values not allowed\n")
    out <- rep(NA, length(y))
    out[is.na(y)] <- seq(1, sum(is.na(y)))
    factor(out)
}
asreml.negative.binomial <-
function (link = "log", dispersion = 1, phi = 1) 
{
    link <- as.character(substitute(link))
    misnames <- c("log", "identity", "inverse", "reciprocal", 
        "Log", "Identity", "Inverse", "Reciprocal")
    corresp <- c(1, 2, 3, 3, 1, 2, 3, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Negative binomial links are \"log\", \"identity\", or \"inverse\"")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("negative.binomial", link = link, 
        phi = phi)
    fam$dispersion <- dispersion
    fam$phi <- phi
    fam
}
asreml.niceGammas <-
function (obj, use.names = TRUE) 
{
    asreml.vList <- function(GR, use.names) {
        if (is.null(GR)) 
            return(NULL)
        out <- lapply(GR, function(z, nana) {
            asreml.unlist(lapply(z, function(x) {
                if (is.null(x$model) && length(x$s2)) 
                  return({
                    if (x$s2 != 1) x$s2 else NULL
                  })
                switch(x$model, id = NULL, ped = {
                  if (all(is.na(x$initial))) y <- NULL else y <- x$initial
                  y
                }, ide = {
                  if (all(is.na(x$initial))) y <- NULL else y <- x$initial
                  y
                }, giv = {
                  if (all(is.na(x$initial))) y <- NULL else y <- x$initial
                  y
                }, us = {
                  y <- asreml.ltri2mat(x$initial)
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, ante = {
                  y <- asreml.ltri2mat(x$initial)
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, chol = {
                  y <- asreml.ltri2mat(x$initial)
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, corv = {
                  n <- length(x$levels) * (length(x$levels) - 
                    1)/2
                  vec <- c(rep(x$initial[1], n), rep(x$initial[2], 
                    length(x$levels)))
                  y <- asreml.ltri2mat(vec, seq(n + 1, length(vec)))
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, corh = {
                  n <- length(x$levels) * (length(x$levels) - 
                    1)/2
                  vec <- c(rep(x$initial[1], n), x$initial[-1])
                  y <- asreml.ltri2mat(vec, seq(n + 1, length(vec)))
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, corg = {
                  y <- asreml.ltri2mat(x$initial, diag = FALSE)
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, corgv = {
                  y <- asreml.ltri2mat(x$initial, diag = length(x$init))
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, corgh = {
                  y <- asreml.ltri2mat(x$initial, seq(length(x$initial) - 
                    length(x$levels) + 1, length(x$initial)))
                  dimnames(y) <- list(x$levels, x$levels)
                  y
                }, x$initial)
            }), nana)
        }, use.names)
        which <- unlist(lapply(out, function(x) is.null(x)))
        out[!which]
    }
    asreml.unlist <- function(LL, use.names = TRUE) {
        back <- list()
        nana <- character(0)
        k <- 0
        for (i in 1:length(LL)) {
            if (is.null(LL[[i]])) 
                next
            k <- k + 1
            if (is.vector(LL[[i]])) {
                if (length(LL[[i]]) == 1) 
                  names(LL[[i]]) <- NULL
                else names(LL[[i]]) <- {
                  if (use.names) 
                    unlist(lapply(strsplit(names(LL[[i]]), "!"), 
                      function(x) {
                        x[length(x)]
                      }))
                  else NULL
                }
                back[[k]] <- LL[[i]]
                nana <- c(nana, names(LL)[i])
            }
            else if (is.matrix(LL[[i]])) {
                back[[k]] <- LL[[i]]
                dimnames(back[[k]]) <- {
                  if (use.names) 
                    dimnames(LL[[i]])
                  else list(NULL, NULL)
                }
            }
        }
        return({
            if (length(back) == 0) NULL else if (length(back) == 
                1) back[[1]] else {
                names(back) <- nana
                back
            }
        })
    }
    out <- c(asreml.vList(obj$G.param, use.names), asreml.vList(obj$R.param, 
        use.names))
    out
}
asreml.ped <-
function (obj, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to ped() must be a factor)\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "ped"
    obj <- as.character(substitute(obj))
    out[[2]] <- asreml.spCall(sys.call())
    out[[3]] <- obj
    v <- match(obj, names(data))
    if (!is.factor(data[[obj]])) 
        stop("Argument to ped() must be a factor)\n")
    gv <- attr(data, "Control")$ginverse
    if (is.null(gv$rownames)) 
        stop("Deprecated 'ginverse' component of the fitted object, rerun with a newer version of asreml()")
    glvls <- gv$rownames(obj, gv)
    if (length(flvls <- levels(data[[obj]])) == 0) 
        flvls <- unique(data[[obj]])
    if (any(is.na(match(flvls, glvls)))) {
        print(flvls[is.na(match(flvls, glvls))])
        stop(paste(obj, "has levels in data missing in ginverse\n"))
    }
    out[[4]] <- glvls
    IniFlag <- FALSE
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- 0.1
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != 1) 
            stop("ped(): Wrong number of initial values\n")
    }
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) 
        out[[6]] <- "P"
    else out[[6]] <- con
    out[[7]] <- "!ped"
    out[[8]] <- 2
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, -7, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0)
    names(out[[9]]) <- c("", "", obj, rep("", 13))
    out[[10]] <- c(-1, -1, v, NA)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- obj
    out[[13]] <- ifelse(IniFlag, TRUE, FALSE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.poisson <-
function (link = "log", dispersion = 1) 
{
    link <- as.character(substitute(link))
    misnames <- c("log", "identity", "sqrt", "Log", "Identity", 
        "Sqrt")
    corresp <- c(1, 2, 3, 1, 2, 3)
    lmatch <- pmatch(link, misnames, FALSE)
    if (!lmatch) 
        stop("Poisson links are \"log\", \"identity\" or \"sqrt\"")
    link <- misnames[corresp[lmatch]]
    fam <- asreml.makeFamily("poisson", link = link)
    fam$dispersion <- dispersion
    fam
}
asreml.pol <-
function (obj, t = 1, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to pol() must be a simple object (variate or factor)\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "pol"
    out[[2]] <- asreml.spCall(sys.call())
    obj <- as.character(substitute(obj))
    out[[3]] <- obj
    splXtras <- attr(data, "Control")$splXtras
    knots <- splXtras$knots
    inter <- c(-7, 0, 0, t)
    ibv <- inter
    jbv <- inter
    if (!is.null(splXtras$IBV)) {
        which <- match(obj, dimnames(splXtras$IBV)[[2]])
        if (!is.na(which)) {
            ibv <- splXtras$IBV[, which]
            ibv[1] <- 1
            ibv[3] <- 1
        }
    }
    nux <- length(unique(data[[obj]]))
    ncz <- min(abs(t), nux) + 1
    N <- min(splXtras$splstp[7], nux) + ncz + ibv[4] + splXtras$nsppoints
    lx <- length(data[[obj]])
    one <- 1
    knotpoints <- vector(mode = "double", length = N)
    fail <- 0
    temp <- .C("splinek", as.double(data[[obj]]), as.integer(one), 
        as.integer(lx), as.integer(inter), as.integer(ibv), as.integer(jbv), 
        as.double(splXtras$splstp), as.integer(knots), as.integer(splXtras$nsppoints), 
        as.double(splXtras$splinescale), as.double(splXtras$points), 
        as.integer(N), ncz = as.integer(ncz), nux = as.integer(nux), 
        knotpoints = as.double(knotpoints), fail = as.integer(fail), 
        NAOK = TRUE)
    if (temp$fail > 0) 
        stop()
    lvls <- temp$ncz
    if (t <= 0) 
        out[[4]] <- paste("order", seq(1, lvls), sep = "")
    else out[[4]] <- paste("order", seq(0, lvls - 1), sep = "")
    n <- length(out[[4]])
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- 0.1
    }
    else if (length(init) != 1) 
        stop("Wrong number of initial values\n")
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) 
        out[[6]] <- "P"
    else out[[6]] <- con
    out[[7]] <- paste("!", out[[2]], sep = "")
    out[[8]] <- 2
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(-7, 0, NA, t)
    out[[11]] <- list(faconst = NULL, nuxPoints = temp$knotpoints[seq(1, 
        temp$nux)], coords = NULL)
    out[[12]] <- obj
    out[[13]] <- TRUE
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.pow <-
function (obj, k = 1, offset = 0, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to pow() must be a simple object\n")
    out <- asreml.spc("pow", "mdl", 0, 0, NA, substitute(obj), 
        NA, data, match.call()$Rcov)
    out[[4]] <- out[[2]]
    out[[10]] <- c(-30, 0, match(out$Obj, names(data)), k * 1000)
    out$Coords$faconst <- offset
    out
}
asreml.powerFac <-
function (formula, data, Rcov) 
{
    if (length(formula[[2]]) == 0) 
        return(data)
    pow <- asreml.Spcls$Fun[asreml.Spcls$Type == 2]
    terms.var <- dimnames(attr(terms(formula, specials = asreml.Spcls$Fun), 
        "factors"))[[1]]
    for (v in terms.var) {
        X <- eval(asreml.Eval(v, data, Rcov))
        if (inherits(X, "asreml.special") && !is.null(X$Argv[2])) {
            if (any(!is.na(match(pow, X$Fun)))) {
                newfac <- paste(as.character(data[[X$Argv[1]]]), 
                  as.character(data[[X$Argv[2]]]), sep = ",")
                newfac[newfac == "NA,NA"] <- "NA"
                data[[X$Obj]] <- factor(newfac, levels = unique(paste(as.character(data[[X$Argv[1]]]), 
                  as.character(data[[X$Argv[2]]]), sep = ",")))
            }
        }
    }
    data
}
asreml.prdList <-
function (predict, prdStruc, inter, prdvals, vpvmat, avsed, pvtext) 
{
    typeConvert <- function(x) {
        if (is.numeric(x) || is.factor(x)) 
            return(x)
        if (asreml.Rsys) 
            return(type.convert(x, as.is = TRUE))
        else {
            ok = 1
            storage.mode(ok) <- "integer"
            n <- length(x)
            storage.mode(n) <- "integer"
            result <- .C("canBeNumeric", x = x, n, ok = ok)
            if (result$ok == 1) 
                return(as.numeric(x))
            else return(as.character(x))
        }
    }
    if (is.null(predict)) 
        return(NULL)
    prdLabels <- prdStruc$prdLabels
    npred <- 1
    prd <- vector(mode = "list", length = npred)
    estim <- c("NA", "Estimable", "Not Estimable", "Aliased")
    NPcol <- ncol(prdvals)
    baseNames <- dimnames(prdStruc$baseFac)[[2]]
    vpvBase <- 0
    Pnames <- names(prdLabels$var)
    ij <- match(Pnames, baseNames) + prdStruc$nfinmd
    swap <- order(rank(prdStruc$npfact[8, ij]))
    Pnames <- Pnames[swap]
    PA <- rep(FALSE, length(Pnames))
    if (!is.null(parList <- prdLabels$parallel)) 
        PA <- !is.na(match(Pnames, parList))
    names(PA) <- Pnames
    Pgrid <- vector(length = length(Pnames[!PA]), mode = "list")
    names(Pgrid) <- Pnames[!PA]
    Ppar <- vector(length = length(Pnames[PA]), mode = "list")
    names(Ppar) <- Pnames[PA]
    sed <- prdStruc$nptabl[6]
    margin <- prdStruc$nptabl[5] == 1
    for (Pfac in Pnames) {
        if (!PA[Pfac]) 
            Pgrid[[Pfac]] <- typeConvert(dimnames(prdLabels$var[[Pfac]])[[1]])
        else Ppar[[Pfac]] <- typeConvert(dimnames(prdLabels$var[[Pfac]])[[1]])
    }
    prd <- list()
    if (margin) {
        if (any(PA)) 
            stop("\nAttempt to fit margins with parallel option\n")
        out <- NULL
        for (mf in rev(names(Pgrid))) {
            temp <- vector(mode = "list", length = length(Pgrid))
            names(temp) <- names(Pgrid)
            for (tt in 1:length(temp)) temp[[tt]] <- NA
            temp[[mf]] <- Pgrid[[mf]]
            out <- rbind(out, as.matrix(expand.grid(temp)))
        }
        df <- expand.grid(Pgrid)
        nn <- names(df)
        df <- data.frame(df[do.call("order", df), ], row.names = NULL)
        names(df) <- nn
        prd$pvals <- data.frame(rbind(as.data.frame(out), df), 
            row.names = NULL)
    }
    else prd$pvals <- asreml.expandGrid(Pgrid, Ppar)
    loc <- prdStruc$ptrPrdvals[1]
    len <- prdStruc$ptrPrdvals[2]
    prd$pvals$predicted.value <- prdvals[seq(loc, loc + len - 
        1), 1]
    prd$pvals$standard.error <- prdvals[seq(loc, loc + len - 
        1), 2]
    prd$pvals$est.status <- estim[prdvals[seq(loc, loc + len - 
        1), 3] + 1]
    if (NPcol > 3) {
        prd$pvals$transformed.value <- prdvals[seq(loc, loc + 
            len - 1), 4]
        prd$pvals$approx.se <- prdvals[seq(loc, loc + len - 1), 
            5]
    }
    attr(prd$pvals, "heading") <- pvtext
    oldClass(prd$pvals) <- c("asremlPredict", "data.frame")
    Npv <- c(0, 0)
    if (margin) 
        Npv <- prdStruc$ptrPrdvals[2] - sum(sapply(Pgrid, function(x) {
            length(x)
        }))
    else Npv <- prdStruc$ptrPrdvals[2]
    Nvpv <- 0
    if (sed == 1 | sed == 3) {
        diag <- 1
        Nvpv <- Npv * (Npv + 1)/2
        prd$vcov <- asreml.ltri2mat(vpvmat[seq(vpvBase + 1, vpvBase + 
            Nvpv)], as.logical(diag))
        vpvBase <- vpvBase + Nvpv
    }
    if (sed == 2 | sed == 3) {
        diag <- 0
        Nvpv <- Npv * (Npv - 1)/2
        prd$sed <- asreml.ltri2mat(vpvmat[seq(vpvBase + 1, vpvBase + 
            Nvpv)], as.logical(diag), cor = FALSE)
        vpvBase <- vpvBase + Nvpv
    }
    if (margin) {
        prd$avsed <- matrix(avsed[2:10], ncol = 3, byrow = TRUE)
        dimnames(prd$avsed) <- list(c(paste("Same level of", 
            names(Pgrid)[1]), paste("Same level of", names(Pgrid)[2]), 
            paste("Different levels of", names(Pgrid)[1], "and", 
                names(Pgrid[2]))), c("min", "mean", "max"))
    }
    else if (sed == 0) 
        prd$avsed <- avsed[1]
    else if (sed > 0) {
        prd$avsed <- avsed[2:4]
        names(prd$avsed) <- c("min", "mean", "max")
    }
    return(prd)
}
asreml.prdStruc <-
function (predict, dataFrame, inter, struc, glm) 
{
    if (is.null(predict)) {
        npred <- 0
        nfinmd <- 0
        baseFac <- 0
        npfact <- 0
        nptabl <- 0
        ivals <- 0
        pvals <- 0
        lpwts <- 0
        lenPV <- 0
        lenPrdvals <- 1
        ptrPrdvals <- matrix(c(0, 1), nrow = 1, ncol = 2)
        lenVpvals <- 1
        ptrVpvals <- matrix(c(0, 1, 0, 1), nrow = 1, ncol = 4)
        prdLabels <- list()
        return(list(npred = npred, prdLabels = prdLabels, nfinmd = nfinmd, 
            baseFac = baseFac, npfact = npfact, nptabl = nptabl, 
            ivals = ivals, pvals = pvals, lpwts = lpwts, lenPrdvals = lenPrdvals, 
            ptrPrdvals = ptrPrdvals, lenVpvals = lenVpvals, ptrVpvals = ptrVpvals))
    }
    npred <- 1
    bf <- inter$baseFac$facNum
    names(bf) <- as.character(inter$baseFac$baseObj)
    axb <- asreml.grep(":", inter$facnam)
    for (ff in axb) {
        if (asreml.Rsys) {
            fax <- unlist((strsplit(inter$facnam[ff], split = ":")))
        }
        else fax <- unlist((unpaste(inter$facnam[ff], sep = ":")))
        i <- seq(1, length(bf))
        k <- match(fax, names(bf))
        if (any(is.na(k))) 
            break
        i[i[sort(k)]] <- k
        bf <- bf[i]
    }
    bFac <- bf
    bFun <- bf
    names(bFun) <- as.character(inter$baseFac$baseFun[match(bf, 
        inter$baseFac$facNum)])
    which <- (names(bFac) == names(bFun))
    baseOrd <- data.frame(bf = bf, bFac = names(bFac), bFun = names(bFun), 
        row.names = NULL)
    nfinmd <- length(unique(baseOrd$bFac))
    baseFac <- matrix(0, nrow = 8, ncol = nfinmd)
    colNames <- as.character(unique(baseOrd$bFac))
    buff <- vector(length = 8)
    buff16 <- vector(length = 16)
    for (i in seq(1, length(baseOrd$bf))) {
        ij <- baseOrd$bf[i]
        if (inter$inter[1, ij] == -8) {
            buff[1] <- ij
            buff[3] <- inter$flev[ij]
            buff[c(2, 4, 5, 6, 7, 8)] <- rep(0, 6)
        }
        else {
            addr <- inter$inter[3, ij]
            buff[c(3, 4, 5, 6, 7, 8)] <- rep(0, 6)
            buff[1] <- ij
            buff[2] <- addr
            if (inter$nlev[ij] > 1) {
                if (abs(inter$inform[ij]) == 2) 
                  buff[5] <- -inter$nlev[ij]
                else buff[4] <- -inter$nlev[ij]
            }
            else if (inter$inter[1, ij] == -3 | inter$inter[1, 
                ij] == -7) 
                buff[6] <- ij
            else if (inter$inter[1, ij] == -5) 
                buff[5] <- ij
            else if (inter$flev[ij] == 1 & inter$inter[3, ij] > 
                0 & inter$inter[1, ij] != -15) 
                buff[3] <- ij
            else if (abs(inter$inform[ij]) == 2) 
                buff[5] <- ij
            else buff[4] <- ij
            if (buff[3] == 0 && buff[6] == 0 && buff[5] == 0 && 
                buff[4] == 0) {
                for (ff in inter$facnam) {
                  if (length(asreml.grep(":", ff)) > 0) {
                    if (!is.na(match(baseOrd$bFac[i], unlist(asreml.unpaste(ff, 
                      sep = ":")))) | !is.na(match(baseOrd$bFun[i], 
                      unlist(asreml.unpaste(ff, sep = ":"))))) {
                      where <- match(ff, inter$facnam)
                      if (abs(inter$inform[where]) != 2) 
                        buff[4] <- ij
                      if (abs(inter$inform[where]) == 2) 
                        buff[5] <- ij
                    }
                  }
                }
            }
        }
        colNo <- match(baseOrd$bFac[i], colNames)
        temp <- baseFac[, colNo]
        temp[temp == 0] <- buff[temp == 0]
        baseFac[, colNo] <- temp
    }
    dimnames(baseFac) <- list(NULL, colNames)
    nptabl <- matrix(0, nrow = 16, ncol = npred)
    npfact <- matrix(0, nrow = 8, ncol = nfinmd * 2)
    pvals <- vector(mode = "double")
    prdLabels <- list()
    howMany <- 0
    ivals <- 0
    lpvals <- 1
    livals <- 1
    lenPrdvals <- 0
    numnest <- 0
    ptrPrdvals <- c(0, 0)
    vPtr <- 1
    lenVpvals <- 0
    ptrVpvals <- c(0, 0, 0, 0)
    names(ptrVpvals) <- c("locVar", "lenVar", "locSed", "lenSed")
    mvFlag <- FALSE
    Pf <- predict$classify
    Pnames <- unlist((asreml.unpaste(Pf, sep = ":")))
    buff16[16] <- nfinmd
    buff16[1:15] <- rep(0, 15)
    buff16[1] <- -1
    if (length(predict$present) > 0) 
        buff16[2] <- 0
    if (length(predict$estimable) > 0) 
        buff16[2] <- 1
    buff16[3] <- 0
    buff16[7] <- 0
    buff16[11] <- glm$link
    if (buff16[11] == 1) 
        buff16[11] <- 0
    buff16[12] <- 0
    buff16[15] <- 0
    buff16[2] <- as.numeric(predict$estimable)
    if (buff16[2] == 0) {
        a <- predict$aliased
        buff16[2] <- asreml.ie(is.logical(a), -2 * as.numeric(a), 
            a)
    }
    buff16[5] <- as.numeric(predict$margin)
    prdLabels$margin <- as.numeric(predict$margin)
    buff16[6] <- 2 * as.numeric(predict$sed) + as.numeric(predict$vcov)
    buff16[13] <- as.numeric(predict$inrandom)
    if (length(predict$exrandom) > 0 && predict$exrandom == TRUE) 
        buff16[13] <- -1
    parList <- NULL
    if (predict$parallel) {
        if (buff16[5] > 0) 
            stop("\nCannot mix margin and parallel\n")
        parList <- Pnames
        px <- match(parList, dimnames(baseFac)[[2]])
        if (any(is.na(match(parList, Pnames)))) 
            stop("\nTerms in parallel list not in classify set\n")
        if (any(is.na(px))) 
            stop("\nTerm(s) in parallel list not in model\n")
        if (length(px) < 2) 
            stop("\nParallel list must be at least length 2\n")
        prdLabels$parallel <- parList
        buff16[4] <- px[1]
        flag <- 1
        for (ipx in px) {
            npfact[7, nfinmd + ipx] <- flag
            flag <- -1
        }
        if (length(predict$levels) && any(is.na(match(parList, 
            names(predict$levels))))) 
            stop("\nMust specify levels for all parallel terms\n")
    }
    tmpIvals <- NULL
    igList <- predict$ignore
    if (!is.na(match("mv", baseOrd$bFac))) 
        igList <- c(igList, "mv")
    ufac <- asreml.grep("units", inter$facnam)
    if (length(ufac) > 0) 
        igList <- c(igList, inter$facnam[ufac])
    if (length(igList) > 0) {
        if (any(is.na(match(igList, inter$facnam)))) 
            stop("\nTerm in ignore list not in model\n")
        buff16[14] <- livals
        tmpIvals <- match(inter$facnam[inter$neword], igList, 
            nomatch = -1)
        tmpIvals[tmpIvals < 0] <- 0
        tmpIvals[tmpIvals > 0] <- -1
    }
    if (length(uList <- predict$use) > 0) {
        if (any(is.na(match(uList, inter$facnam)))) 
            stop("\nTerm in use list not in model\n")
        buff16[14] <- livals
        tmpIvals <- match(inter$facnam[inter$neword], uList, 
            nomatch = 0)
        tmpIvals[tmpIvals > 0] <- 1
    }
    if (length(eList <- predict$except) > 0) {
        if (any(is.na(match(eList, inter$facnam)))) 
            stop("\nTerm in except list not in model\n")
        buff16[14] <- livals
        tmpIvals <- match(inter$facnam[inter$neword], eList, 
            nomatch = NA)
        tmpIvals[!is.na(tmpIvals)] <- -1
        tmpIvals[is.na(tmpIvals)] <- 1
    }
    if (length(oList <- predict$only) > 0) {
        if (any(is.na(match(oList, inter$facnam)))) 
            stop("\nTerm in only list not in model\n")
        buff16[14] <- livals
        tmpIvals <- match(inter$facnam[inter$neword], oList, 
            nomatch = -1)
        tmpIvals[tmpIvals > 0] <- 1
    }
    if (!is.null(tmpIvals)) {
        ivals <- c(ivals, tmpIvals)
        livals <- length(ivals)
    }
    if (length(nList <- predict$associate) > 0) {
        buff16[9] <- livals + 1
        if (length(nList[[1]]) == 0) 
            stop("Empty associate list\n")
        tmpNvals <- match(nList[[1]], dimnames(baseFac)[[2]])
        if (any(is.na(tmpNvals))) 
            stop(paste("\n", nList[[1]][is.na(tmpNvals)], "is not in the model\n"))
        numnest <- numnest * 10 + length(tmpNvals)
        ivals <- c(ivals, numnest, tmpNvals)
        if (length(nList[[2]]) > 0) {
            tmpNvals <- match(nList[[2]], dimnames(baseFac)[[2]])
            if (any(is.na(tmpNvals))) 
                stop(paste("\n", nList[[2]][is.na(tmpNvals)], 
                  "is not in the model\n"))
            numnest <- numnest * 10 + length(tmpNvals)
            ivals[livals + 1] <- numnest
            ivals <- c(ivals, tmpNvals)
        }
        livals <- length(ivals)
    }
    nptabl <- buff16
    buff <- rep(0, 8)
    parFlag <- 1
    if (buff16[4]) 
        parFlag <- -1
    prOrd <- 0
    for (fac in Pnames) {
        if (is.na(ji <- match(fac, baseOrd$bFac))) {
            if (is.na(ji <- match(fac, baseOrd$bFun))) 
                stop(paste("Predict: ", fac, "not in model\n"))
            else Pfac <- as.character(baseOrd$bFac[ji])
        }
        else Pfac <- fac
        if (is.na(match(baseOrd$bFac[ji], names(inter$varLevels)))) {
            if (is.na(match(baseOrd$bFun[ji], names(inter$varLevels)))) 
                stop(paste("Predict: Cannot get levels for", 
                  Pfac))
            else PVfac <- as.character(baseOrd$bFun[ji])
        }
        else PVfac <- as.character(baseOrd$bFac[ji])
        if (length(predict$levels) > 0) {
            nn <- names(predict$levels)
            if (!is.na(match(fac, nn))) {
                nn[nn == fac] <- Pfac
                names(predict$levels) <- nn
            }
        }
        Pp <- inter$facnam[baseOrd$bf[ji]]
        ij <- match(Pp, inter$facnam)
        buff[1] <- ij
        npcol <- match(Pfac, dimnames(baseFac)[[2]])
        if (is.na(npcol)) 
            stop(paste(sQuote(Pfac), "not found in base factor table."))
        spline <- baseFac[6, npcol]
        noskip <- TRUE
        if (spline > 0) 
            noskip <- !is.element(inter$inter[1, spline], c(-3, 
                -7))
        pp <- asreml.getPpoints(dataFrame, Pfac)
        if (length(pp) > 0) {
            is2d <- asreml.getPdim(dataFrame, Pfac)
            if (is.null(predict$levels[[Pfac]]) && noskip) {
                predict$levels[[Pfac]] <- pp
                if (is2d > 1) {
                  parList <- c(parList, Pfac)
                  prdLabels$parallel <- c(prdLabels$parallel, 
                    Pfac)
                  npfact[7, nfinmd + npcol] <- parFlag
                  parFlag <- -1
                }
            }
        }
        fixfac <- baseFac[4, npcol]
        ranfac <- baseFac[5, npcol]
        if (length(predict$levels[[Pfac]]) > 0) {
            if (ranfac < 0 || fixfac < 0) {
                which <- inter$varLevels[[inter$facnam[ij]]]
                Plevs <- predict$levels[[Pfac]]
                if (is.matrix(Plevs)) {
                  if (ncol(Plevs) != length(which)) 
                    stop(paste("Matrix of values for", Pfac, 
                      "must have", length(which), "columns"))
                  prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                    nrow = 2)
                  prdLabels$var[[Pfac]] <- Plevs
                  if (is.null(dimnames(Plevs)[[1]])) 
                    dimnames(prdLabels$var[[Pfac]]) <- list(seq(1, 
                      nrow(Plevs)), which)
                  Plevs <- as.vector(t(Plevs))
                }
                else {
                  if (length(Plevs)%%length(which) != 0) 
                    stop(paste("Number of values for", Pfac, 
                      "must be a multiple of", length(which)))
                  nr <- length(Plevs)/length(which)
                  prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                    nrow = 2)
                  prdLabels$var[[Pfac]] <- matrix(Plevs, nrow = nr, 
                    byrow = TRUE)
                  dimnames(prdLabels$var[[Pfac]]) <- list(seq(1, 
                    nr), which)
                  Plevs <- as.vector(Plevs)
                }
            }
            else if (is.numeric(predict$levels[[Pfac]])) {
                Plevs <- predict$levels[[Pfac]]
                if (inter$flev[match(PVfac, inter$facnam)] > 
                  1) 
                  Plabs <- inter$varLevels[[PVfac]][Plevs]
                else Plabs <- Plevs
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(Plevs, ncol = 1)
                dimnames(prdLabels$var[[Pfac]]) <- list(Plabs, 
                  Pfac)
            }
            else {
                Plevs <- match(predict$levels[[Pfac]], inter$varLevels[[PVfac]])
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(Plevs, ncol = 1)
                dimnames(prdLabels$var[[Pfac]]) <- list(predict$levels[[Pfac]], 
                  Pfac)
            }
            pvals <- c(pvals, Plevs)
        }
        else {
            if (spline > 0) {
                if (inter$inter[1, spline] == -3 | inter$inter[1, 
                  spline] == -7) {
                  if (length(asreml.getPpoints(dataFrame, Pfac)) == 
                    0) {
                    pv <- inter$coords$nuxPoints[[Pfac]]
                    pvals <- c(pvals, pv)
                    prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                      nrow = 2)
                    prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                    dimnames(prdLabels$var[[Pfac]]) <- list(pv, 
                      Pfac)
                  }
                  else {
                    pv <- inter$coords$nuxPoints[[Pfac]]
                    pvals <- c(pvals, pv)
                    prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                      nrow = 2)
                    prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                    dimnames(prdLabels$var[[Pfac]]) <- list(pv, 
                      Pfac)
                  }
                }
                else {
                  dat <- dataFrame[, inter$inter[3, spline]]
                  pv <- c(min(dat, na.rm = TRUE), mean(dat, na.rm = TRUE), 
                    max(dat, na.rm = TRUE))
                  pvals <- c(pvals, pv)
                  prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                    nrow = 2)
                  prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                  dimnames(prdLabels$var[[Pfac]]) <- list(pv, 
                    Pfac)
                }
            }
            else if (inter$inter[1, ij] == -8) {
                pvals <- c(pvals, 1)
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(1, ncol = 1)
                dimnames(prdLabels$var[[Pfac]]) <- list("(Intercept)", 
                  "(Intercept)")
            }
            else if (inter$flev[ij] == 1 & inter$inter[4, ij] <= 
                0) {
                if (length(levels(dataFrame[, inter$inter[3, 
                  ij]])) > 0) {
                  stop(inter$facnam[ij], "has only one level", 
                    "\n")
                }
                pv <- mean(dataFrame[, inter$inter[3, ij]], na.rm = TRUE)
                pvals <- c(pvals, pv)
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                dimnames(prdLabels$var[[Pfac]]) <- list(pv, Pfac)
            }
            else if (inter$inter[4, ij] == -1) {
                pv <- seq(1, inter$flev[ij] + 1)
                pvals <- c(pvals, pv)
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                lvCon <- levels(dataFrame[, inter$inter[3, ij]])
                dimnames(prdLabels$var[[Pfac]]) <- list(lvCon, 
                  Pfac)
            }
            else if (ranfac < 0 || fixfac < 0) {
                if (length(inter$coords[[Pfac]]) > 0) {
                  which <- inter$varLevels[[inter$facnam[ij]]]
                  what <- matrix(inter$coords[[Pfac]][1, ], nrow = 1)
                }
                else {
                  which <- inter$varLevels[[inter$facnam[ij]]]
                  what <- matrix(apply(dataFrame[, which], 2, 
                    mean, na.rm = TRUE), nrow = 1)
                }
                dimnames(what) <- list("Average", which)
                pvals <- c(pvals, what[1, ])
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- what
            }
            else {
                pv <- seq(1, length(inter$varLevels[[PVfac]]))
                pvals <- c(pvals, pv)
                prdLabels$var[[Pfac]] <- matrix(0, ncol = 2, 
                  nrow = 2)
                prdLabels$var[[Pfac]] <- matrix(pv, ncol = 1)
                dimnames(prdLabels$var[[Pfac]]) <- list(inter$varLevels[[PVfac]], 
                  Pfac)
            }
        }
        prOrd <- prOrd + 1
        buff[8] <- prOrd
        buff[7] <- npfact[7, nfinmd + npcol]
        buff[2] <- lpvals
        buff[3] <- length(pvals) - lpvals + 1
        lpvals <- length(pvals) + 1
        npfact[, nfinmd + npcol] <- buff
    }
    xLevs <- sapply(prdLabels$var, function(x) {
        if (is.vector(x)) 
            return(length(x))
        else return(nrow(x))
    })
    if (is.null(prdLabels$parallel)) 
        pll <- rep(FALSE, length(xLevs))
    else pll <- !is.na(match(names(prdLabels$var), prdLabels$parallel))
    howMany <- asreml.ie(sum(xLevs[!pll]) == 0, 0, prod(xLevs[!pll])) + 
        sum(xLevs[!pll]) * prdLabels$margin + asreml.ie(length(xLevs[pll]) == 
        0, 0, xLevs[pll][1])
    if (length(parList) > 0 & (nptabl[4] == 0)) 
        nptabl[4] <- match(parList[1], dimnames(baseFac)[[2]])
    ptrPrdvals[1] <- lenPrdvals + 1
    ptrPrdvals[2] <- howMany
    lenPrdvals <- lenPrdvals + howMany
    if (buff16[6] == 1 | buff16[6] == 3) {
        ptrVpvals[1] <- vPtr
        ptrVpvals[2] <- howMany * (howMany + 1)/2
        vPtr <- vPtr + ptrVpvals[2]
        lenVpvals <- vPtr - 1
    }
    if (buff16[6] == 2 | buff16[6] == 3) {
        ptrVpvals[3] <- vPtr
        ptrVpvals[4] <- howMany * (howMany - 1)/2
        vPtr <- vPtr + ptrVpvals[4]
        lenVpvals <- vPtr - 1
    }
    if (buff16[6] == 0) {
        ptrVpvals[1] <- vPtr
        ptrVpvals[2] <- 0
        ptrVpvals[3] <- vPtr
        ptrVpvals[4] <- 0
    }
    buff <- rep(0, 8)
    prsList <- predict$present
    if (length(prsList) > 0) {
        tmpList <- list()
        if (is.list(prsList)) 
            tmpList <- prsList
        else {
            tmpList[[1]] <- prsList
            names(tmpList) <- Pf
        }
        for (L in 1:length(tmpList)) {
            if (casefold(names(tmpList)[L]) == "prwts") {
                which <- seq(1, length(tmpList))[-L][1]
                ij <- match(tmpList[[which]], inter$facnam)
                if (any(is.na(ij))) 
                  stop("Not all factors in present list are in the model\n")
                kk <- prod(inter$flev[ij])
                if (length(tmpList[[L]]) != kk) 
                  stop(paste("Expected", kk, "weights for present list", 
                    "\n"))
                pvals <- c(pvals, tmpList[[L]])
                nptabl[8] <- lpvals
                lpvals <- length(pvals) + 1
            }
            else {
                for (Pfac in tmpList[[L]]) {
                  if (is.na(match(Pfac, inter$facnam))) 
                    stop(paste("\nPresent: ", Pfac, " is not in the model\n"))
                  Plevs <- NULL
                  if (!is.null(predict$levels[[Pfac]])) {
                    if (is.numeric(predict$levels[[Pfac]])) 
                      Plevs <- predict$levels[[Pfac]]
                    else Plevs <- match(predict$levels[[Pfac]], 
                      inter$varLevels[[Pfac]])
                  }
                  npcol <- match(Pfac, dimnames(baseFac)[[2]])
                  where <- nfinmd + npcol
                  buff <- npfact[, where]
                  buff[4] <- -match(Pfac, tmpList[[L]]) - ((L - 
                    1) * 100)
                  if (buff[1] == 0) {
                    ij <- match(Pfac, inter$facnam)
                    if (is.null(Plevs)) 
                      Plevs <- seq(1, inter$flev[ij])
                    pvals <- c(pvals, Plevs)
                    buff[2] <- lpvals
                    buff[3] <- length(pvals) - lpvals + 1
                    lpvals <- length(pvals) + 1
                  }
                  npfact[, where] <- buff
                }
            }
        }
    }
    prsList <- predict$average
    lvec <- 0
    if (length(prsList) > 0) {
        for (L in 1:length(prsList)) {
            Plevs <- NULL
            if (is.list(prsList)) {
                if (is.list(prsList[[L]])) 
                  stop("Average must be a vector of variate names or a list of weight vectors\n")
                Pfac <- names(prsList)[L]
                Plevs <- prsList[[L]]
                lvec <- length(Plevs)
            }
            else Pfac <- prsList[L]
            if (is.na(match(Pfac, inter$facnam))) 
                stop(paste("\nAverage: ", Pfac, " is not in the model\n"))
            ij <- match(Pfac, inter$facnam)
            npcol <- match(Pfac, dimnames(baseFac)[[2]])
            if (is.na(npcol)) {
                npcol <- match(baseOrd$bFac[match(Pfac, baseOrd$bFun)], 
                  dimnames(baseFac)[[2]])
                ijBase <- match(baseOrd$bFac[match(Pfac, baseOrd$bFun)], 
                  inter$facnam)
            }
            else ijBase <- ij
            where <- nfinmd + npcol
            buff <- npfact[, where]
            if (is.null(Plevs)) {
                Plevs <- c(rep(1, inter$flev[ijBase])/inter$flev[ijBase], 
                  rep(0, max(0, inter$flev[ij] - inter$flev[ijBase])))
            }
            if (length(Plevs) != inter$flev[ij]) {
                warning(paste(Pfac, "has", inter$flev[ij], "levels; extra weights ignored", 
                  "or missing weights treated as zero.\n"))
                Plevs <- c(Plevs, rep(0, max(0, inter$flev[ij] - 
                  length(Plevs))))
            }
            pvals <- c(pvals, Plevs[1:inter$flev[ij]])
            buff[4] <- lpvals
            buff[5] <- lvec
            lpvals <- length(pvals) + 1
            npfact[, where] <- buff
        }
    }
    who <- match(colNames, Pnames)
    miss <- (seq(1, length(colNames)) + i * nfinmd)[is.na(who)]
    pres <- (seq(1, length(colNames)) + i * nfinmd)[!is.na(who)]
    npfact[8, seq(nfinmd + 1, 2 * nfinmd)] <- seq(1, nfinmd)
    names(ptrPrdvals) <- c("loc", "len")
    for (i in 1:nfinmd) npfact[, i] <- baseFac[, i]
    npfact[8, 1:nfinmd] <- 1:nfinmd
    if (length(ivals) == 0) 
        ivals <- 0
    lpwts <- length(pvals) + 1
    pvals <- c(pvals, rep(0, (2 * nfinmd + 1)))
    TT <- list(npred = npred, prdLabels = prdLabels, nfinmd = nfinmd, 
        baseFac = baseFac, npfact = npfact, nptabl = nptabl, 
        ivals = ivals, pvals = pvals, lpwts = lpwts, lenPrdvals = lenPrdvals, 
        ptrPrdvals = ptrPrdvals, lenVpvals = lenVpvals, ptrVpvals = ptrVpvals, 
        predict = predict)
    TT
}
asreml.randomStr <-
function (random, data) 
{
    if (length(random[[2]]) == 0) 
        return(list(form = random, dummyTerms = character(0)))
    tt <- terms(random, specials = asreml.Spcls$Fun, keep.order = TRUE)
    which <- attr(tt, "specials")$str
    if (all(is.null(which))) 
        return(list(form = random, dummyTerms = character(0)))
    terms.fac <- attr(tt, "factors")
    terms.var <- dimnames(terms.fac)[[1]]
    terms.lab <- dimnames(terms.fac)[[2]]
    str <- numeric(0)
    for (w in which) str <- c(str, seq(1, length(terms.lab))[terms.fac[w, 
        ] > 0])
    terms.use <- character(0)
    dummy <- character(0)
    nLevels <- numeric(0)
    kk <- 0
    for (i in 1:length(terms.lab)) {
        if (any(str == i)) {
            tts <- terms(formula(paste("~", terms.lab[i])), specials = asreml.Spcls$Fun, 
                keep.order = TRUE)
            trm <- dimnames(attr(tts, "factors"))[[1]]
            ww <- attr(tts, "specials")$str
            if (length(ww) > 1) 
                stop("Direct product of 'str' terms not implemented")
            wa <- attr(tts, "specials")$at
            if (length(ww) + length(wa) != length(trm)) 
                stop("Direct product with 'str' not implemented")
            atTrm <- character(0)
            if (length(wa)) {
                if (length(wa) > 1) 
                  stop("Multiple at' terms not implemented")
                Y <- eval(asreml.Eval((dimnames(attr(tts, "factors"))[[1]])[wa], 
                  data))
                atTrm <- paste("(", paste(paste("at(", Y$Obj, 
                  ", ", "'", Y$Lvls, "'", ")", sep = ""), collapse = "+"), 
                  ")")
            }
            Y <- eval(asreml.Eval((dimnames(attr(tts, "factors"))[[1]])[ww], 
                data))
            tt.use <- Y$Obj
            if (any(!is.na(fa <- match("fa", names(Y$Lvls))))) {
                if (length(fa) > 1) 
                  stop("fa x fa not implemented")
                k <- asreml.fak(Y$Lvls[[fa]], Y$Lab[[fa]])
                for (j in 1:k) {
                  lev <- prod(unlist(lapply(Y$Lvls, function(x) length(x)))[-fa])
                  dum <- paste(paste("fa", lev, sep = ""), kk, 
                    j, sep = ".")
                  tt.use <- c(tt.use, dum)
                  dummy <- c(dummy, dum)
                  nLevels <- c(nLevels, lev)
                }
            }
            kk <- kk + 1
            if (length(atTrm)) 
                tt.use <- paste(atTrm, paste("(", paste(tt.use, 
                  collapse = "+"), ")", sep = ""), sep = ":")
            else tt.use <- paste(tt.use, collapse = "+")
            terms.use <- c(terms.use, tt.use)
        }
        else terms.use <- c(terms.use, terms.lab[i])
    }
    nLevels <- as.list(nLevels)
    names(nLevels) <- dummy
    return(list(form = asreml.formula(parse(text = paste("~", 
        paste(terms.use, collapse = "+")))), dummyTerms = nLevels))
}
asreml.rdflt <-
function (y, rcov = ~units, dataFrame = sys.parent, dispersion = NA) 
{
    asreml.rdflt2 <- function(Y, Var, sectioname, section, dataFrame, 
        y, dispersion) {
        if (length(names(dispersion)) > 0 && names(dispersion) == 
            "multinomial") 
            Nresid <- prod(sapply(Var, function(v, Y, data) {
                length(levels(data[[Y[[v]]$Obj]]))
            }, Y, dataFrame))
        else Nresid <- prod(sapply(Var, function(v, Y, data) {
            length(unique(data[[Y[[v]]$Obj]]))
        }, Y, dataFrame))
        if (Nresid != nrow(dataFrame)) 
            stop(paste("Error model specifies", Nresid, "observations, but data section has", 
                nrow(dataFrame), "rows\n"))
        asreml.chkOrd(Y, Var, dataFrame)
        ndim <- length(Var)
        level2 <- vector(mode = "list", length = ndim + 1)
        if (length(y) == 1) 
            scale <- var(dataFrame[[y]][!is.na(dataFrame[[y]])])/2
        else scale <- 1
        Term <- paste(unique(c(sectioname, section)), collapse = "_")
        isVar <- 0
        temp <- NULL
        for (i in seq(1, ndim)) {
            X <- Y[[Var[i]]]
            z <- ifelse((is.na(X$isVariance) || !X$isVariance), 
                1, scale)
            X$Initial[X$Tgamma == 2] <- X$Initial[X$Tgamma == 
                2] * z
            if (length(names(dispersion)) > 0 && names(dispersion) == 
                "multinomial") {
                if (X$Fun == "us") {
                  nn <- (sqrt(8 * length(X$Lab) + 1) - 1)/2
                  X$Initial <- diag(1, nn)[t(lower.tri(diag(1, 
                    nn), TRUE))]
                }
                else if (X$Fun == "idv") {
                  X$Initial <- 1
                }
                X$Con <- rep("F", length(X$Con))
                names(X$Initial) <- names(X$Con) <- X$Lab
            }
            level2[[i]]$levels <- X$Lvls
            level2[[i]]$initial <- X$Initial
            level2[[i]]$con <- X$Con
            level2[[i]]$model <- X$Fun
            isVar <- isVar + ifelse(is.na(X$isVariance), 0, as.numeric(X$isVariance))
            names(level2[[i]]$initial) <- paste(Term, X$Lab, 
                sep = "")
            temp <- c(temp, X$Call)
        }
        if (is.na(dispersion)) {
            if (isVar > 0) {
                level2[[ndim + 1]]$s2 <- 1
                level2[[ndim + 1]]$con <- "F"
            }
            else {
                level2[[ndim + 1]]$s2 <- scale
                level2[[ndim + 1]]$con <- "P"
            }
        }
        else if (dispersion == 0 | dispersion == -1) {
            level2[[ndim + 1]]$s2 <- 1
            level2[[ndim + 1]]$con <- ifelse(isVar == 0, "P", 
                "F")
        }
        else {
            level2[[ndim + 1]]$s2 <- abs(dispersion)
            level2[[ndim + 1]]$con <- "F"
        }
        names(level2[[ndim + 1]]$s2) <- paste(Term, "!", "variance", 
            sep = "")
        temp <- c(temp, "variance")
        names(level2) <- temp
        return(level2)
    }
    asreml.chkOrd <- function(Y, Var, data) {
        objs <- sapply(Var, function(v, Y) {
            Y[[v]]$Obj
        }, Y)
        dd <- list()
        for (x in objs) {
            if (is.factor(data[, x])) {
                ok <- 0
                out <- .C("canBeNumeric", as.character(data[, 
                  x]), as.integer(length(data[, x])), ok = as.integer(ok))
                if (out$ok == 1) 
                  dd[[x]] <- as.numeric(as.character(data[, x]))
                else dd[[x]] <- match(as.character(data[, x]), 
                  unique(as.character(data[, x])))
            }
            else dd[[x]] <- data[, x]
        }
        idx <- do.call("order", dd)
        if (any(diff(idx) < 0)) 
            stop("Data frame order does not match that specified by the R model formula.")
        invisible(0)
    }
    if (!inherits(rcov, "formula")) 
        stop("\nrcov must be a formula")
    if (length(rcov) != 2) 
        stop("\nRcov model formula must be of the form \" ~ pred\"")
    if (asreml.Rsys) {
        tt <- lapply(asreml.subForm(paste(deparse(rcov[[2]]), 
            collapse = "")), function(x) terms(x, specials = asreml.Spcls$Fun))
        terms.var <- unique(unlist(lapply(tt, function(x) dimnames(attr(x, 
            "factors"))[[1]])))
        terms.lab <- unique(unlist(lapply(tt, function(x) dimnames(attr(x, 
            "factors"))[[2]])))
    }
    else {
        tt <- terms(rcov, specials = asreml.Spcls$Fun)
        terms.fac <- attr(tt, "factors")
        terms.var <- dimnames(terms.fac)[[1]]
        terms.lab <- dimnames(terms.fac)[[2]]
    }
    nlab <- length(terms.lab)
    nvar <- length(terms.var)
    Y <- lapply(terms.var, function(x, dataFrame) {
        y <- eval(asreml.Eval(x, dataFrame, Rcov = 1))
        if (is.na(match("Fun", names(y)))) {
            y <- list(Fun = NULL, isVariance = NA)
        }
        return(y)
    }, dataFrame)
    names(Y) <- terms.var
    section <- NULL
    for (z in terms.var) {
        X <- Y[[z]]
        if (inherits(X, "asreml.special") && X$Fun == "at") {
            section <- c(section, X$Obj)
        }
    }
    if (is.null(section)) {
        nsect <- 1
        level1 <- vector(mode = "list", length = 1)
        section <- "R"
        names(level1) <- section
    }
    else if (length(unique(section)) > 1) 
        stop("More than one section named by at()\n")
    else {
        section <- section[1]
        nsect <- length(unique(dataFrame[[section]]))
        level1 <- vector(mode = "list", length = nsect)
        names(level1) <- unique(as.character(dataFrame[[section]]))
    }
    if (nsect > 1) {
        dd <- match(as.character(dataFrame[, section]), unique(as.character(dataFrame[, 
            section])))
        if (any(diff(order(dd)) < 0)) 
            stop(paste("Data frame not grouped by", section))
    }
    for (z in seq(1, nlab)) {
        if (asreml.Rsys) 
            Var <- dimnames(attr(tt[[z]], "factors"))[[1]]
        else Var <- dimnames(attr(tt[z], "factors"))[[1]]
        which <- sapply(Var, function(v, y) {
            ifelse(!is.null(y[[v]]$Fun) && y[[v]]$Fun == "at", 
                TRUE, FALSE)
        }, Y)
        which.id <- sapply(Var, function(v, y) {
            ifelse(is.null(y[[v]]$Fun) || is.na(y[[v]]$isVariance), 
                TRUE, FALSE)
        }, Y)
        if (sum(which.id) > 0) {
            v <- Var
            v[which.id] <- paste("id(", Var[which.id], ")", sep = "")
            for (i in seq(1, length(Var))[which.id]) Y[[Var[i]]] <- eval(asreml.Eval(v[i], 
                dataFrame, Rcov = 1))
        }
        if (sum(which) > 1) 
            stop("Illegal use of at()\n")
        else if (sum(which) == 1) {
            X <- Y[[Var[which]]]
            for (lvl in X$Lvls) {
                for (what in Var[!which]) {
                  if (Y[[what]]$Fun != "id") {
                    simple <- dataFrame[dataFrame[[X$Obj]] == 
                      lvl, ]
                    Y[[what]] <- eval(asreml.Eval(what, simple, 
                      Rcov = 1))
                  }
                }
                which.row <- (dataFrame[[X$Obj]] == lvl)
                level1[[lvl]] <- asreml.rdflt2(Y, Var[!which], 
                  section, lvl, dataFrame[which.row, ], y, dispersion)
            }
        }
        else {
            if (nlab > 1) 
                stop("Missing at() to identify error sections\n")
            if (is.na(dispersion)) 
                dispersion <- 0
            level1[[1]] <- asreml.rdflt2(Y, Var, section, section, 
                dataFrame, y, dispersion)
        }
    }
    ldim <- unlist(lapply(level1, length))
    nz <- ldim[ldim != 0]
    which <- seq(1, nsect)[ldim != 0]
    if (length(which) == 0) 
        stop("No error dimensions\n")
    if (any(nz != nz[1])) 
        stop("Unequal numbers of dimensions\n")
    if (any(ldim == 0)) {
        Var <- names(level1[[which[1]]])
        Var <- paste("id(", Var[-length(Var)], ")", sep = "")
        YY <- lapply(Var, function(x, dataFrame) {
            eval(asreml.Eval(x, dataFrame, Rcov = 1))
        }, dataFrame)
        names(YY) <- Var
        for (i in seq(1, nsect)[ldim == 0]) {
            which.row <- (dataFrame[[section]] == names(level1)[i])
            level1[[i]] <- asreml.rdflt2(YY, Var, section, names(level1)[i], 
                dataFrame[which.row, ], y, dispersion)
        }
    }
    return(level1)
}
asreml.read.table <-
function (...) 
{
    x <- read.table(...)
    if (any(dim(x) == 0)) 
        stop("Empty data frame")
    which <- .Call("a2i", x)
    if (length(which) > 0) {
        for (i in which) {
            x[[i]] <- factor(x[[i]])
        }
    }
    invisible(x)
}
asreml.rlist <-
function (rcov, dataFrame) 
{
    asreml.rlist2 <- function(Y, Var, sectioname, section, dataFrame) {
        ndim <- length(Var)
        level2 <- vector(mode = "list", length = ndim + 1)
        temp <- NULL
        for (i in seq(1, ndim)) {
            X <- Y[[Var[i]]]
            level2[[i]]$length <- length(unique(dataFrame[[X$Obj]]))
            level2[[i]]$struc <- X$Struc
            level2[[i]]$tgamma <- abs(X$Tgamma)
            if (is.matrix(X$Coords$coords) && nrow(X$Coords$coords) > 
                1) 
                level2[[i]]$coords <- X$Coords$coords[, rank(dataFrame$units)]
            else level2[[i]]$coords <- X$Coords$coords
            temp <- c(temp, X$Call)
        }
        level2[[ndim + 1]]$size <- nrow(dataFrame)
        temp <- c(temp, "size")
        names(level2) <- temp
        return(level2)
    }
    if (asreml.Rsys) {
        tt <- lapply(asreml.subForm(paste(deparse(rcov[[2]]), 
            collapse = "")), function(x) terms(x, specials = asreml.Spcls$Fun))
        terms.var <- unique(unlist(lapply(tt, function(x) dimnames(attr(x, 
            "factors"))[[1]])))
        terms.lab <- unique(unlist(lapply(tt, function(x) dimnames(attr(x, 
            "factors"))[[2]])))
    }
    else {
        tt <- terms(rcov, specials = asreml.Spcls$Fun)
        terms.fac <- attr(tt, "factors")
        terms.var <- dimnames(terms.fac)[[1]]
        terms.lab <- dimnames(terms.fac)[[2]]
    }
    nlab <- length(terms.lab)
    nvar <- length(terms.var)
    Y <- lapply(terms.var, function(x, dataFrame) {
        y <- eval(asreml.Eval(x, dataFrame, Rcov = 1))
        if (is.na(match("Fun", names(y)))) {
            y <- list(Fun = NULL, isVariance = NA)
        }
        return(y)
    }, dataFrame)
    names(Y) <- terms.var
    section <- NULL
    for (z in terms.var) {
        X <- Y[[z]]
        if (inherits(X, "asreml.special") && X$Fun == "at") {
            section <- c(section, X$Obj)
        }
    }
    if (is.null(section)) {
        nsect <- 1
        level1 <- vector(mode = "list", length = 1)
        section <- "R"
        names(level1) <- section
    }
    else if (length(unique(section)) > 1) 
        stop("More than one section named by at()\n")
    else {
        section <- section[1]
        nsect <- length(unique(dataFrame[[section]]))
        level1 <- vector(mode = "list", length = nsect)
        names(level1) <- unique(as.character(dataFrame[[section]]))
    }
    for (z in seq(1, nlab)) {
        if (asreml.Rsys) 
            Var <- dimnames(attr(tt[[z]], "factors"))[[1]]
        else Var <- dimnames(attr(tt[z], "factors"))[[1]]
        which <- sapply(Var, function(v, y) {
            ifelse(!is.null(y[[v]]$Fun) && y[[v]]$Fun == "at", 
                TRUE, FALSE)
        }, Y)
        which.id <- sapply(Var, function(v, y) {
            ifelse(is.null(y[[v]]$Fun) || is.na(y[[v]]$isVariance), 
                TRUE, FALSE)
        }, Y)
        if (sum(which.id) > 0) {
            v <- Var
            v[which.id] <- paste("id(", Var[which.id], ")", sep = "")
            for (i in seq(1, length(Var))[which.id]) Y[[Var[i]]] <- eval(asreml.Eval(v[i], 
                dataFrame, Rcov = 1))
        }
        if (sum(which) > 1) 
            stop("Ilegal use of at()\n")
        else if (sum(which) == 1) {
            X <- Y[[Var[which]]]
            for (lvl in X$Lvls) {
                for (what in Var[!which]) {
                  if (Y[[what]]$Fun != "id") {
                    simple <- dataFrame[dataFrame[[X$Obj]] == 
                      lvl, ]
                    Y[[what]] <- eval(asreml.Eval(what, simple, 
                      Rcov = 1))
                  }
                }
                which.row <- (dataFrame[[X$Obj]] == lvl)
                level1[[lvl]] <- asreml.rlist2(Y, Var[!which], 
                  section, lvl, dataFrame[which.row, ])
            }
        }
        else {
            if (nlab > 1) 
                stop("Missing at() to identify error sections\n")
            level1[[1]] <- asreml.rlist2(Y, Var, section, section, 
                dataFrame)
        }
    }
    ldim <- unlist(lapply(level1, length))
    nz <- ldim[ldim != 0]
    which <- seq(1, nsect)[ldim != 0]
    if (length(which) == 0) 
        stop("No error dimensions\n")
    if (any(nz != nz[1])) 
        stop("Unequal numbers of dimensions\n")
    if (any(ldim == 0)) {
        Var <- names(level1[[which[1]]])
        Var <- paste("id(", Var[-length(Var)], ")", sep = "")
        YY <- lapply(Var, function(x, dataFrame) {
            eval(asreml.Eval(x, dataFrame, Rcov = 1))
        }, dataFrame)
        names(YY) <- Var
        for (i in seq(1, nsect)[ldim == 0]) {
            which.row <- (dataFrame[[section]] == names(level1)[i])
            level1[[i]] <- asreml.rlist2(YY, Var, section, names(level1[i]), 
                dataFrame[which.row, ])
        }
    }
    return(level1)
}
asreml.rr <-
function (obj, k = 1, init = NA, data, ...) 
{
    out <- vector(mode = "list", length = 13)
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) {
        out[[1]] <- "rr"
        out[[2]] <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- c(obj[[4]], paste("Comp", seq(1, k), sep = ""))
        n <- length(out[[4]]) - k
        N <- k * n + n
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0, n), rep(0.1, n * k))
        }
        if (!missing(init)) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("rr(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("F", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", temp, ".var", sep = ""), paste(paste(rep(paste("!", 
            temp, sep = ""), k), "rr", sep = "."), rep(seq(1, 
            k), rep(n, k)), sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        temp[seq(1, length(con))[con == "P"]] <- 1
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n + k, 0, NA, 0)
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
        out[[13]] <- TRUE
    }
    else if (is.numeric(substitute(obj))) {
        out[[1]] <- "rr"
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        n <- obj
        N <- k * n + n
        lvls <- seq(1, obj + k)
        out[[4]] <- lvls
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0, n), rep(0.1, n * k))
        }
        if (IniFlag) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("rr(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("F", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", temp, ".var", sep = ""), paste(paste(rep(paste("!", 
            temp, sep = ""), k), "rr", sep = "."), rep(seq(1, 
            k), rep(n, k)), sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(obj, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 
            0, 0, 0, 0)
        out[[10]] <- c(obj, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    else {
        out[[1]] <- "rr"
        obj <- as.character(substitute(obj))
        out[[2]] <- asreml.spCall(sys.call())
        out[[3]] <- obj
        drop.unused.levels <- attr(data, "Control")$drop.unused.levels
        data <- data[[obj]]
        if (is.factor(data)) 
            out[[4]] <- c(asreml.levels(data, drop.unused.levels), 
                paste("Comp", seq(1, k), sep = ""))
        else stop(paste(obj, "must be a factor\n"))
        n <- length(out[[4]]) - k
        N <- k * n + n
        IniFlag <- TRUE
        if (missing(init)) {
            IniFlag <- FALSE
            init <- c(rep(0, n), rep(0.1, n * k))
        }
        if (IniFlag) {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != N) 
                stop("rr(): Wrong number of initial values\n")
        }
        if (k > 1 & !IniFlag) 
            init[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- 0
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            con <- c(rep("F", n), rep("U", k * n))
        if (k > 1) 
            con[c(rep(FALSE, n), !lower.tri(matrix(nrow = n, 
                ncol = k), diag = TRUE))] <- "F"
        out[[6]] <- con
        temp <- out[[4]][1:n]
        out[[7]] <- c(paste("!", obj, ".", temp, ".var", sep = ""), 
            paste(paste(rep(paste("!", obj, ".", temp, sep = ""), 
                k), "rr", sep = "."), rep(seq(1, k), rep(n, k)), 
                sep = ""))
        temp <- c(rep(1, n), rep(6, k * n))
        out[[8]] <- temp
        if (IniFlag) 
            out[[8]] <- -out[[8]]
        out[[9]] <- c(0, 0, 15, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(n + k, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- TRUE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.Rsys <-
TRUE
asreml.rupdt <-
function (rcov.param, gammas, gmcstr = NULL, VCov = FALSE) 
{
    if (length(rcov.param) == 0) 
        return(NULL)
    rupdt <- lapply(rcov.param, function(x, gammas, gmcstr, VCov) {
        lapply(x, function(y, gammas, gmcstr, VCov) {
            S2 <- FALSE
            if (j <- match("initial", names(y), nomatch = 0)) 
                tmp <- names(y[[j]])
            else if (j <- match("s2", names(y), nomatch = 0)) {
                S2 <- TRUE
                tmp <- names(y[[j]])
            }
            else stop("Cannot match initial or s2\n")
            i <- match(tmp, names(gammas), nomatch = NA)
            ii <- !is.na(i)
            i <- i[ii]
            if (length(i) > 0) {
                y[[j]][ii] <- gammas[i]
                if (!S2 && VCov) {
                  y[[j]] <- switch(y$model, ante = asreml.udu(y[[j]]), 
                    chol = asreml.ldl(y[[j]]), cholc = asreml.ldl(y[[j]]), 
                    y[[j]])
                  if (!is.null(gmcstr)) {
                    if (is.element(y$model, c("ante", "chol", 
                      "cholc"))) 
                      y$con <- y$con
                    else y$con[ii] <- casefold(substr(gmcstr[i], 
                      1, 1), upper = TRUE)
                  }
                }
                else if (!is.null(gmcstr)) 
                  y$con[ii] <- casefold(substr(gmcstr[i], 1, 
                    1), upper = TRUE)
                names(y[[j]]) <- tmp
            }
            y
        }, gammas, gmcstr, VCov)
    }, gammas, gmcstr, VCov)
    rupdt
}
asreml.sar <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sar", "cor", 0.1, numeric(0), NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out
}
asreml.sar2 <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sar2", "cor", "rep(0.1,2)", numeric(0), 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out
}
asreml.sar2h <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sar2h", "cor", "rep(0.1,2)", "rep(0.1,n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out
}
asreml.sar2v <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sar2v", "cor", "rep(0.1,2)", 0.1, NA, 
        obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 
        0, 0)
    out
}
asreml.sarh <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sarh", "cor", 0.1, "rep(0.1,n)", NA, obj, 
        init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out
}
asreml.sarv <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("sarv", "cor", 0.1, 0.1, NA, obj, init, 
        data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 14, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 
        0, 0)
    out
}
asreml.sectionName <-
function (rcov, data) 
{
    tt <- terms(rcov, specials = asreml.Spcls$Fun)
    terms.fac <- attr(tt, "factors")
    terms.var <- dimnames(terms.fac)[[1]]
    terms.lab <- dimnames(terms.fac)[[2]]
    nvar <- length(terms.var)
    section <- NULL
    for (z in seq(1, nvar)) {
        X <- eval(asreml.Eval(terms.var[z], data))
        if (inherits(X, "asreml.special")) {
            if (X$Fun == "at") 
                section <- c(section, X$Obj)
        }
    }
    if (!is.null(section)) {
        if (length(unique(section)) > 1) 
            stop("More than one section named by at()\n")
        else section <- section[1]
    }
    return(section)
}
asreml.setG <-
function (Gdflt, Gparam, update.con) 
{
    if (mode(Gparam) == "character") {
        if (asreml.Rsys) 
            Gtab <- read.table(Gparam, header = TRUE, sep = ",", 
                as.is = TRUE)
        else Gtab <- importData(Gparam, type = "ASCII", delimiter = ",", 
            stringsAsFactors = FALSE)
        if (any(is.na(match(c("Gamma", "Value", "Constraint"), 
            names(Gtab))))) 
            stop(paste(Gparam, "must have components 'Gamma','Value','Constraint'"))
        gammas <- Gtab$Value
        gmcstr <- Gtab$Constraint
        names(gammas) <- Gtab$Gamma
        names(gmcstr) <- Gtab$Gamma
        return(asreml.gupdt(Gdflt, gammas, gmcstr))
    }
    if (is.data.frame(Gparam)) {
        if (any(is.na(match(c("Gamma", "Value", "Constraint"), 
            names(Gparam))))) 
            stop(paste(Gparam, "must have components 'Gamma','Value','Constraint'"))
        gammas <- Gparam$Value
        gmcstr <- Gparam$Constraint
        names(gammas) <- Gparam$Gamma
        names(gmcstr) <- Gparam$Gamma
        return(asreml.gupdt(Gdflt, gammas, gmcstr))
    }
    if (mode(Gparam) == "list") {
        if (!is.na(match("G.param", names(Gparam)))) 
            Gparam <- Gparam[["G.param"]]
    }
    gammas <- vector(mode = "numeric")
    gmcstr <- vector(mode = "character")
    if (length(Gparam) > 0) {
        nterm <- length(Gparam)
        for (n in seq(1, nterm)) {
            nfact <- length(Gparam[[n]])
            for (j in 1:nfact) {
                gammas <- c(gammas, Gparam[[n]][[j]]$initial)
                gmcstr <- c(gmcstr, Gparam[[n]][[j]]$con)
            }
        }
    }
    if (length(gammas) == 0) 
        return(Gdflt)
    gammas <- gammas[!is.na(gammas)]
    gmcstr <- gmcstr[gmcstr != ""]
    if (update.con) 
        return(asreml.gupdt(Gdflt, gammas, gmcstr, FALSE))
    else return(asreml.gupdt(Gdflt, gammas, NULL, FALSE))
}
asreml.setR <-
function (Rdflt, Rparam, update.con) 
{
    if (mode(Rparam) == "character") {
        if (asreml.Rsys) 
            Rtab <- read.table(Rparam, header = TRUE, sep = ",", 
                as.is = TRUE)
        else Rtab <- importData(Rparam, type = "ASCII", delimiter = ",", 
            stringsAsFactors = FALSE)
        if (any(is.na(match(c("Gamma", "Value", "Constraint"), 
            names(Rtab))))) 
            stop(paste(Rparam, "must have components 'Gamma','Value','Constraint'"))
        gammas <- Rtab$Value
        gmcstr <- Rtab$Constraint
        names(gammas) <- Rtab$Gamma
        names(gmcstr) <- Rtab$Gamma
        return(asreml.rupdt(Rdflt, gammas, gmcstr))
    }
    if (is.data.frame(Rparam)) {
        if (any(is.na(match(c("Gamma", "Value", "Constraint"), 
            names(Rparam))))) 
            stop(paste(Rparam, "must have components 'Gamma','Value','Constraint'"))
        gammas <- Rparam$Value
        gmcstr <- Rparam$Constraint
        names(gammas) <- Rparam$Gamma
        names(gmcstr) <- Rparam$Gamma
        return(asreml.rupdt(Rdflt, gammas, gmcstr))
    }
    if (mode(Rparam) == "list") {
        if (!is.na(match("R.param", names(Rparam)))) 
            Rparam <- Rparam[["R.param"]]
    }
    gammas <- vector(mode = "numeric")
    gmcstr <- vector(mode = "character")
    if (length(Rparam) > 0) {
        nsect <- length(Rparam)
        nspat <- max(sapply(Rparam, length)) - 1
        for (n in seq(1, nsect)) {
            ndim <- length(Rparam[[n]]) - 1
            gammas <- c(gammas, Rparam[[n]]$variance$s2)
            gmcstr <- c(gmcstr, Rparam[[n]]$variance$con)
            for (j in 1:ndim) {
                gammas <- c(gammas, Rparam[[n]][[j]]$initial)
                gmcstr <- c(gmcstr, Rparam[[n]][[j]]$con)
            }
        }
    }
    if (length(gammas) == 0) 
        stop("No variance components found.\n")
    gammas <- gammas[!is.na(gammas)]
    gmcstr <- gmcstr[gmcstr != ""]
    if (update.con) 
        return(asreml.rupdt(Rdflt, gammas, gmcstr, FALSE))
    else return(asreml.rupdt(Rdflt, gammas, NULL, FALSE))
}
asreml.sfa <-
function (obj, k = 1, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to sfa() must be a factor\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "sfa"
    obj <- as.character(substitute(obj))
    out[[2]] <- obj
    out[[3]] <- obj
    drop.unused.levels <- attr(data, "Control")$drop.unused.levels
    data <- data[[obj]]
    if (is.factor(data)) 
        out[[4]] <- asreml.levels(data, drop.unused.levels)
    else stop(paste(obj, "must be a factor\n"))
    n <- length(out[[4]])
    N <- k * n + n
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- c(rep(0.7, n * k), rep(0.1, n))
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != N) 
            stop("sfa(): Wrong number of initial values\n")
    }
    if (k > 1) 
        init[c(!lower.tri(matrix(nrow = n, ncol = k), diag = TRUE), 
            rep(FALSE, n))] <- 0
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) {
        temp <- rep("P", N)
        if (k > 1) 
            temp[c(!lower.tri(matrix(nrow = n, ncol = k), diag = TRUE), 
                rep(FALSE, n))] <- "F"
        out[[6]] <- temp
    }
    else out[[6]] <- con
    out[[7]] <- c(paste(paste(rep(paste("!", obj, ".", out[[4]], 
        sep = ""), k), "sfa", sep = "."), rep(seq(1, k), rep(n, 
        k)), sep = ""), paste("!", obj, ".", out[[4]], ".var", 
        sep = ""))
    out[[8]] <- c(rep(3, n * k), rep(2, n))
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 10, 0, 0, k, 0, 1, 0, 0, 0, 0, 0, 0, 
        0, 0)
    out[[10]] <- c(n, 0, NA, 0)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- obj
    out[[13]] <- TRUE
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.sparse2mat <-
function (x) 
{
    nrow <- max(x[, 1])
    ncol <- max(x[, 2])
    y <- rep(0, nrow * ncol)
    y[(x[, 2] - 1) * nrow + x[, 1]] <- x[, 3]
    y[(x[, 1] - 1) * nrow + x[, 2]] <- x[, 3]
    matrix(y, nrow = nrow, ncol = ncol, byrow = FALSE)
}
asreml.spc <-
function (fun, fam, ncor, nvar, k, obj, init, data, Rcov) 
{
    out <- vector(mode = "list", length = 13)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    emflag <- attr(data, "Control")$emflag
    if (fam == "mdl") {
        obj <- as.character(obj)
        out[[1]] <- fun
        out[[2]] <- asreml.spCall(sys.call(-1))
        out[[3]] <- obj
        if (is.factor(data[[obj]])) 
            out[[4]] <- asreml.levels(data[[obj]], as.logical(Rcov))
        else out[[4]] <- obj
        out[[5]] <- NA
        out[[6]] <- ""
        out[[7]] <- ""
        out[[8]] <- 2
        out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0)
        out[[10]] <- c(length(out[[4]]), 0, NA, -1)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
        out[[13]] <- NA
        oldClass(out) <- "asreml.special"
        return(out)
    }
    nest <- inherits(obj, "asreml.special")
    if (nest) {
        out[[1]] <- obj[[1]]
        out[[2]] <- obj[[2]]
        lab <- obj[[2]]
        out[[3]] <- obj[[3]]
        out[[4]] <- obj[[4]]
        out[[10]] <- obj[[10]]
        out[[11]] <- obj[[11]]
        out[[12]] <- obj[[12]]
    }
    else {
        if (is.numeric(obj)) {
            lvls <- seq(1, obj)
            lab <- paste(fun, "(", obj, ")", sep = "")
            cal <- asreml.spCall(sys.call(-1))
        }
        else {
            obj <- as.character(obj)
            cal <- obj
            if (is.factor(data[[obj]])) 
                lvls <- asreml.levels(data[[obj]], as.logical(Rcov))
            else if (fun == "id" || fun == "idv") 
                lvls <- obj
            else stop(paste(obj, "must be a factor"))
            lab <- obj
        }
        out[[1]] <- fun
        out[[2]] <- cal
        out[[3]] <- obj
        out[[4]] <- lvls
    }
    n <- length(out[[4]])
    what.var <- eval(parse(text = as.character(nvar)))
    if (is.vector(what.var)) 
        Nv <- length(what.var)
    else if (is.matrix(what.var)) 
        Nv <- n * (n + 1)/2
    else if (is.null(what.var)) 
        Nv <- 0
    else stop("Unrecognised mode for nvar")
    what.cor <- eval(parse(text = as.character(ncor)))
    if (is.vector(what.cor)) 
        Nc <- length(what.cor)
    else if (is.matrix(what.cor)) 
        Nc <- n * (n - 1)/2
    else if (is.null(what.cor)) 
        Nc <- 0
    else stop("Unrecognised mode for ncor")
    IniFlag <- ifelse(all(is.na(init)), FALSE, TRUE)
    if (IniFlag) {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != Nc + Nv) 
            stop(paste(fun, "- Wrong number of initial values\n"))
    }
    if (fam == "cor") {
        if (!IniFlag) {
            init <- c(rep(0.1, Nc), rep(0.1, Nv))
            if (length(init) == 0) 
                init <- NA
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- c(rep("U", Nc), rep("P", Nv))
        else out[[6]] <- con
        if (is.matrix(what.cor)) {
            x <- (row(what.cor) < col(what.cor))
            corLab <- c(paste(paste("!", lab, ".", out[[4]][col(x)[x]], 
                sep = ""), ":", paste("!", lab, ".", out[[4]][row(x)[x]], 
                sep = ""), ".cor", sep = ""))
        }
        else if (Nc > 1) 
            corLab <- paste("!", lab, paste(".cor", seq(1, Nc), 
                sep = ""), sep = "")
        else if (Nc == 1) 
            corLab <- paste("!", lab, ".cor", sep = "")
        else corLab <- character(0)
        if (Nv > 1) 
            varLab <- paste("!", lab, ".", out[[4]], sep = "")
        else if (Nv == 1) 
            varLab <- paste("!", lab, ".var", sep = "")
        else if (Nv + Nc == 0) 
            varLab <- paste("!", lab, ".id", sep = "")
        else varLab <- character(0)
        out[[7]] <- c(corLab, varLab)
        out[[8]] <- c(rep(3, Nc), rep(2, Nv))
        out[[13]] <- ifelse(Nv == 0, FALSE, TRUE)
    }
    else if (fam == "var") {
        if (is.matrix(what.var)) {
            if (!IniFlag) 
                init <- what.var[t(lower.tri(what.var, diag = TRUE))]
            out[[5]] <- init
            con <- asreml.matchCon(init)
            if (length(con) == 0) 
                out[[6]] <- rep("U", Nv)
            else out[[6]] <- con
            if (fun == "us" && emflag) 
                out[[6]] <- rep("P", Nv)
            labels <- switch(fun, us = {
                x <- t(lower.tri(matrix(nrow = n, ncol = n), 
                  diag = TRUE))
                paste("!", lab, ".", paste(out[[4]][col(x)[x]], 
                  out[[4]][row(x)[x]], sep = ":"), sep = "")
            }, chol = {
                x <- t(lower.tri(what.var, diag = TRUE))
                y <- col(x) - row(x) > k
                temp <- paste("!", out[[2]], ".", paste(out[[4]][col(x)[x]], 
                  out[[4]][row(x)[x]], sep = ":"), sep = "")
                temp[y[!lower.tri(y)]] <- paste(temp[y[!lower.tri(y)]], 
                  "<NotEstimated>", sep = "")
                temp
            }, cholc = {
                x <- t(lower.tri(what.var, diag = TRUE))
                y <- (col(x) - row(x) >= k) & row(x) > k
                temp <- paste("!", out[[2]], ".", paste(out[[4]][col(x)[x]], 
                  out[[4]][row(x)[x]], sep = ":"), sep = "")
                temp[y[!lower.tri(y)]] <- paste(temp[y[!lower.tri(y)]], 
                  "<NotEstimated>", sep = "")
                temp
            }, ante = {
                x <- t(lower.tri(what.var, diag = TRUE))
                y <- col(x) - row(x) > k
                temp <- paste("!", out[[2]], ".", paste(out[[4]][col(x)[x]], 
                  out[[4]][row(x)[x]], sep = ":"), sep = "")
                temp[y[!lower.tri(y)]] <- paste(temp[y[!lower.tri(y)]], 
                  "<NotEstimated>", sep = "")
                temp
            })
            out[[7]] <- labels
        }
        else {
            if (!IniFlag) 
                init <- what.var
            out[[5]] <- init
            con <- asreml.matchCon(init)
            if (length(con) == 0) 
                out[[6]] <- rep("P", Nv)
            else out[[6]] <- con
            out[[7]] <- paste("!", lab, ".", out[[4]], ".var", 
                sep = "")
        }
        out[[8]] <- rep(2, Nv)
        out[[13]] <- TRUE
    }
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    if (!nest) {
        out[[10]] <- c(n, 0, NA, 0)
        out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
        out[[12]] <- obj
    }
    if (length(out[[6]]) == 0) 
        out[[6]] <- ""
    if (length(out[[8]]) == 0) 
        out[[8]] <- NA
    oldClass(out) <- "asreml.special"
    out
}
asreml.spCall <-
function (sc) 
{
    sc$init <- sc$data <- sc$Rcov <- NULL
    return(substring(paste(deparse(sc), collapse = ""), 8))
}
asreml.Spcls <-
structure(list(Fun = c("con", "lin", "pow", "pol", "spl", "dev", 
"ped", "ide", "giv", "ma", "at", "and", "grp", "mbf", "id", "idv", 
"idh", "ar1", "ar1v", "ar1h", "ar2", "ar2v", "ar2h", "ar3", "ar3v", 
"ar3h", "sar", "sarv", "sarh", "sar2", "sar2v", "sar2h", "ma1", 
"ma1v", "ma1h", "ma2", "ma2v", "ma2h", "arma", "armav", "armah", 
"cor", "corv", "corh", "corb", "corbv", "corbh", "corg", "corgv", 
"corgh", "diag", "us", "sfa", "chol", "cholc", "ante", "exp", 
"expv", "exph", "iexp", "iexpv", "iexph", "aexp", "aexpv", "bexpv", 
"aexph", "gau", "gauv", "gauh", "igau", "igauv", "igauh", "agau", 
"agauv", "agauh", "ieuc", "ieucv", "ieuch", "isp", "sph", "sphv", 
"sphh", "ispv", "isph", "cir", "cirv", "cirh", "mtrn", "mtrnv", 
"mtrnh", "facv", "fa", "rr", "link", "str"), Code = c(1L, 1L, 
1L, -7L, -3L, -5L, -1L, -2L, -2L, -6L, -15L, -12L, -99L, -24L, 
0L, 0L, 8L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 8L, 9L, 10L, 11L, 11L, 12L, 6L, 6L, 6L, 6L, 6L, 
6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 
6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 13L, 15L, 
15L, -100L, -101L), Type = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 1L, 1L, 1L, 0L, 0L), ObjArgs = c(1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 0L, -1L)), .Names = c("Fun", 
"Code", "Type", "ObjArgs"), class = "data.frame", row.names = c(NA, 
-95L))
asreml.specialsOned <-
function (fun, obj, type, data, struc, dist, init = NA, p = 1) 
{
    out <- vector(mode = "list", length = 13)
    pwrpoints <- attr(data, "Control")$pwrpoints
    if (inherits(obj, "asreml.special")) {
        out[[1]] <- obj[[1]]
        out[[2]] <- obj[[2]]
        object <- out[[3]] <- obj[[3]]
        out[[10]] <- obj[[10]]
        out[[12]] <- obj[[12]]
    }
    else {
        object <- as.character(obj)
        out[[1]] <- fun
        out[[2]] <- object
        out[[3]] <- object
        out[[12]] <- object
    }
    if (length(pwrpoints[[object]]) > 0) 
        dist <- pwrpoints[[object]]
    else if (length(dist) == 1 && is.na(dist)) 
        dist <- data[[object]]
    else if (is.character(dist)) 
        dist <- eval(parse(text = dist))
    if (is.factor(dist)) 
        dist <- as.numeric(as.character(dist))
    if (asreml.getPdim(data, object) == 1) {
        PP <- asreml.getPpoints(data, object)
        pmat <- c(dist, PP)
        names(pmat) <- c(rep("DATA", length(dist)), rep("PP", 
            length(pmat) - length(dist)))
        pmat <- pmat[!duplicated(pmat)]
    }
    else if (asreml.getPdim(data, object) == 0) {
        PP <- list()
        pmat <- dist
    }
    out[[4]] <- unique(pmat)
    n <- length(out[[4]])
    ni <- switch(type, i = 1, v = 2, h = n + 1)
    IniFlag <- TRUE
    if (all(is.na(init))) {
        IniFlag <- FALSE
        init <- rep(0.1, ni)
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != ni) 
            stop(paste(fun, "Wrong number of initial values\n"))
    }
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) 
        out[[6]] <- switch(type, i = "U", v = c("U", "P"), h = c("U", 
            rep("P", n)))
    else out[[6]] <- con
    out[[7]] <- switch(type, i = paste("!", object, ".pow", sep = ""), 
        v = paste("!", object, c(".pow", ".var"), sep = ""), 
        h = c(paste("!", object, ".pow", sep = ""), paste("!", 
            object, ".", out[[4]], sep = "")), paste("!", object, 
            ".", out[[4]], sep = ""))
    out[[8]] <- switch(type, i = 3, v = c(3, 2), h = c(3, rep(2, 
        n)))
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- struc
    if (!inherits(obj, "asreml.special")) 
        out[[10]] <- c(n, 0, NA, 0)
    temp <- as.numeric(out[[4]])
    if (missing(dist)) {
        if (any(rep(temp, length.out = length(dist)) - asreml.ie(is.factor(dist), 
            as.numeric(as.character(dist)), dist)) != 0) 
            stop(paste("\nData not sorted in ", object, "order\n"))
    }
    out[[11]] <- list(coords = matrix(temp, nrow = 1))
    out[[13]] <- switch(type, i = FALSE, v = TRUE, h = TRUE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.specialsTwod <-
function (fun, type, xx, yy, data, struc, init = NA, ...) 
{
    call <- match.call(expand.dots = TRUE)
    if (is.null(call$p)) 
        p <- 1
    else p <- call$p
    if (is.null(call$Rcov)) 
        Rcov <- 0
    else Rcov <- call$Rcov
    if (is.null(call$dist)) 
        dist <- list()
    else dist <- call$dist
    if (is.element(type, c("mi", "mv", "mh"))) {
        phi <- call$phi
        nu <- call$nu
        delta <- call$delta
        alpha <- call$alpha
        lambda <- call$lambda
    }
    out <- vector(mode = "list", length = 13)
    pwrpoints <- attr(data, "Control")$pwrpoints
    if (inherits(xx, "asreml.special")) {
        out[[1]] <- xx[[1]]
        out[[2]] <- xx[[2]]
        obj <- out[[3]] <- xx[[3]]
        out[[10]] <- xx[[10]]
        out[[12]] <- xx[[12]]
        if (length(pwrpoints[[obj]]) > 0) 
            dist <- pwrpoints[[obj]]
        else if (length(dist) == 0) 
            stop(paste("No distances given to", fun))
        if (is.list(dist)) {
            if (length(dist) > 1) {
                xvals <- dist[[1]]
                yvals <- dist[[2]]
            }
            else {
                xvals <- dist[[1]]
                yvals <- asreml.ie(length(yy) == 1, rep(yy, length(xvals)), 
                  yy)
            }
        }
        else {
            xvals <- dist
            yvals <- asreml.ie(length(yy) == 1, rep(yy, length(xvals)), 
                yy)
        }
        if ((asreml.getPdim(data, obj) == 2)) 
            PP <- list(asreml.getPpoints(data, obj), asreml.getPpoints(data, 
                "y"))
        if ((asreml.getPdim(data, obj) == 1)) 
            PP <- list(asreml.getPpoints(data, obj), rep(yy, 
                length(asreml.getPpoints(data, obj))))
        else if ((asreml.getPdim(data, obj) == 0)) {
            PP <- list()
            pmat <- cbind(xvals, yvals)
        }
        if ((asreml.getPdim(data, obj) > 0)) {
            names(PP) <- c(obj, "yy")
            pmat <- rbind(cbind(xvals, yvals), cbind(PP[[obj]], 
                PP[["yy"]]))
            dimnames(pmat) <- list(c(rep("DATA", length(xvals)), 
                rep("PP", nrow(pmat) - length(xvals))), c(obj, 
                "yy"))
            pmat <- pmat[!duplicated(do.call("paste", as.data.frame(pmat))), 
                ]
        }
        if (Rcov == 0) 
            idx <- order(pmat[, 2], pmat[, 1])
        else idx <- seq(1, nrow(pmat))
        out[[4]] <- unique(paste(pmat[idx, 1], pmat[idx, 2], 
            sep = ","))
        n <- length(out[[4]])
    }
    else {
        out[[1]] <- fun
        obj <- paste(fun, xx, yy, sep = ".")
        out[[2]] <- obj
        out[[3]] <- obj
        if (is.na(match(obj, names(data)))) {
            out[[12]] <- c(xx, yy)
            out[[13]] <- TRUE
            if (fun == "mtrn") 
                out[[13]] <- FALSE
            names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
                "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", 
                "Argv", "isVariance")
            oldClass(out) <- "asreml.special"
            return(out)
        }
        if ((asreml.getPdim(data, xx) == 2) && (asreml.getPdim(data, 
            yy) == 2)) {
            PP <- list(asreml.getPpoints(data, xx), asreml.getPpoints(data, 
                yy))
            names(PP) <- c(xx, yy)
            pmat <- rbind(cbind(as.numeric(data[, xx]), as.numeric(data[, 
                yy])), cbind(PP[[xx]], PP[[yy]]))
            dimnames(pmat) <- list(c(rep("DATA", nrow(data)), 
                rep("PP", nrow(pmat) - nrow(data))), c(xx, yy))
            pmat <- pmat[!duplicated(do.call("paste", as.data.frame(pmat))), 
                ]
        }
        else if ((asreml.getPdim(data, xx) == 0) && (asreml.getPdim(data, 
            yy) == 0)) {
            PP <- list()
            pmat <- cbind(as.numeric(data[, xx]), as.numeric(data[, 
                yy]))
        }
        else {
            stop("Must supply predictpoints for both dimensions\n")
        }
        if (Rcov == 0) 
            idx <- order(pmat[, 2], pmat[, 1])
        else idx <- seq(1, nrow(pmat))
        out[[4]] <- unique(paste(pmat[idx, 1], pmat[idx, 2], 
            sep = ","))
        n <- length(out[[4]])
    }
    if (is.element(type, c("mi", "mv", "mh"))) {
        IniFlag <- switch(type, mi = rep(FALSE, 5), mv = {
            x <- rep(FALSE, 6)
            x[6] <- !attr(init, "missing")
            x
        }, mh = {
            x <- rep(FALSE, (5 + n))
            x[6:(5 + n)] <- rep(!attr(init, "missing"), n)
            x
        })
        initv <- switch(type, mi = vector(mode = "numeric", length = 5), 
            mv = vector(mode = "numeric", length = 6), mh = vector(mode = "numeric", 
                length = 5 + n))
        cons <- switch(type, mi = vector(mode = "character", 
            length = 5), mv = vector(mode = "character", length = 6), 
            mh = vector(mode = "character", length = 5 + n))
        if (!attr(phi, "missing")) {
            IniFlag[1] <- TRUE
            if (is.numeric(phi)) {
                initv[1] <- phi
                cons[1] <- "P"
            }
            else {
                initv[1] <- asreml.iniCon(phi)$value
                cons[1] <- asreml.iniCon(phi)$constraint
            }
        }
        else {
            initv[1] <- 0.1
            cons[1] <- "P"
        }
        if (!attr(nu, "missing")) {
            IniFlag[2] <- TRUE
            if (is.numeric(nu)) {
                initv[2] <- nu
                cons[2] <- "P"
            }
            else {
                initv[2] <- asreml.iniCon(nu)$value
                cons[2] <- asreml.iniCon(nu)$constraint
            }
        }
        else {
            initv[2] <- nu
            cons[2] <- "F"
        }
        if (!attr(delta, "missing")) {
            IniFlag[3] <- TRUE
            if (is.numeric(delta)) {
                initv[3] <- delta
                cons[3] <- "P"
            }
            else {
                initv[3] <- asreml.iniCon(delta)$value
                cons[3] <- asreml.iniCon(delta)$constraint
            }
        }
        else {
            initv[3] <- delta
            cons[3] <- "F"
        }
        if (!attr(alpha, "missing")) {
            IniFlag[4] <- TRUE
            if (is.numeric(alpha)) {
                initv[4] <- alpha
                cons[4] <- "P"
            }
            else {
                initv[4] <- asreml.iniCon(alpha)$value
                cons[4] <- asreml.iniCon(alpha)$constraint
            }
        }
        else {
            initv[4] <- alpha
            cons[4] <- "F"
        }
        if (!attr(lambda, "missing")) {
            IniFlag[5] <- TRUE
            if (is.numeric(lambda)) {
                initv[5] <- lambda
                cons[5] <- "F"
            }
            else {
                initv[5] <- asreml.iniCon(lambda)$value
                cons[5] <- "F"
            }
        }
        else {
            initv[5] <- lambda
            cons[5] <- "F"
        }
        out[[5]] <- initv
        out[[6]] <- cons
        out[[5]] <- switch(type, mi = initv, mv = {
            if (IniFlag[6] == FALSE) initv[6] <- 0.1 else {
                if (is.character(init)) init <- eval(parse(text = init))
                if (length(init) != 1) stop("mtrnv(): Wrong number of initial values\n")
                initv[6] <- init
            }
            initv
        }, mh = {
            if (IniFlag[6] == FALSE) initv[6:(n + 5)] <- rep(0.1, 
                n) else {
                if (is.character(init)) init <- eval(parse(text = init))
                if (length(init) != n) stop("mtrnh(): Wrong number of initial values\n")
                initv[6:(n + 5)] <- init
            }
            initv
        })
        out[[6]] <- switch(type, mi = cons, mv = {
            tmp <- asreml.matchCon(initv[6])
            if (length(tmp) == 0) cons[6] <- "P" else cons[6] <- tmp
            cons
        }, mh = {
            tmp <- asreml.matchCon(init[6:(n + 5)])
            if (length(tmp) == 0) cons[6:(n + 5)] <- rep("P", 
                n) else cons[6:(n + 5)] <- tmp
            cons
        })
    }
    else {
        ni <- switch(type, ii = 1, iv = 2, ih = n, ai = 2, av = 3, 
            ah = n + 2)
        IniFlag <- TRUE
        if (all(is.na(init))) {
            IniFlag <- FALSE
            init <- rep(0.1, ni)
        }
        else {
            if (is.character(init)) 
                init <- eval(parse(text = init))
            if (length(init) != ni) 
                stop(paste(fun, "Wrong number of initial values\n"))
        }
        out[[5]] <- init
        con <- asreml.matchCon(init)
        if (length(con) == 0) 
            out[[6]] <- switch(type, ii = "U", iv = c("U", "P"), 
                ih = c("U", rep("P", n)), ai = c("U", "U"), av = c("U", 
                  "U", "P"), ah = c("U", "U", rep("P", n)))
        else out[[6]] <- con
    }
    out[[7]] <- switch(type, ii = paste("!", "pow", sep = ""), 
        iv = paste("!", c("pow", "var"), sep = ""), ih = c(paste("!", 
            "pow", sep = ""), paste("!", out[[4]], sep = "")), 
        ai = paste("!", c(xx, yy), c(".pow", ".pow"), sep = ""), 
        av = paste("!", c(xx, yy, ""), c(".pow", ".pow", "var"), 
            sep = ""), ah = c(paste("!", c(xx, yy), c(".pow", 
            ".pow"), sep = ""), paste("!", out[[4]], sep = "")), 
        mi = paste("!", c("phi", "nu", "delta", "alpha", "lambda"), 
            sep = ""), mv = paste("!", c("phi", "nu", "delta", 
            "alpha", "lambda", "var"), sep = ""), mh = c(paste("!", 
            obj, c(".phi", ".nu", ".delta", ".alpha", ".lambda"), 
            sep = ""), paste("!", out[[4]], sep = "")))
    out[[8]] <- switch(type, ii = 3 * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), iv = c(3, 2) * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), ih = c(3, rep(2, n)) * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), ai = c(3, 3) * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), av = c(3, 3, 2) * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), ah = c(3, 3, rep(2, n)) * as.numeric(asreml.ie(IniFlag, 
        -1, 1)), mi = {
        x <- c(3, 3, 3, 3, 3)
        what <- as.numeric(IniFlag)
        what[what == 1] <- -1
        what[what == 0] <- 1
        x <- what * x
        x
    }, mv = {
        x <- c(3, 3, 3, 3, 3, 2)
        what <- as.numeric(IniFlag)
        what[what == 1] <- -1
        what[what == 0] <- 1
        x <- what * x
        x
    }, mh = {
        x <- c(3, 3, 3, 3, 3, rep(2, n))
        what <- as.numeric(IniFlag)
        what[what == 1] <- -1
        what[what == 0] <- 1
        x <- what * x
        x
    })
    out[[9]] <- struc
    if (!inherits(xx, "asreml.special")) 
        out[[10]] <- c(-5, 0, match(xx, names(data)), match(yy, 
            names(data)))
    out[[11]] <- list(coords = t(pmat[idx, ]), predict = PP, 
        points = PP)
    if (!inherits(xx, "asreml.special")) 
        out[[12]] <- c(obj, xx, yy)
    out[[13]] <- switch(type, ii = FALSE, iv = TRUE, ih = TRUE, 
        ai = FALSE, av = TRUE, ah = TRUE, mi = FALSE, mv = TRUE, 
        mh = TRUE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.sph <-
function (x, y, init = NA, data, Rcov) 
{
    fun <- "isp"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ii"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- do.call("asreml.specialsTwod", list(fun = fun, type = type, 
        xx = xx, yy = yy, data = data, struc = struc, init = init, 
        Rcov = Rcov))
    out
}
asreml.sphh <-
function (x, y, init = NA, data, ...) 
{
    fun <- "sphh"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "ih"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.sphv <-
function (x, y, init = NA, data, ...) 
{
    fun <- "sphv"
    xx <- as.character(substitute(x))
    yy <- as.character(substitute(y))
    type <- "iv"
    struc <- c(0, 0, 6, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    out <- asreml.specialsTwod(fun, type, xx, yy, data, struc, 
        init)
    out
}
asreml.spl <-
function (obj, k = 0, init = NA, data, ...) 
{
    if (mode(substitute(obj)) == "call" && inherits(obj, "asreml.special")) 
        stop("Argument to spl() must be a simple object (variate or factor)\n")
    out <- vector(mode = "list", length = 13)
    out[[1]] <- "spl"
    out[[2]] <- asreml.spCall(sys.call())
    obj <- as.character(substitute(obj))
    out[[3]] <- obj
    splXtras <- attr(data, "Control")$splXtras
    knots <- splXtras$knots
    kl <- 0
    kptr <- 0
    if (k > 0) {
        if (k < 3) 
            stop("Insufficient knots for spline (must have at least 3)\n")
        knots <- k
    }
    if (!is.null(splXtras$sPtr)) {
        which <- match(obj, dimnames(splXtras$sPtr)[[1]])
        if (!is.na(which)) {
            kptr <- splXtras$sPtr[which, 1]
            kl <- splXtras$sPtr[which, 2]
        }
    }
    if (kptr == 0) 
        kl <- k
    inter <- c(-3, kptr, 0, kl)
    ibv <- inter
    jbv <- inter
    if (!is.null(splXtras$IBV)) {
        which <- match(obj, dimnames(splXtras$IBV)[[2]])
        if (!is.na(which)) {
            ibv <- splXtras$IBV[, which]
            ibv[1] <- 1
            ibv[3] <- 1
        }
    }
    if (length(unique(data[[obj]])) < 3) 
        stop("Insufficient knots for spline (must have at least 3)\n")
    ncz <- inter[4] - 2
    N <- min(splXtras$splstp[3], length(unique(data[[obj]]))) + 
        ncz + splXtras$nsppoints + ibv[4]
    lx <- length(data[[obj]])
    one <- 1
    knotpoints <- vector(mode = "double", length = N)
    nux <- vector(mode = "integer", length = 1)
    storage.mode(ncz) <- "integer"
    fail <- 0
    temp <- .C("splinek", as.double(data[[obj]]), as.integer(one), 
        as.integer(lx), as.integer(inter), as.integer(ibv), as.integer(jbv), 
        as.double(splXtras$splstp), as.integer(knots), as.integer(splXtras$nsppoints), 
        as.double(splXtras$splinescale), as.double(splXtras$points), 
        as.integer(N), ncz = as.integer(ncz), nux = as.integer(nux), 
        knotpoints = as.double(knotpoints), fail = as.integer(fail), 
        NAOK = TRUE)
    if (temp$fail > 0) 
        stop()
    out[[4]] <- seq(1, temp$ncz)
    IniFlag <- TRUE
    if (missing(init)) {
        IniFlag <- FALSE
        init <- 0.1
    }
    else {
        if (is.character(init)) 
            init <- eval(parse(text = init))
        if (length(init) != 1) 
            stop("Wrong number of initial values\n")
    }
    out[[5]] <- init
    con <- asreml.matchCon(init)
    if (length(con) == 0) 
        out[[6]] <- "P"
    else out[[6]] <- con
    out[[7]] <- ""
    out[[8]] <- 2
    if (IniFlag) 
        out[[8]] <- -out[[8]]
    out[[9]] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out[[10]] <- c(-3, kptr, NA, kl)
    out[[11]] <- list(faconst = NULL, nuxPoints = temp$knotpoints[seq(1, 
        temp$nux)], coords = NULL)
    out[[12]] <- obj
    out[[13]] <- ifelse(IniFlag, TRUE, FALSE)
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.splinek <-
function (x, k = 0, points = NULL, predictpoints = NULL, Ppoints = asreml.control()$nsppoints, 
    splstp = asreml.control()$splstp, knots = asreml.control()$knots, 
    scale = asreml.control()$splinescale) 
{
    kl <- 0
    kptr <- 0
    kPptr <- 0
    if ((kl <- length(points)) > 0) 
        kptr <- 1
    else if (k > 0) 
        knots <- k
    kp <- 0
    if (ndim <- length(predictpoints[[1]]) > 0) {
        kPtr <- vector(mode = "numeric", length = ndim)
        for (i in seq(1, ndim)) {
            kPtr[i] <- length(points) + 1
            points <- c(points, predictpoints[[1]][[i]])
        }
        kp <- length(predictpoints[[1]][[1]])
    }
    N <- length(unique(x)) + max(k, knots) + max(Ppoints, kp)
    lx <- length(x)
    if (kptr == 0) 
        k1 <- k
    inter <- c(-3, kptr, 0, kl)
    if (kp > 0) {
        ibv <- c(1, kPptr[1], 1, kp)
        if (ndim == 2) 
            jbv <- c(1, kPptr[2], 1, kp)
    }
    else {
        ibv <- inter
        jbv <- inter
    }
    knotpoints <- vector(mode = "double", length = N)
    nux <- vector(mode = "integer", length = 1)
    ncz <- vector(mode = "integer", length = 1)
    lbin <- 1
    fail <- 0
    out <- .C("splinek", as.double(x), as.integer(lbin), as.integer(lx), 
        as.integer(inter), as.integer(ibv), as.integer(jbv), 
        as.double(splstp), as.integer(knots), as.integer(Ppoints), 
        as.double(scale), as.double(points), as.integer(N), ncz = as.integer(ncz), 
        nux = as.integer(nux), knotpoints = as.double(knotpoints), 
        fail = as.integer(fail), NAOK = TRUE)
    if (out$fail > 0) 
        stop()
    nux <- out$nux
    ncz <- out$ncz
    knotpoints <- out$knotpoints[seq(1, nux)]
    list(knotpoints = knotpoints, ncz = ncz, nux = nux)
}
asreml.splinez <-
function (x, k = 0, points = NULL, predictpoints = NULL, Ppoints = asreml.control()$nsppoints, 
    splstp = asreml.control()$splstp, knots = asreml.control()$knots, 
    scale = asreml.control()$splinescale) 
{
    kl <- 0
    kptr <- 0
    kPptr <- 0
    if ((kl <- length(points)) > 0) 
        kptr <- 1
    else if (k > 0) 
        knots <- k
    kp <- 0
    if (ndim <- length(predictpoints[[1]]) > 0) {
        kPtr <- vector(mode = "numeric", length = ndim)
        for (i in seq(1, ndim)) {
            kPtr[i] <- length(points) + 1
            points <- c(points, predictpoints[[1]][[i]])
        }
        kp <- length(predictpoints[[1]][[1]])
    }
    N <- length(unique(x)) + max(k, knots) + max(Ppoints, kp)
    lx <- length(x)
    if (kptr == 0) 
        k1 <- k
    inter <- c(-3, kptr, 0, kl)
    if (kp > 0) {
        ibv <- c(1, kPptr[1], 1, kp)
        if (ndim == 2) 
            jbv <- c(1, kPptr[2], 1, kp)
    }
    else {
        ibv <- inter
        jbv <- inter
    }
    z <- vector(mode = "double", length = N * N)
    zknots <- vector(mode = "double", length = N * N)
    knotpoints <- vector(mode = "double", length = N)
    designpoints <- vector(mode = "double", length = N)
    nux <- vector(mode = "integer", length = 1)
    ncz <- vector(mode = "integer", length = 1)
    out <- .C("splinez", as.double(x), as.integer(lx), as.integer(inter), 
        as.integer(ibv), as.integer(jbv), as.double(splstp), 
        as.integer(knots), as.integer(Ppoints), as.double(scale), 
        as.double(points), as.integer(N), z = as.double(z), ncz = as.integer(ncz), 
        nux = as.integer(nux), knotpoints = as.double(knotpoints), 
        designpoints = as.double(designpoints), zknots = as.double(zknots), 
        NAOK = TRUE)
    ncz <- out$ncz
    nux <- out$nux
    knotpoints <- out$knotpoints[seq(1, nux)]
    designpoints <- out$designpoints[seq(1, (ncz + 2))]
    z <- matrix(out$z[seq(1, nux * ncz)], nrow = nux, ncol = ncz, 
        byrow = TRUE)
    zknots <- matrix(out$zknots[seq(1, (ncz + 2) * ncz)], nrow = ncz + 
        2, ncol = ncz, byrow = TRUE)
    list(z = z, knotpoints = knotpoints, ncz = ncz, nux = nux, 
        designpoints = designpoints, zknots = zknots)
}
asreml.str <-
function (form, vmodel, data, ...) 
{
    terms.labels <- function(obj, keep.order) {
        if (is.list(obj)) 
            unlist(lapply(obj, function(x) attr(x, "term.labels")))
        else if (inherits(obj, "formula")) 
            unlist(lapply(terms.order(obj, keep.order), function(x) attr(x, 
                "term.labels")))
        else stop("Programming error")
    }
    out <- vector(mode = "list", length = 13)
    ttf <- terms(as.formula(form))
    ttm <- terms(as.formula(vmodel))
    form.lab <- terms.labels(form, TRUE)
    form.var <- lapply(form.lab, function(x) {
        dimnames(attr(terms(formula(paste("~", x))), "factors"))[[1]]
    })
    vmodel.lab <- dimnames(attr(ttm, "factors"))[[2]]
    if (length(vmodel.lab) > 1) 
        stop("Invalid vmodel formula in str()")
    vmodel.var <- dimnames(attr(ttm, "factors"))[[1]]
    ndim <- length(vmodel.var)
    Y <- lapply(vmodel.var, function(x, data) {
        X <- eval(asreml.Eval(x, data))
    }, data)
    out[[1]] <- "str"
    out[[2]] <- lapply(Y, function(x) x$Call)
    out[[3]] <- paste(unlist(lapply(form.var, function(x, data) {
        w <- sapply(x, function(z, data) {
            X <- eval(asreml.Eval(z, data))
            if (inherits(X, "asreml.special")) 
                X$Call
            else X
        }, data)
        paste(w, collapse = ":")
    }, data)), collapse = "+")
    lcheck <- lapply(Y, function(x) {
        if (x$Fun == "fa") {
            k <- asreml.fak(x$Lvls, x$Lab)
            seq(1:(length(x$Lvls) - k))
        }
        else x$Lvls
    })
    out[[4]] <- lapply(Y, function(x) x$Lvls)
    vsize <- prod(sapply(lcheck, function(x) length(x)))
    fsize <- sum(unlist(lapply(form.var, function(x, data) {
        w <- sapply(x, function(z, data) {
            X <- eval(asreml.Eval(z, data))
            if (inherits(X, "asreml.special")) {
                if (X$Fun == "and") NA else length(X$Lvls)
            } else ifelse(is.factor(data[, z]), length(levels(data[, 
                z])), 1)
        }, data)
        prod(w)
    }, data)), na.rm = TRUE)
    if (vsize != fsize) 
        stop(paste("Size of direct product (", vsize, ") does not conform with total size of included terms (", 
            fsize, ")", sep = ""))
    names(out[[4]]) <- sapply(Y, function(x) x$Fun)
    out[[5]] <- lapply(Y, function(x) x$Initial)
    out[[6]] <- lapply(Y, function(x) x$Con)
    names(out[[5]]) <- sapply(Y, function(x) x$Fun)
    names(out[[6]]) <- sapply(Y, function(x) x$Fun)
    out[[7]] <- lapply(Y, function(x) x$Lab)
    names(out[[7]]) <- sapply(Y, function(x) x$Fun)
    out[[8]] <- lapply(Y, function(x) x$Tgamma)
    names(out[[8]]) <- sapply(Y, function(x) x$Fun)
    out[[9]] <- lapply(Y, function(x) x$Struc)
    names(out[[9]]) <- sapply(Y, function(x) x$Fun)
    yy <- sapply(Y, function(x) x$Fun)
    out[[10]] <- vector(mode = "list", length = length(yy))
    names(out[[10]]) <- yy
    for (i in 1:length(yy)) out[[10]][[i]] <- c(-101, 0, NA, 
        0)
    out[[11]] <- list(faconst = NULL, nuxPoints = NULL, coords = NULL)
    out[[12]] <- paste(unlist(lapply(form.var, function(x, data) {
        w <- sapply(x, function(z, data) {
            X <- eval(asreml.Eval(z, data))
            if (inherits(X, "asreml.special") && X$Fun == "grp") 
                X$Obj
            else if (inherits(X, "asreml.special")) 
                X$Call
            else X
        }, data)
        paste(w, collapse = ":")
    }, data)), collapse = "+")
    out[[13]] <- {
        x <- unlist(lapply(Y, function(x) x$isVariance))
        if (any(is.na(x))) 
            NA
        else if (any(x)) 
            TRUE
        else FALSE
    }
    names(out) <- c("Fun", "Call", "Obj", "Lvls", "Initial", 
        "Con", "Lab", "Tgamma", "Struc", "Inter", "Coords", "Argv", 
        "isVariance")
    oldClass(out) <- "asreml.special"
    out
}
asreml.struc <-
function (data.frame, asr.inter, gcov.form, G.param, rcov.form, 
    R.param, fixgammas, Ainv) 
{
    match.term <- function(struc, inter) {
        which <- 0
        x <- strsplit(struc, ":", fixed = TRUE)
        for (i in seq(along = x)) {
            for (j in seq(along = inter)) {
                y <- strsplit(inter[j], ":", fixed = TRUE)[[1]]
                if (ifelse(any(is.na(match(x[[i]], y))), FALSE, 
                  TRUE)) {
                  which <- j
                  break
                }
            }
        }
        return(which)
    }
    arg.tyfact <- 10000
    neword <- asr.inter$neword
    nfactp <- length(neword)
    nfactd <- sum(as.numeric(asr.inter$fixedFlag))
    ngam <- nfactp
    facnam <- asr.inter$facnam
    minus12 <- (asr.inter$inter[1, ] == -12)
    gcov.mdl <- asreml.glist(gcov.form, data.frame)
    gcncode <- c(1, 1, 3, 4, 5, 6, 7)
    names(gcncode) <- c("P", "?", "U", "F", "C", "S", "B")
    lstruc <- 16
    struc <- matrix(NA, nrow = lstruc)
    struc.lab <- vector(mode = "character")
    tyfact <- rep(-1, nfactp)
    tyfact[asr.inter$fixedFlag | asr.inter$sparseFlag] <- 0
    tyfact[asr.inter$submds == 1] <- 0
    tyfact[minus12[neword]] <- 0
    gammas <- rep(0, ngam)
    names(gammas) <- rep("<NotEstimated>", ngam)
    tgamma <- rep(0, ngam)
    gmcstr <- rep(0, ngam)
    gpos <- length(gammas) + 1
    struc.point <- 1
    for (Gn in names(gcov.mdl)) {
        which <- match(Gn, facnam[neword])
        if (!is.na(which) && minus12[neword][which]) 
            next
        if (is.na(which)) {
            which <- match.term(Gn, facnam[neword])
            if (which == 0) {
                trms <- strsplit(Gn, " *\\+ *")[[1]]
                who <- seq(along = neword)[sapply(facnam[neword], 
                  function(x, y) {
                    as.logical(match.term(y, x))
                  }, trms)]
                if (length(who) > 0) {
                  if (match.term(trms, facnam[neword[who]]) == 
                    0) 
                    stop(paste("Cannot match", Gn, "in facnam\n"))
                  else which <- who[1]
                }
                else stop(paste("Cannot match", Gn, "in facnam\n"))
            }
        }
        gp <- match(Gn, names(G.param))
        if (is.na(gp)) 
            stop(paste("Cannot match", Gn, "in G.param\n"))
        strord <- length(gcov.mdl[[Gn]])
        strucCode <- sapply(gcov.mdl[[Gn]], function(x) {
            x$struc[3]
        })
        if (all(strucCode == 0) && !all(sapply(gcov.mdl[[Gn]], 
            function(x) {
                x$is.str
            }))) {
            initialValue <- unlist(sapply(G.param[[gp]], function(x) {
                x$initial
            }))
            gamnames <- unlist(sapply(G.param[[gp]], function(x) {
                names(x$initial)
            }))
            names(initialValue) <- gamnames
            ok <- !is.na(initialValue)
            gammas <- c(gammas, initialValue[ok])
            if (length(gammas) > (arg.tyfact - 1)) 
                stop(paste("A simple variance component has a position greater than", 
                  arg.tyfact, "in gammas"))
            tgamma <- c(tgamma, sapply(gcov.mdl[[Gn]], function(x) {
                x$tgamma
            })[ok])
            gmcstr <- c(gmcstr, gcncode[unlist(sapply(G.param[[gp]], 
                function(x) {
                  substr(x$con, 1, 1)
                })[ok])])
            tyfact[which] <- length(gammas)
            gpos <- length(gammas) + 1
        }
        else {
            for (i in seq(1, strord)) {
                buffer <- gcov.mdl[[Gn]][[i]]$struc
                buffer[1] <- gcov.mdl[[Gn]][[i]]$order
                if (sum(!is.na(G.param[[gp]][[i]]$initial)) == 
                  0) 
                  buffer[4] <- 0
                else buffer[4] <- gpos
                buffer[2] <- 0
                if (all(!is.na(G.param[[gp]][[i]]$initial))) {
                  buffer[2] <- length(G.param[[gp]][[i]]$initial)
                  gammas <- c(gammas, G.param[[gp]][[i]]$initial)
                  tgamma <- c(tgamma, gcov.mdl[[Gn]][[i]]$tgamma)
                  gmcstr <- c(gmcstr, gcncode[substr(G.param[[gp]][[i]]$con, 
                    1, 1)])
                }
                s2point <- buffer[13]
                buffer[13] <- ifelse(s2point == 0, 0, gpos + 
                  s2point)
                gpos <- gpos + buffer[2]
                if (buffer[3] == -7) 
                  buffer[3] <- -Ainv$locgiv[match(names(gcov.mdl[[Gn]][[i]]$struc)[3], 
                    names(Ainv$locgiv))]
                if (is.na(struc[1, 1])) 
                  struc <- cbind(buffer)
                else struc <- cbind(struc, buffer)
                if (buffer[3] == 6) 
                  Ainv <- asreml.ainvPow(Ainv, buffer[6], gcov.mdl[[Gn]][[i]]$coords, 
                    dim(struc)[[2]])
                struc.lab <- c(struc.lab, paste(Gn, names(gcov.mdl[[Gn]])[i], 
                  sep = "."))
            }
            tyfact[which] <- strord * arg.tyfact + struc.point
            struc.point <- struc.point + strord
        }
    }
    rcov.mdl <- asreml.rlist(rcov.form, data.frame)
    nsect <- length(rcov.mdl)
    nspat <- max(sapply(rcov.mdl, length)) - 1
    nrsect <- matrix(nrow = 3, ncol = nsect)
    for (n in seq(1, nsect)) {
        nrsect[1, n] <- rcov.mdl[[n]]$size$size
        ndim <- length(rcov.mdl[[n]]) - 1
        rp <- match(names(rcov.mdl)[n], names(R.param))
        if (is.na(rp)) 
            stop("\nR-model term not in initial value list")
        gammas <- c(gammas, R.param[[rp]]$variance$s2)
        tgamma <- c(tgamma, 1)
        gmcstr <- c(gmcstr, gcncode[substr(R.param[[rp]]$variance$con, 
            1, 1)])
        nrsect[2, n] <- gpos
        gpos <- gpos + 1
        nrsect[3, n] <- struc.point
        for (j in 1:ndim) {
            buffer <- rcov.mdl[[n]][[j]]$struc
            buffer[1] <- rcov.mdl[[n]][[j]]$length
            buffer[2] <- sum(!is.na(R.param[[rp]][[j]]$initial))
            buffer[4] <- ifelse(sum(!is.na(R.param[[rp]][[j]]$initial)) == 
                0, 0, gpos)
            if (all(!is.na(R.param[[rp]][[j]]$initial))) {
                gammas <- c(gammas, R.param[[rp]][[j]]$initial)
                tgamma <- c(tgamma, rcov.mdl[[n]][[j]]$tgamma)
                gmcstr <- c(gmcstr, gcncode[substr(R.param[[rp]][[j]]$con, 
                  1, 1)])
            }
            s2point <- buffer[13]
            buffer[13] <- ifelse(s2point == 0, 0, gpos + s2point)
            gpos <- gpos + buffer[2]
            if (buffer[3] == 0 && Ainv$locwts != 0) 
                buffer[3] <- -1
            if (is.na(struc[1, 1])) 
                struc <- cbind(buffer)
            else struc <- cbind(struc, buffer)
            if (buffer[3] == 6) 
                Ainv <- asreml.ainvPow(Ainv, buffer[6], rcov.mdl[[n]][[j]]$coords, 
                  dim(struc)[[2]])
            struc.lab <- c(struc.lab, paste(names(rcov.mdl[n]), 
                names(rcov.mdl[[n]])[j], sep = ":"))
        }
        struc.point <- struc.point + ndim
    }
    dimnames(struc) <- list(NULL, struc.lab)
    if (fixgammas) {
        gmcstr <- rep(4, length(gammas))
        names(gmcstr) <- rep("F", length(gammas))
    }
    names(tgamma) <- names(gammas)
    list(tyfact = tyfact, gammas = gammas, gmcstr = gmcstr, nfrstg = ngam + 
        1, tgamma = tgamma, struc = struc, nspat = nspat, nrsect = nrsect, 
        ainv = Ainv)
}
asreml.subForm <-
function (charVec) 
{
    ff <- unlist(strsplit(charVec, split = "\\+"))
    formVec <- vector(length = length(ff), mode = "list")
    for (i in seq(1, length(ff))) formVec[[i]] <- formula(paste("~", 
        ff[i]))
    formVec
}
asreml.terms <-
function (data, formula, facnam, varLevels, baseFac, inter, IBV, 
    JBV, coords, rndm, keep.order, locgiv) 
{
    rn <- function(str) {
        match(str, c("inter1", "inter2", "inter3", "inter4", 
            "flev", "nlev", "inform", "facord", "faconst"))
    }
    tt <- terms(formula, specials = asreml.Spcls$Fun, keep.order = keep.order)
    terms.fac <- attr(tt, "factors")
    terms.spc <- attr(tt, "specials")
    terms.var <- dimnames(terms.fac)[[1]]
    terms.lab <- dimnames(terms.fac)[[2]]
    nvar <- length(terms.var)
    which.spc <- unique(unlist(lapply(attr(tt, "specials"), function(x) x)))
    is.spc <- rep(FALSE, length(terms.var))
    is.spc[which.spc] <- TRUE
    names(is.spc) <- terms.var
    Y <- lapply(terms.var, function(x, is.spc, data) {
        if (is.spc[x]) 
            eval(asreml.Eval(x, data))
        else x
    }, is.spc, data)
    names(Y) <- terms.var
    for (z in seq(1, nvar)) {
        X <- Y[[terms.var[z]]]
        if (!inherits(X, "asreml.special")) {
            v <- match(X, names(data))
            if (!is.na(v)) {
                a <- max(1, length(levels(data[, v])))
                if (is.na(match(X, facnam))) {
                  inter <- cbind(inter, c(a, 0, v, 0, a, 1, -rndm, 
                    0, 0))
                  facnam <- c(facnam, X)
                  baseFac$facNum <- c(baseFac$facNum, length(facnam))
                  baseFac$baseFun <- c(baseFac$baseFun, X)
                  baseFac$baseObj <- c(baseFac$baseObj, X)
                }
                if (is.na(match(X, names(varLevels)))) {
                  if (a == 1) 
                    varLevels[[X]] <- X
                  else varLevels[[X]] <- levels(data[, v])
                }
            }
        }
        else if (X$Fun == "at") {
            varLevels[[terms.var[z]]] <- terms.var[z]
            v <- match(X$Obj, names(data))
            if (is.na(match(X$Obj, facnam))) {
                a <- levels(data[[X$Obj]])
                inter <- cbind(inter, c(length(a), 0, v, 0, length(a), 
                  1, -rndm, 0, 0))
                facnam <- c(facnam, X$Obj)
                baseFac$facNum <- c(baseFac$facNum, length(facnam))
                baseFac$baseFun <- c(baseFac$baseFun, X$Obj)
                baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
                varLevels[[X$Obj]] <- a
            }
            for (lvl in X$Lvls) {
                atnam <- paste("at(", X$Obj, ", ", lvl, ")", 
                  sep = "")
                where <- match(atnam, facnam)
                if (is.na(where)) {
                  a <- levels(data[[X$Obj]])
                  inter <- cbind(inter, c(as.numeric(X$Inter[1]), 
                    length(a), v, match(lvl, a), 1, 1, -rndm, 
                    0, 0))
                  facnam <- c(facnam, atnam)
                  baseFac$facNum <- c(baseFac$facNum, length(facnam))
                  baseFac$baseFun <- c(baseFac$baseFun, atnam)
                  baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
                  varLevels[[atnam]] <- paste(X$Obj, "_", lvl, 
                    sep = "")
                }
            }
        }
        else if (X$Fun == "and") {
            for (vv in X$Argv) {
                ee <- eval(asreml.Eval(vv, data))
                if (!inherits(ee, "asreml.special")) {
                  v <- match(vv, facnam)
                  if (is.na(v)) {
                    facnam <- c(facnam, vv)
                    inter <- cbind(inter, c(rep(0, 4), asreml.ie(is.factor(data[[vv]]), 
                      length(levels(data[[vv]])), 1), 1, 0, 0, 
                      0))
                    b <- dim(inter)[2]
                    inter[1, b] <- asreml.ie(is.factor(data[[vv]]), 
                      length(levels(data[[vv]])), 1)
                    inter[3, b] <- match(vv, names(data))
                    baseFac$facNum <- c(baseFac$facNum, length(facnam))
                    baseFac$baseFun <- c(baseFac$baseFun, vv)
                    baseFac$baseObj <- c(baseFac$baseObj, vv)
                    varLevels[[X$Obj]] <- asreml.ie(is.factor(data[[vv]]), 
                      levels(data[[vv]]), vv)
                  }
                }
                else if (is.element(ee$Fun, c("ped", "giv", "ide"))) {
                  if (!is.element(ee$Call, facnam)) {
                    inter <- cbind(inter, c(as.numeric(ee$Inter), 
                      length(ee$Lvls), 1, -rndm, 0, 0))
                    b <- dim(inter)[2]
                    inter[4, b] <- locgiv[match(ee$Obj, names(locgiv))]
                    facnam <- c(facnam, ee$Call)
                    baseFac$facNum <- c(baseFac$facNum, length(facnam))
                    baseFac$baseFun <- c(baseFac$baseFun, ee$Call)
                    baseFac$baseObj <- c(baseFac$baseObj, ee$Obj)
                    if (is.na(match(ee$Call, names(varLevels)))) 
                      varLevels[[ee$Call]] <- ee$Lvls
                  }
                }
                else if (is.element(ee$Fun, "fa")) {
                  if (!is.element(ee$Call, facnam)) {
                    inter <- cbind(inter, c(as.numeric(ee$Inter), 
                      length(ee$Lvls), 1, -rndm, 0, 0))
                    b <- dim(inter)[2]
                    inter[4, b] <- locgiv[match(ee$Obj, names(locgiv))]
                    facnam <- c(facnam, ee$Call)
                    baseFac$facNum <- c(baseFac$facNum, length(facnam))
                    baseFac$baseFun <- c(baseFac$baseFun, ee$Call)
                    baseFac$baseObj <- c(baseFac$baseObj, ee$Obj)
                    if (is.na(match(ee$Call, names(varLevels)))) 
                      varLevels[[ee$Call]] <- ee$Lvls
                  }
                }
                else stop(paste(ee$Fun, "not yet implemented in and()"))
            }
            if (length(X$Argv) == 1) {
                facnam <- c(facnam, X$Call)
                inter <- cbind(inter, c(rep(0, 4), 0, 0, -rndm, 
                  0, 0))
                b <- dim(inter)[2]
                inter[1:4, b] <- as.numeric(X$Inter)
                inter[2, b] <- match(X$Obj, facnam)
                baseFac$facNum <- c(baseFac$facNum, length(facnam))
                baseFac$baseFun <- c(baseFac$baseFun, X$Call)
                baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
            }
        }
        else if (X$Fun == "giv" | X$Fun == "ped") {
            inter <- cbind(inter, c(as.numeric(X$Inter), length(X$Lvls), 
                1, -rndm, 0, 0))
            b <- dim(inter)[2]
            inter[4, b] <- locgiv[match(X$Obj, names(locgiv))]
            facnam <- c(facnam, X$Call)
            baseFac$facNum <- c(baseFac$facNum, length(facnam))
            baseFac$baseFun <- c(baseFac$baseFun, X$Call)
            baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
            if (is.na(match(X$Call, names(varLevels)))) 
                varLevels[[X$Call]] <- X$Lvls
            at <- attr(X$Inter, "Hack")
            if (!is.null(at) && at == "grp") {
                inter[rn("nlev"), b] <- length(X$Lvls)
                if (is.na(inter[3, b])) 
                  inter[3, b] <- match(X$Argv[1], names(data))
            }
        }
        else if (X$Fun == "mbf") {
            inter <- cbind(inter, c(X$Inter, X$Inter[4], 1, -rndm, 
                0, 0))
            facnam <- c(facnam, X$Call)
            baseFac$facNum <- c(baseFac$facNum, length(facnam))
            baseFac$baseFun <- c(baseFac$baseFun, X$Call)
            baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
            if (is.na(match(X$Call, names(varLevels)))) 
                varLevels[[X$Call]] <- X$Lvls
            v <- X$Inter[3]
            a <- length(levels(as.factor(data[, v])))
            inter <- cbind(inter, c(a, 0, v, 0, a, 1, -rndm, 
                0, 0))
            facnam <- c(facnam, names(data)[v])
            baseFac$facNum <- c(baseFac$facNum, length(facnam))
            baseFac$baseFun <- c(baseFac$baseFun, names(data)[v])
            baseFac$baseObj <- c(baseFac$baseObj, names(data)[v])
            varLevels[[names(data)[v]]] <- levels(as.factor(data[, 
                v]))
        }
        else if (X$Fun == "grp") {
            v <- match(X$Argv[1], names(data))
            inter <- cbind(inter, c(rep(0, 4), length(X$Lvls), 
                0, -rndm, 0, 0))
            b <- dim(inter)[2]
            inter[1:4, b] <- as.numeric(X$Inter)
            if (is.na(inter[3, b])) 
                inter[3, b] <- v
            inter[rn("nlev"), b] <- length(X$Lvls)
            facnam <- c(facnam, X$Obj)
            baseFac$facNum <- c(baseFac$facNum, length(facnam))
            baseFac$baseFun <- c(baseFac$baseFun, X$Call)
            baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
            if (is.na(match(X$Call, names(varLevels)))) 
                varLevels[[X$Obj]] <- X$Lvls
            if (length(X$Coords$coords) > 0) 
                coords[[X$Obj]] <- X$Coords$coords
        }
        else {
            if (is.na(match(X$Call, facnam))) {
                v <- match(X$Obj, names(data))
                inter <- cbind(inter, c(rep(0, 4), length(X$Lvls), 
                  0, 0, 0, 0))
                b <- dim(inter)[2]
                inter[1:4, b] <- as.numeric(X$Inter)
                if (is.na(inter[3, b])) 
                  inter[3, b] <- v
                if (X$Fun == "pow") 
                  inter[rn("faconst"), b] <- X$Coords$faconst
                if (X$Fun == "spl" | X$Fun == "pol") {
                  if (!is.null(IBV)) {
                    which <- match(X$Obj, dimnames(IBV)[[2]])
                    if (!is.na(which)) 
                      IBV[3, which] <- v
                  }
                  coords$nuxPoints[[X$Obj]] <- X$Coords$nuxPoints
                }
                else if (!is.null(IBV)) {
                  which <- match(X$Obj, dimnames(IBV)[[2]])
                  if (!is.na(which)) 
                    IBV[3, which] <- X$Inter[3]
                  if (!is.null(JBV)) {
                    which <- match(X$Obj, dimnames(JBV)[[2]])
                    if (!is.na(which)) 
                      JBV[3, which] <- X$Inter[3]
                  }
                }
                inter[rn("nlev"), b] <- 1
                facnam <- c(facnam, X$Call)
                inter[rn("inform"), b] <- -rndm
                inter[rn("facord"), b] <- 0
                inter[rn("faconst"), b] <- 0
                baseFac$facNum <- c(baseFac$facNum, length(facnam))
                baseFac$baseFun <- c(baseFac$baseFun, X$Call)
                baseFac$baseObj <- c(baseFac$baseObj, X$Obj)
            }
            if (is.na(match(X$Call, names(varLevels)))) {
                varLevels[[X$Call]] <- X$Lvls
            }
        }
    }
    genFactors <- function(trm, lvl, facnam, last.fac, inter, 
        rndm, rn) {
        xL <- trm[1]
        xRlab <- paste(trm[-1], collapse = ":")
        xR <- trm[-1]
        ipA <- match(xL, facnam)
        if (is.na(ipB <- asreml.matchInteraction(xRlab, facnam))) {
            y <- genFactors(xR, lvl, facnam, last.fac, inter, 
                0, rn)
            facnam <- c(y$facnam, paste(xL, xRlab, sep = ":"))
            ipB <- asreml.matchInteraction(xRlab, facnam)
            last.fac <- y$last.fac
            inter <- cbind(y$inter, c(ipA, ipB, y$inter[rn("flev"), 
                ipA], y$inter[rn("flev"), ipB], prod(lvl[xL], 
                lvl[xR]), 0, rndm, last.fac, 0))
        }
        else {
            facnam <- c(facnam, paste(xL, xRlab, sep = ":"))
            inter <- cbind(inter, c(ipA, ipB, inter[rn("flev"), 
                ipA], inter[rn("flev"), ipB], prod(lvl[xL], lvl[xR]), 
                0, rndm, last.fac, 0))
        }
        last.fac <- last.fac + 1
        return(list(facnam = facnam, last.fac = last.fac, inter = inter))
    }
    tt <- terms.order(formula, keep.order)
    last.fac <- max(inter[rn("facord"), ])
    for (z in seq(along = terms.lab)) {
        Var <- dimnames(attr(tt[[z]], "factors"))[[1]]
        terms.use <- lapply(Var, function(z, y, data) {
            X <- y[[z]]
            if (!inherits(X, "asreml.special")) 
                return(list(Fun = "", Obj = X, Lvls = NA, Call = NULL))
            else if (X$Fun == "at") 
                return(list(Fun = "at", Obj = X$Obj, Lvls = X$Lvls, 
                  Call = paste("at(", X$Obj, ", ", X$Lvls, ")", 
                    sep = "")))
            else if (X$Fun == "grp") 
                return(list(Fun = X$Fun, Obj = X$Obj, Lvls = NA, 
                  Call = NULL))
            else if (X$Fun == "and") 
                return(list(Fun = X$Fun, Obj = {
                  if (length(X$Argv) == 1) X$Call else unlist(lapply(X$Argv, 
                    function(x, data) {
                      ee <- eval(asreml.Eval(x, data))
                      if (inherits(ee, "asreml.special")) return(ee$Call) else return(ee)
                    }, data))
                }, Lvls = X$Lvls, Call = X$Call, Inter = X$Inter))
            else return(list(Fun = X$Fun, Obj = X$Call, Lvls = NA, 
                Call = NULL))
        }, Y, data)
        Var <- sapply(terms.use, function(x) {
            x$Obj
        })
        if (length(Var) == 1) {
            inter[rn("inform"), match(Var, facnam)] <- rndm
            last.fac <- last.fac + 1
            inter[rn("facord"), match(Var, facnam)] <- last.fac
            next
        }
        which <- sapply(terms.use, function(x) {
            ifelse(x$Fun == "at", TRUE, FALSE)
        })
        ats <- seq(1, length(terms.use))[which]
        if (sum(which) > 1) {
            for (i in ats) {
                if (length(terms.use[[i]]$Lvls) > 1) 
                  stop("\nIllegal use of at()\n")
            }
            for (i in ats) Var[i] <- terms.use[[i]]$Call
            newLvl <- inter[rn("flev"), match(Var, facnam)]
            names(newLvl) <- Var
            gen.lst <- genFactors(Var, newLvl, facnam, last.fac, 
                inter, rndm, rn)
            facnam <- gen.lst$facnam
            inter <- gen.lst$inter
            last.fac <- gen.lst$last.fac
        }
        else if (sum(which) == 1) {
            for (lvl in terms.use[[ats]]$Call) {
                Var[which] <- lvl
                newLvl <- inter[rn("flev"), match(Var, facnam)]
                names(newLvl) <- Var
                gen.lst <- genFactors(Var, newLvl, facnam, last.fac, 
                  inter, rndm, rn)
                facnam <- gen.lst$facnam
                inter <- gen.lst$inter
                last.fac <- gen.lst$last.fac
            }
        }
        else {
            newLvl <- inter[rn("flev"), match(Var, facnam)]
            names(newLvl) <- Var
            gen.lst <- genFactors(Var, newLvl, facnam, last.fac, 
                inter, rndm, rn)
            facnam <- gen.lst$facnam
            inter <- gen.lst$inter
            last.fac <- gen.lst$last.fac
            if (length(terms.use) == 1 && terms.use[[1]]$Fun == 
                "and" && length(terms.use[[1]]$Obj) > 1) {
                facnam <- c(facnam, terms.use[[1]]$Call)
                inter <- cbind(inter, c(rep(0, 4), 0, 0, -rndm, 
                  0, 0))
                b <- dim(inter)[2]
                inter[1:4, b] <- as.numeric(terms.use[[1]]$Inter)
                oo <- paste(terms.use[[1]]$Obj, collapse = ":")
                inter[2, b] <- match(oo, facnam)
                inter[rn("inform"), match(oo, facnam)] <- 0
                baseFac$facNum <- c(baseFac$facNum, length(facnam))
                baseFac$baseFun <- c(baseFac$baseFun, terms.use[[1]]$Call)
                baseFac$baseObj <- c(baseFac$baseObj, oo)
                last.fac <- last.fac + 1
                inter[rn("facord"), b] <- last.fac
            }
        }
    }
    return(list(inter = inter, facnam = facnam, baseFac = baseFac, 
        varLevels = varLevels, IBV = IBV, JBV = JBV, coords = coords))
}
asreml.typeConvert <-
function (x) 
{
    if (is.numeric(x)) 
        return(x)
    if (asreml.Rsys) 
        return(type.convert(x, as.is = TRUE))
    else {
        ok = 1
        storage.mode(ok) <- "integer"
        n <- length(x)
        storage.mode(n) <- "integer"
        result <- .C("canBeNumeric", x = x, n, ok = ok)
        if (result$ok == 1) 
            return(as.numeric(x))
        else return(as.character(x))
    }
}
asreml.udu <-
function (gammas) 
{
    N <- length(gammas)
    which <- rep(1, N)
    which[asreml.grep("<NotEstimated>$", names(gammas))] <- 0
    n <- (sqrt(8 * N + 1) - 1)/2
    n <- floor(n + 0.5)
    q <- ((2 * n - 1) - sqrt((2 * n + 1)^2 - 8 * length(gammas[which == 
        1])))/2
    q <- floor(q + 0.5)
    U <- D <- matrix(0, nrow = n, ncol = n)
    xx <- 0 <= (col(U) - row(U)) & (col(U) - row(U)) <= q
    U[xx] <- gammas[which == 1]
    diag(D) <- diag(U)
    diag(U) <- 1
    V <- solve(U %*% D %*% t(U))
    V[!lower.tri(V)]
}
asreml.unpaste <-
function (char, sep) 
{
    if (asreml.Rsys) 
        return(strsplit(char, split = sep))
    else return(unpaste(char, sep = sep))
}
asreml.us <-
function (obj, init = NA, data, ...) 
{
    if (!(mode(substitute(obj)) == "call" && inherits(obj, "asreml.special"))) 
        obj <- substitute(obj)
    out <- asreml.spc("us", "var", numeric(0), "matrix(0.1,nrow=n,ncol=n)+diag(0.05,nrow=n)", 
        NA, obj, init, data, match.call()$Rcov)
    out[[9]] <- c(0, 0, 9, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0)
    out
}
asreml.us2fa <-
function (vcov, q = 1) 
{
    if (is.matrix(vcov)) 
        n <- nrow(vcov)
    else n <- (-1 + sqrt(8 * length(vcov) + 1))/2
    method <- 0
    ifault <- 0
    storage.mode(n) <- "integer"
    storage.mode(q) <- "integer"
    storage.mode(method) <- "integer"
    storage.mode(ifault) <- "integer"
    storage.mode(vcov) <- "double"
    fa <- .C("us2fa", fa = vcov, n, q, method, ifault = ifault)
    if (fa$ifault == -1) 
        stop("q must not be greater than (n-1)/2 \n")
    if (fa$ifault == 1) 
        stop("Diagonal elements of vcov must be positive \n")
    fa$fa[seq(1, n * q + n)]
}
asreml.utri2mat <-
function (vec, diag = TRUE, cor = TRUE) 
{
    makeMat <- function(v, m, diagonal) {
        mat <- matrix(nrow = m, ncol = m)
        mat[lower.tri(mat, diagonal)] <- v
        mat[!lower.tri(mat, TRUE)] <- t(mat)[!lower.tri(mat, 
            TRUE)]
        mat
    }
    if (is.logical(diag)) {
        n <- length(vec)
        dii <- numeric(0)
        if (diag) 
            m <- (sqrt(8 * n + 1) - 1)/2
        else {
            m <- (sqrt(8 * n + 1) + 1)/2
            dii <- rep(ifelse(cor, 1, 0), m)
        }
        mat <- makeMat(vec, m, diag)
    }
    else if (is.numeric(diag)) {
        nd <- length(diag)
        n <- length(vec) - nd
        dii <- {
            if (nd == 1) 
                rep(vec[diag], n)
            else vec[diag]
        }
        m <- (sqrt(8 * n + 1) + 1)/2
        mat <- makeMat(vec[1:n], m, FALSE)
    }
    if (length(dii)) 
        diag(mat) <- dii
    return(mat)
}
asreml.variogram <-
function (x, y, z, composite = TRUE, model = c("empirical"), 
    metric = c("euclidean", "manhattan"), angle = 0, angle.tol = 180, 
    nlag = 20, maxdist = 0.5, xlag = NA, lag.tol = 0.5, grid = TRUE) 
{
    model <- match.arg(model)
    metric <- match.arg(metric)
    if (missing(x)) 
        stop("x is a required argument\n")
    if (is.data.frame(x) || is.matrix(x)) {
        if (ncol(x) == 3) {
            y <- as.numeric(x[, 2])
            if (missing(z)) 
                z <- x[, 3]
            x <- as.numeric(x[, 1])
            n <- length(z)
        }
        else if (ncol(x) == 2) {
            if (missing(z)) {
                z <- x[, 2]
                x <- as.numeric(x[, 1])
                n <- length(z)
                y <- rep(1, n)
            }
            else {
                y <- as.numeric(x[, 2])
                x <- as.numeric(x[, 1])
                n <- length(z)
            }
        }
        else stop("x must be a data frame or matrix with 2 or 3 columns\n")
    }
    else {
        x <- as.numeric(x)
        n <- length(z)
        if (missing(y)) 
            y <- rep(1, n)
    }
    nlag <- min(n, nlag)
    ux <- sort(unique(x))
    uy <- sort(unique(y))
    nux <- length(ux)
    nuy <- length(uy)
    regular <- FALSE
    if (nux > 1 & nuy > 1) {
        if (nux * nuy == n) 
            regular <- TRUE
    }
    else if (nuy == 1) {
        if (sum(diff(diff(ux))) == 0) 
            regular <- TRUE
    }
    else if (nux == 1) {
        if (sum(diff(diff(uy))) == 0) 
            regular <- TRUE
    }
    regular <- (regular & grid)
    if (regular) {
        u <- 0
        type <- 1
        if (!composite) 
            type <- 2
        idx <- order(y, x)
        x <- x[idx]
        y <- y[idx]
        z <- z[idx]
    }
    else {
        if (asreml.Rsys) 
            u <- dist(cbind(x, y), method = metric)
        else u <- dist(cbind(x, y), metric = metric)
        if (any(u == 0)) 
            warning("Apparently duplicate points")
        type <- 0
        maxdist <- max(u) * maxdist
        if (is.na(xlag)) 
            xlag <- maxdist/nlag
        else nlag <- min(nlag, floor(maxdist/xlag))
        lag.tol <- xlag * lag.tol
        nu <- attr(u, "Size")
        full <- matrix(0, nu, nu)
        full[lower.tri(full)] <- u
        u <- full[lower.tri(full, diag = TRUE)]
    }
    nang <- length(angle)
    if (length(angle.tol) == 1) 
        angle.tol <- rep(angle.tol, nang)
    y <- as.double(y)
    x <- as.double(x)
    z <- as.double(z)
    u <- as.double(u)
    n <- as.integer(n)
    lbuf <- max(n, nlag) * nang
    vg <- double(length = lbuf)
    np <- integer(length = lbuf)
    distX <- double(length = lbuf)
    distY <- double(length = lbuf)
    direction <- double(length = lbuf)
    if (type == 2) {
        vg <- double(length = 2 * n)
        np <- integer(length = 2 * n)
        distX <- double(length = 2 * n)
        distY <- double(length = 2 * n)
    }
    type <- as.integer(type)
    maxdist <- as.double(maxdist)
    nlag <- as.integer(nlag)
    xlag <- as.double(xlag)
    lag.tol <- as.double(lag.tol)
    mode(nang) <- "integer"
    mode(angle) <- "double"
    mode(angle.tol) <- "double"
    nout <- n
    nout <- as.integer(nout)
    ifail <- 0
    ifail <- as.integer(ifail)
    vgram <- .C("variogram", n = n, x = x, y = y, z = z, distX = distX, 
        distY = distY, vg = vg, u = u, dist = distX, np = np, 
        type = type, nang, angle, angle.tol, angle = direction, 
        maxdist = maxdist, nlag = nlag, xlag = xlag, lag.tol = lag.tol, 
        nout = nout, NAOK = TRUE)
    if (type == 0) 
        return(data.frame(angle = vgram$angle[1:vgram$nout], 
            distance = vgram$dist[1:vgram$nout], gamma = vgram$vg[1:vgram$nout], 
            np = vgram$np[1:vgram$nout]))
    else return(data.frame(x = vgram$distX, y = vgram$distY, 
        gamma = vgram$vg, np = vgram$np))
}
coef.asreml <-
function (object, list = FALSE, pattern = character(0)) 
{
    if (list) {
        x <- c(object$coefficients$fixed, object$coefficients$random, 
            object$coefficients$sparse)
        labs <- attr(object$coefficients, "Terms")
        effects <- vector(mode = "list", length = length(unique(labs)))
        names(effects) <- unique(labs)
        for (la in unique(labs)) {
            which <- (la == labs)
            effects[[la]] <- matrix(x[which], ncol = 1)
            dimnames(effects[[la]]) <- list(names(x[which]), 
                "effect")
        }
    }
    else if (length(pattern) == 0) {
        effects <- vector(mode = "list", length = 3)
        names(effects) <- c("fixed", "random", "sparse")
        for (i in names(effects)) {
            if (length(object$coefficients[[i]]) > 0) {
                effects[[i]] <- matrix(object$coefficients[[i]], 
                  ncol = 1)
                dimnames(effects[[i]]) <- list(names(object$coefficients[[i]]), 
                  "effect")
            }
        }
    }
    else {
        trms <- strsplit(pattern, ":")[[1]]
        if (length(trms) == 1) 
            pattern <- paste(trms[1], ".*[^:]", sep = "")
        else pattern <- paste(trms[1], ".*", sep = "")
        if (length(trms) > 1) {
            for (i in 2:length(trms)) pattern <- paste(pattern, 
                ":", trms[i], ".*", sep = "")
        }
        effects <- lapply(object$coefficients, function(x, pattern) {
            which <- grep(pattern, names(x))
            x[which]
        }, pattern)
        effects <- c(effects$fixed, effects$random, effects$sparse)
        nn <- names(effects)
        effects <- matrix(effects, ncol = 1)
        dimnames(effects) <- list(nn, "effect")
    }
    effects
}
family.asreml <-
function (object) 
{
    obj <- object$call$family
    cl <- data.class(obj)[1]
    return(switch(cl, asreml.family = obj, `function` = obj(), 
        character = eval(parse(text = obj)), name = eval(parse(text = as.character(obj))), 
        call = eval(obj), `NULL` = asreml.gaussian(), stop("The object argument is invalid")))
}
fitted.asreml <-
function (object, type = c("response", "link")) 
{
    family <- asreml.getFamily(object$call$family)
    type <- match.arg(type)
    switch(type, link = object$linear.predictors, response = family$inverse(object$linear.predictors))
}
lrt.asreml <-
function (..., one.sided = TRUE) 
{
    call <- match.call(expand.dots = TRUE)
    if (!is.null(call$one.sided)) 
        call <- call[-match("one.sided", names(call))]
    if (length(call) < 3) 
        stop("Must specify at least 2 objects of class 'asreml'")
    args <- as.character(call[-1])
    call <- lapply(call[-1], function(x) eval(x))
    lapply(call, function(x) {
        if (!inherits(x, "asreml")) 
            stop("Arguments must be of class 'asreml'")
    })
    n <- length(call)
    fixed.labels <- lapply(call, function(x) {
        attr(terms(x$fixed.formula), "term.labels")
    })
    sparse.labels <- lapply(call, function(x) {
        attr(terms(x$sparse.formula), "term.labels")
    })
    mu <- sapply(call, function(x) {
        attr(terms(x$fixed.formula), "intercept")
    })
    if (!all(mu == mu[1])) 
        stop("fixed models must be identical - (Intercept)")
    if (!all(diff(sapply(fixed.labels, function(x) length(x))) == 
        0)) 
        stop("fixed models differ in length")
    if (all(sapply(fixed.labels, function(x) length(x)) > 0)) {
        for (i in 2:length(fixed.labels)) {
            if (!all(is.element(fixed.labels[[1]], fixed.labels[[i]]))) 
                stop("fixed models differ")
        }
    }
    if (!all(diff(sapply(sparse.labels, function(x) length(x))) == 
        0)) 
        stop("sparse models differ in length")
    if (all(sapply(sparse.labels, function(x) length(x)) > 0)) {
        for (i in 2:length(sparse.labels)) {
            if (!all(is.element(sparse.labels[[1]], sparse.labels[[i]]))) 
                stop("sparse models differ")
        }
    }
    qchisq.mixture <- function(prob, ntrait = 2, trace = F, maxiter = 10) {
        df <- 0:ntrait
        mixprobs <- dbinom(df, size = ntrait, prob = 0.5)
        cv <- qchisq(prob, df = ntrait)
        if (trace) 
            cat(" Starting value: ", cv, "\n")
        obj.fn <- function(cv = cv, df = cv, mixprobs = mixprobs) {
            obj <- ifelse(df[1] == 0, mixprobs[1] + sum(pchisq(cv, 
                df = df[-1]) * mixprobs[-1]), sum(pchisq(cv, 
                df = df) * mixprobs)) - prob
            d.obj <- ifelse(df[1] == 0, sum(dchisq(cv, df = df[-1]) * 
                mixprobs[-1]), sum(dchisq(cv, df = df) * mixprobs))
            list(f = obj, df = d.obj)
        }
        convergence <- F
        iteration <- 0
        while (!convergence & iteration <= maxiter) {
            obj <- obj.fn(cv = cv, df = df, mixprobs = mixprobs)
            corr <- (obj$f/obj$df)
            cv <- cv - corr
            iteration <- iteration + 1
            if (trace) 
                cat(" Iteration ", iteration, ": cv = ", cv, 
                  "\n")
            convergence <- abs(corr) < 1e-05
        }
        if (!convergence) 
            warning(" non-convergence: result incorrect\n")
        cv
    }
    pchisq.mixture <- function(x, ntrait = 2) {
        df <- 0:ntrait
        mixprobs <- dbinom(df, size = ntrait, prob = 0.5)
        p <- c()
        for (i in 1:length(x)) {
            p[i] <- sum(mixprobs * pchisq(x[i], df))
        }
        p
    }
    LL <- sapply(call, function(x) x$loglik)
    np <- sapply(call, function(x) {
        length(x$gammas) - sum(as.numeric(!(x$gammas.con == 1 | 
            x$gammas.con == 3)))
    })
    idx <- order(np)
    D <- df <- prob <- numeric(n - 1)
    models <- character(n - 1)
    for (i in seq(1, n - 1)) {
        ii <- idx[i]
        jj <- idx[i + 1]
        D[i] <- -2 * (LL[ii] - LL[jj])
        df[i] <- np[jj] - np[ii]
        prob[i] <- asreml.ie(one.sided, 1 - pchisq.mixture(D[i], 
            df[i]), 1 - pchisq(D[i], df[i]))
        models[i] <- paste(args[jj], "/", args[ii], sep = "")
    }
    tab <- matrix(nrow = n - 1, ncol = 3)
    tab[, 1] <- df
    tab[, 2] <- D
    tab[, 3] <- prob
    tab <- data.frame(tab)
    heading <- "Likelihood ratio test(s) assuming nested random models.\n"
    if (one.sided) 
        heading <- paste(heading, "Chisq probability adjusted using Stram & Lee, 1994.\n", 
            sep = "")
    attr(tab, "heading") <- heading
    attr(tab, "names") <- c("Df", "LR statistic", "Pr(Chisq)")
    attr(tab, "row.names") <- models
    oldClass(tab) <- c("anova", "data.frame")
    return(tab)
}
plot.asreml <-
function (object, formula = ~NULL, fun = NULL, res = "default", 
    spatial = "plot", npanels = NA, ...) 
{
    if (!inherits(object, "asreml")) 
        stop("\nObject must be of class asreml\n")
    call <- object$call
    control <- object$control
    drop.unused.levels <- control$drop.unused.levels
    fixed <- as.formula(object$fixed.formula)
    if (missingRcov <- is.null(call$rcov)) 
        rcov <- ~NULL
    else rcov <- as.formula(call$rcov)
    if (is.null(call$random)) 
        random <- ~NULL
    else random <- as.formula(object$random.formula)
    if (is.null(call$sparse)) 
        sparse <- ~NULL
    else sparse <- as.formula(object$sparse.formula)
    varNames <- names(eval(call$data, parent.frame()))
    groups <- control$group
    if (length(groups) > 0) {
        for (i in names(groups)) {
            if (is.numeric(groups[[i]])) 
                groups[[i]] <- varNames[groups[[i]]]
            else {
                if (any(is.na(match(groups[[i]], varNames)))) 
                  stop(paste("Object in group", names(groups[1]), 
                    "not in data."))
            }
        }
    }
    mbfList <- list()
    if (length(control$mbf) > 0) {
        mbfList <- asreml.getMbf(control$mbf, eval(call$data))
    }
    weights <- as.character(call$weights)
    offset <- as.character(call$offset)
    ran.order <- call$ran.order
    if (is.null(ran.order)) 
        ran.order <- "user"
    if (is.null(call$as.multivariate)) 
        as.multivariate <- NULL
    else as.multivariate <- call$as.multivariate
    if (is.null(call$na.method.Y)) 
        na.method.Y <- "include"
    else na.method.Y <- call$na.method.Y
    if (is.null(call$na.method.X)) 
        na.method.X <- "include"
    else na.method.X <- call$na.method.X
    form <- asreml.ModelFormula(fixed, random, sparse, rcov, 
        weights, offset, groups, mbfList$mbfAttr, ignore = c("mv", 
            "trait", "units"), varNames)
    model.y <- attr(form, "model.y")
    multivariate <- (length(model.y) > 1) || (!is.null(as.multivariate))
    if (is.null(call$subset)) 
        mf <- as.call(list(as.name("model.frame"), formula = form, 
            data = as.name(call$data), drop.unused.levels = drop.unused.levels, 
            na.action = "na.pass"))
    else mf <- as.call(list(as.name("model.frame"), formula = form, 
        data = as.name(call$data), subset = call$subset, drop.unused.levels = drop.unused.levels, 
        na.action = "na.pass"))
    mf[[1]] <- as.name("model.frame")
    attr(mf$formula, ".Environment") <- sys.frame(sys.parent())
    data <- eval(mf, sys.frame(sys.parent()))
    if (missingRcov) 
        rcov <- asreml.formula(parse(text = "~ units"))
    family <- asreml.getFamily(call$family)
    asr.glm <- asreml.glm(family)
    asrAttrib <- list(Control = control, GROUP = groups, MBF = mbfList)
    for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    data <- asreml.data(data, model.y, fixed, random, sparse, 
        rcov, asr.glm, asrAttrib, na.method.X, na.method.Y, as.multivariate, 
        weights, offset, ran.order)
    thrlevels <- attr(data, "thrlevels")
    keep <- attr(data, "keep")
    ntrt <- attr(data, "ntrt")
    model.y <- attr(data, "model.y")
    mvtrait <- attr(data, "mvtrait")
    args <- list(...)
    section <- asreml.sectionName(rcov, data)
    if (is.null(section)) 
        nsect <- 1
    else nsect <- length(unique(data[[section]]))
    if (!is.null(call$aom)) 
        aom <- as.logical(as.character(call$aom))
    else if (!is.null(object$control$aom)) 
        aom <- object$control$aom
    else aom <- asreml.control()$aom
    if (asr.glm$id == 1) {
        if (aom && (res == "default" || res == "stdCond")) {
            res <- "stdCond"
            caption <- "Studentised conditional residuals"
        }
        else {
            res <- "working"
            caption <- "Residuals"
        }
    }
    else {
        res <- "deviance"
        caption <- "Deviance residuals"
    }
    rr <- resid(object, type = res, spatial = spatial)
    fv <- fitted(object)
    rl <- rep(NA, length(keep))
    rl[keep] <- as.numeric(data$units)
    trellis.obj <- vector(length = 4, mode = "list")
    names(trellis.obj) <- c("histogram", "qq", "fitted", "rowlabels")
    if (((nsect == 1 & !multivariate) | aom | (asr.glm$id > 1)) & 
        length(formula[[2]]) == 0) {
        ff <- formula("~ rr")
        trellis.obj[["histogram"]] <- histogram(ff, xlab = caption, 
            ...)
        trellis.obj[["qq"]] <- qqmath(ff, ylab = caption, xlab = "Normal quantile", 
            ...)
        ff <- formula("rr ~ fv")
        trellis.obj[["fitted"]] <- xyplot(ff, ylab = caption, 
            xlab = "Fitted", panel = function(x, y) {
                panel.xyplot(x, y, type = c("p", "g"), cex = 0.5)
                panel.abline(h = 0)
            }, ...)
        ff <- formula("rr ~ rl")
        trellis.obj[["rowlabels"]] <- xyplot(ff, ylab = caption, 
            xlab = "Unit number", panel = function(x, y) {
                panel.xyplot(x, y, type = c("p", "g"), cex = 0.5)
                panel.abline(h = 0)
            }, ...)
        onepage <- TRUE
    }
    else if (nsect > 1 | multivariate) {
        if (multivariate) {
            if (nsect == 1) {
                if (is.na(npanels)) 
                  npanels <- length(unique(data[[mvtrait]]))
                trellis.obj[["histogram"]] <- histogram(~rr | 
                  data[[mvtrait]], xlab = caption, layout = c(0, 
                  npanels), ...)
                trellis.obj[["qq"]] <- qqmath(~rr | data[[mvtrait]], 
                  ylab = caption, xlab = "Normal quantile", layout = c(0, 
                    npanels), ...)
                trellis.obj[["fitted"]] <- xyplot(rr ~ fv | data[[mvtrait]], 
                  ylab = caption, xlab = "Fitted", layout = c(0, 
                    npanels), panel = function(x, y) {
                    panel.xyplot(x, y, type = c("p", "g"), cex = 0.5)
                    panel.abline(h = 0)
                  }, ...)
                trellis.obj[["rowlabels"]] <- xyplot(rr ~ rl | 
                  data[[mvtrait]], ylab = caption, xlab = "Unit number", 
                  layout = c(0, npanels), panel = function(x, 
                    y) {
                    panel.xyplot(x, y, type = c("p", "g"), cex = 0.5)
                    panel.abline(h = 0)
                  }, ...)
            }
            else {
                npanels <- length(unique(data[[mvtrait]])) * 
                  nsect
                stop("Multi-trait, multi-section plots not implemented")
            }
        }
        else {
            npanels <- nsect
            trellis.obj[["histogram"]] <- histogram(~rr | data[[section]], 
                xlab = caption, layout = c(0, npanels), ...)
            trellis.obj[["qq"]] <- qqmath(~rr | data[[section]], 
                ylab = caption, xlab = "Normal quantile", layout = c(0, 
                  npanels), ...)
            psym <- c(rep(c(letters, LETTERS), trunc(nsect/52)), 
                c(letters, LETTERS)[seq(1, nsect%%52)])
            trellis.obj[["fitted"]] <- xyplot(rr ~ fv, aspect = 0.75, 
                groups = data[[section]], pch = psym, key = list(text = list(paste(psym, 
                  levels(data[[section]]))), space = "right"), 
                panel = function(x, y, ...) {
                  panel.superpose(x, y, ..., type = c("p", "g"), 
                    cex = 0.5)
                  panel.abline(h = 0)
                }, ylab = caption, xlab = "Fitted", ...)
            trellis.obj[["rowlabels"]] <- xyplot(rr ~ rl, aspect = 0.75, 
                groups = data[[section]], pch = psym, key = list(text = list(paste(psym, 
                  levels(data[[section]]))), space = "right"), 
                panel = function(x, y, ...) {
                  panel.superpose(x, y, ..., type = c("p", "g"), 
                    cex = 0.5)
                  panel.abline(h = 0)
                }, ylab = caption, xlab = "Unit number", ...)
        }
        onepage <- FALSE
    }
    else if (length(formula[[2]]) > 0) {
        onepage <- FALSE
        trellis.obj <- list()
        if (is.null(fun)) 
            stop("'fun' must specify a trellis plot function")
        data <- as.list(c(as.list(data), . = list(object), object))
        args <- c(list(x = formula, data = data), args)
        trellis.obj[[1]] <- do.call(fun, args)
    }
    adt <- trellis.par.get("add.text")
    xlb <- trellis.par.get("par.xlab.text")
    ylb <- trellis.par.get("par.ylab.text")
    zlb <- trellis.par.get("par.zlab.text")
    axt <- trellis.par.get("axis.text")
    syx <- trellis.par.get("plot.symbol")
    trellis.par.set("add.text", list(cex = 0.75))
    trellis.par.set("par.xlab.text", list(cex = 0.75))
    trellis.par.set("par.ylab.text", list(cex = 0.75))
    trellis.par.set("par.zlab.text", list(cex = 0.75))
    trellis.par.set("axis.text", list(cex = 0.75))
    trellis.par.set("plot.symbol", list(cex = 0.6))
    if (onepage) {
        print(trellis.obj[[1]], position = c(0, 0.5, 0.5, 1), 
            more = TRUE)
        print(trellis.obj[[2]], position = c(0.5, 0.5, 1, 1), 
            more = TRUE)
        print(trellis.obj[[3]], position = c(0, 0, 0.5, 0.5), 
            more = TRUE)
        print(trellis.obj[[4]], position = c(0.5, 0, 1, 0.5))
    }
    else {
        grDevices::devAskNewPage(TRUE)
        lapply(trellis.obj, function(x) {
            print(x)
        })
        grDevices::devAskNewPage(FALSE)
    }
    trellis.par.set("add.text", adt)
    trellis.par.set("par.xlab.text", xlb)
    trellis.par.set("par.ylab.text", ylb)
    trellis.par.set("par.zlab.text", zlb)
    trellis.par.set("axis.text", axt)
    trellis.par.set("plot.symbol", syx)
    return(invisible(trellis.obj))
}
predict.asreml <-
function (object = NULL, classify = character(0), levels = list(), 
    present = list(), ignore = character(0), use = character(0), 
    except = character(0), only = character(0), associate = formula("~NULL"), 
    margin = logical(0), average = list(), as.average = list(), 
    vcov = logical(0), sed = logical(0), parallel = logical(0), 
    inrandom = logical(0), exrandom = logical(0), aliased = logical(0), 
    estimable = logical(0), xform = list(), evaluate = TRUE, 
    ...) 
{
    if (length(classify) == 0) 
        return(NULL)
    if (length(classify) > 1) 
        stop("The classify argument must be of length 1\n")
    if (mode(classify) == "list") 
        stop("The 'classify' argument must be a character string")
    Pnames <- asreml.unpaste(classify, ":")[[1]]
    if (length(Pnames) > 1) {
        if (mode(levels) != "list") 
            stop("\nlevels argument must be of mode list\n")
        if (any(is.na(match(names(levels), Pnames)))) 
            stop(paste("\n'levels' must be a list with component names a subset of", 
                paste(Pnames, collapse = ",")))
    }
    else if (!is.list(levels)) {
        levels <- list(x = levels)
        names(levels) <- Pnames
    }
    if (length(present) > 0) {
        if (mode(present) == "list") {
            if (length(present) == 1 && names(present) == classify) 
                present <- present[[1]]
        }
    }
    if (is.data.frame(average)) 
        average <- as.list(average)
    if (is.data.frame(as.average)) 
        as.average <- as.list(as.average)
    if (!inherits(associate, "formula")) 
        stop("'associate' must be a formula\n")
    if (length(associate[[2]]) == 0) 
        associate <- list()
    else {
        ff <- asreml.getOp(associate[[2]], "+")
        if (is.null(ff)) {
            ffL <- dimnames(attr(terms(associate), "factors"))[[1]]
            ffR <- character(0)
        }
        else {
            ffL <- dimnames(attr(terms(formula(paste("~", asreml.deparse(ff[[2]])))), 
                "factors"))[[1]]
            ffR <- dimnames(attr(terms(formula(paste("~", asreml.deparse(ff[[3]])))), 
                "factors"))[[1]]
        }
        associate <- list(ffL, ffR)
    }
    if (length(ignore) > 0 && mode(ignore) == "list") 
        ignore <- ignore[[1]]
    if (length(use) > 0 && mode(use) == "list") 
        use <- use[[1]]
    if (length(except) > 0 && mode(except) == "list") 
        except <- except[[1]]
    if (length(only) > 0 && mode(only) == "list") 
        only <- only[[1]]
    if (length(margin) > 0 && mode(margin) == "list") 
        margin <- margin[[1]]
    if (length(vcov) > 0 && mode(vcov) == "list") 
        vcov <- vcov[[1]]
    if (length(sed) > 0 && mode(sed) == "list") 
        sed <- sed[[1]]
    if (length(parallel) && mode(parallel) != "logical") 
        stop("'parallel' must be a logical variable")
    if (length(parallel) == 0) 
        parallel <- FALSE
    if (length(inrandom) > 0 && inrandom == "list") 
        inrandom <- inrandom[[1]]
    if (length(exrandom) > 0 && mode(exrandom) == "list") 
        exrandom <- exrandom[[1]]
    if (length(aliased) > 0 && mode(aliased) == "list") 
        aliased <- aliased[[1]]
    if (length(estimable) > 0 && mode(estimable) == "list") 
        estimable <- estimable[[1]]
    prd <- list(classify = classify, levels = levels, present = present, 
        ignore = ignore, use = use, except = except, only = only, 
        associate = associate, margin = asreml.ie(length(margin) == 
            0, FALSE, margin), average = average, as.average = as.average, 
        vcov = asreml.ie(length(vcov) == 0, FALSE, vcov), sed = asreml.ie(length(sed) == 
            0, FALSE, sed), parallel = parallel, inrandom = asreml.ie(length(inrandom) == 
            0, TRUE, inrandom), exrandom = exrandom, aliased = asreml.ie(length(aliased) == 
            0, FALSE, aliased), estimable = asreml.ie(length(estimable) == 
            0, FALSE, estimable), xform = xform)
    if (is.null(object) | !evaluate) 
        return(prd)
    update(object, predict = prd, ...)
}
print.asreml.family <-
function (x, ...) 
{
    fm <- x$family
    dist <- fm[1]
    lnk <- substring(fm[2], 1, match(":", substring(fm[2], 1:nchar(fm[2]), 
        1:nchar(fm[2]))) - 1)
    cat("Family:", dist, "\n")
    cat("Link function:", lnk, "\n")
}
print.asremlPredict <-
function (x, digits = .Options$digits, quote = FALSE, ...) 
{
    heading <- attr(x, "heading")
    if (!is.null(heading)) 
        cat("\nNotes:", heading, sep = "\n")
    attr(x, "heading") <- NULL
    d <- dim(x)
    for (i in 1:d[2]) {
        xx <- x[[i]]
        if (!length(levels(xx)) && is.numeric(xx)) {
            xx <- format(zapsmall(xx, digits))
            x[[i]] <- xx
        }
    }
    invisible(NextMethod("print"))
}
residuals.asreml <-
function (object, type = c("stdCond", "working", "deviance", 
    "pearson", "response"), spatial = c("trend", "plot")) 
{
    type <- match.arg(type)
    if (type == "stdCond") {
        if (length(object$aom) > 0) 
            rw <- object$aom$R[, "stdCondRes"]
        else {
            rw <- object$residuals
            type <- "deviance"
        }
    }
    else rw <- object$residuals
    if (casefold(names(object$distribution)) == "gaussian") {
        spatial <- match.arg(spatial)
        rw <- switch(spatial, trend = rw, plot = {
            type <- "working"
            if (is.null(object$call$random)) rw else {
                units <- coef(object)$random
                i <- asreml.grep("units_.*", dimnames(units)[[1]])
                if (length(i) > 0) {
                  ncol <- length(rw)/length(i)
                  (rw + as.vector(matrix(units[i], ncol = ncol, 
                    byrow = FALSE) %*% rep(1, ncol)))
                } else rw
            }
        })
    }
    if (is.null(object$call$data)) 
        y <- eval(object$call$fixed[[2]])
    else if (mode(object$call$data) == "name" || mode(object$call$data) == 
        "character") {
        data <- eval(as.name(object$call$data), sys.frame(sys.parent(2)))
        y <- eval(eval(object$call$fixed)[[2]], envir = data)
    }
    else y <- eval(eval(object$call$fixed)[[2]], envir = object$call$data)
    if (is.matrix(y)) 
        N <- nrow(y)
    else N <- length(y)
    if (is.null(object$call$subset)) 
        subset <- rep(TRUE, N)
    else if (is.null(object$call$data)) 
        subset <- as.vector(eval(object$call$subset))
    else if (asreml.Rsys) 
        subset <- as.vector(eval(object$call$subset, envir = data))
    else subset <- as.vector(eval(object$call$subset, local = data))
    if (is.matrix(y)) 
        y <- as.vector(t(y[subset, ]))
    else y <- y[subset]
    w <- object$prior.weights
    if (is.null(w)) 
        w <- rep(1, length(rw))
    family <- asreml.getFamily(object$call$family)
    Fv <- object$linear.predictors
    mu <- family$inverse(Fv)
    rr <- switch(type, stdCond = rw, working = rw, pearson = sqrt(w) * 
        rw/family$deriv(mu)/sqrt(family$variance(mu)), deviance = {
        where <- is.na(y)
        y[where] <- 0
        rr <- family$deviance(mu, y, w, residuals = TRUE, phi = family$phi)
        rr[where] <- NA
        rr
    }, response = rw/family$deriv(mu))
    rr
}
summary.asreml <-
function (object, nice = FALSE, all = FALSE) 
{
    safeSqrt <- function(x) {
        x[!is.na(x) & x > 0] <- sqrt(x[!is.na(x) & x > 0])
        x[!is.na(x) & x < 0] <- NA
        x
    }
    ng <- length(object$gammas)
    s2 <- object$sigma2
    nsect <- length(object$R.param)
    xx <- seq(1, ng)[object$gammas.type == 1]
    x <- xx[length(xx)]
    ntrt <- 1
    if (ntrt + nsect == 2) 
        convrt <- object$gammas[x] == 1 && object$gammas.con[x] != 
            4
    else convrt <- FALSE
    y <- asreml.ltri2mat(object$ai)
    v <- vector(mode = "numeric", length = ng)
    if (convrt) {
        for (i in seq(1, ng)) {
            if (object$gammas.type[i]%%2 != 0) 
                v[i] <- y[i, i]
            else v[i] <- (object$gammas[i]^2) * y[x, x] + 2 * 
                s2 * object$gammas[i] * y[x, i] + (s2^2) * y[i, 
                i]
        }
        scale <- rep(s2, ng)
        scale[object$gammas.type == 1 | object$gammas.type == 
            3 | object$gammas.type == 5] <- 1
        scale[x] <- s2
    }
    else {
        v <- diag(y)
        scale <- rep(1, ng)
    }
    v[v == 0] <- NA
    v <- safeSqrt(v)
    comp <- object$gammas * scale
    tval <- comp/v
    varcomp <- matrix(nrow = ng, ncol = 5)
    varcomp <- cbind(object$gammas, comp, v, tval)
    dimnames(varcomp) <- list(names(object$gammas), c("gamma", 
        "component", "std.error", "z ratio"))
    constraint <- matrix(asreml.guzpfx(object$gammas.con), ncol = 1)
    dimnames(constraint) <- list(names(object$gammas), "constraint")
    ndense <- length(object$coefficients$fixed)
    fix.nam <- names(object$coefficients$fixed)
    if (is.null(object$Cfixed)) 
        cinv <- NULL
    else {
        cinv <- asreml.ltri2mat(object$Cfixed)
        dimnames(cinv) <- list(fix.nam, fix.nam)
    }
    if (all) {
        if (convrt) 
            ss <- s2
        else ss <- 1
        fixed <- matrix(nrow = ndense, ncol = 3)
        fixed[, 1] <- object$coefficients$fixed
        fixed[, 2] <- safeSqrt(ss * object$vcoeff$fixed)
        fixed[, 2][fixed[, 2] == 0] <- NA
        fixed[, 3] <- fixed[, 1]/fixed[, 2]
        dimnames(fixed) <- list(fix.nam, c("solution", "std error", 
            "z ratio"))
        if (length(object$coefficients$random) > 0) {
            random <- matrix(nrow = length(object$coefficients$random), 
                ncol = 3)
            random[, 1] <- object$coefficients$random
            random[, 2] <- safeSqrt(ss * object$vcoeff$random)
            random[, 2][random[, 2] == 0] <- NA
            random[, 3] <- random[, 1]/random[, 2]
            dimnames(random) <- list(names(object$coefficients$random), 
                c("solution", "std error", "z ratio"))
        }
        else random <- NULL
        if (length(object$coefficients$sparse) > 0) {
            sparse <- matrix(nrow = length(object$coefficients$sparse), 
                ncol = 3)
            sparse[, 1] <- object$coefficients$sparse
            sparse[, 2] <- safeSqrt(ss * object$vcoeff$sparse)
            sparse[, 2][sparse[, 2] == 0] <- NA
            sparse[, 3] <- sparse[, 1]/sparse[, 2]
            dimnames(sparse) <- list(names(object$coefficients$sparse), 
                c("solution", "std error", "z ratio"))
        }
        else sparse <- NULL
    }
    else fixed <- random <- sparse <- NULL
    sum.object <- vector(mode = "list")
    sum.object$call <- object$call
    sum.object$loglik <- object$loglik
    sum.object$nedf <- object$nedf
    sum.object$sigma <- safeSqrt(s2)
    sum.object$varcomp <- cbind(data.frame(varcomp), constraint)
    if (nice) 
        sum.object$nice <- asreml.niceGammas(object)
    if (object$deviance != 0) {
        sum.object$deviance <- object$deviance
        sum.object$heterogeneity <- object$deviance/object$nedf
    }
    sum.object$Cfixed <- cinv
    if (all) {
        sum.object$distribution <- names(object$distribution)
        sum.object$link <- names(object$link)
        sum.object$coef.fixed <- fixed
        sum.object$coef.random <- random
        sum.object$coef.sparse <- sparse
    }
    class(sum.object) <- "summary.asreml"
    sum.object
}
svc.asreml <-
function (object, order = "user") 
{
    ng <- length(object$gammas)
    s2 <- object$sigma2
    x <- seq(1, ng)[object$gammas.type == 1]
    ai <- matrix(0, nrow = ng, ncol = ng)
    convrt <- as.numeric(ifelse(length(x) == 1, object$gammas[x] == 
        1 & object$gammas.con[x] != 4, FALSE))
    object$ai <- zapsmall(object$ai)
    object$gammas[object$gammas.con == 7] <- 0
    storage.mode(object$ai) <- "double"
    storage.mode(object$gammas) <- "double"
    storage.mode(object$gammas.type) <- "integer"
    storage.mode(object$sigma2) <- "double"
    storage.mode(convrt) <- "integer"
    storage.mode(ai) <- "double"
    .Call("aimat", object$ai, object$gammas, object$gammas.type, 
        object$sigma2, convrt, ai)
    if (is.character(order) && length(order) == 1) {
        idx <- switch(order, user = seq(along = object$gammas), 
            R = {
                new <- c(attr(terms(formula(paste("~", paste(names(object$gammas)[-x], 
                  collapse = "+"))), keep.order = FALSE), "term.labels"), 
                  names(object$gammas)[x])
                match(new, names(object$gammas))
            }, noeff = {
                v <- 1:ng
                v[-x] <- (v[-x])[order(object$noeff[match(names(object$gammas)[-x], 
                  names(object$noeff))])]
                match(1:ng, v)
            }, stop("odrer must be one of 'user','R','noeff'"))
    }
    if (is.character(order) && length(order) > 1) 
        idx <- match(order, names(object$gammas))
    if (is.numeric(order)) 
        idx <- order
    if (length(idx) > ng || max(idx) > ng) 
        stop("'order' too long or elements out of range")
    if (convrt == 1) {
        scale <- rep(s2, ng)
        scale[object$gammas.type == 3 | object$gammas.type == 
            5] <- 1
    }
    else scale <- rep(1, ng)
    comp <- (object$gammas * scale)[idx]
    gg <- names(object$gammas)[idx]
    ai <- ai[idx, idx]
    which <- !apply(ai, 1, function(x) {
        all(abs(x) < 1e-10)
    })
    gg <- gg[which]
    comp <- comp[which]
    ff <- ai[which, which]
    F <- solve(ff)
    U <- chol(F)
    Dc <- diag(1/U[, ncol(U)])
    Uc <- Dc %*% U
    eta <- Uc %*% comp
    df <- 2 * (eta/diag(Dc))^2
    svtable <- cbind(df, eta, Uc)
    dimnames(svtable) <- list(gg, c("df", "variance", gg))
    zapsmall(svtable)
}
tr.asreml <-
function (obj, gammas = seq(along = obj$gammas), iter = seq(1, 
    length(obj$monitor) - 1), loglik = FALSE, S2 = FALSE, Rvariance = FALSE) 
{
    if (!inherits(obj, "asreml")) 
        stop("\nObject must be of class asreml\n")
    which <- names(obj$gammas)[gammas]
    if (!is.na(what <- match("R!variance", which)) && Rvariance == 
        FALSE) 
        which <- which[-what]
    frame()
    ask <- ((howMany <- length(which) + sum(as.numeric(c(loglik, 
        S2)))) > 8)
    pp <- par(mfrow = c(min(8, howMany), 1), mar = c(1, 4, 2, 
        2), cex = 0.6, ask = ask)
    if (loglik) {
        plot(iter, obj$monitor["loglik", iter], xaxt = "n", bty = "l", 
            xlab = "", ylab = "", type = "n", main = "loglik")
        text(iter, obj$monitor["loglik", iter], labels = iter)
    }
    if (S2) {
        plot(iter, obj$monitor["S2", iter], xaxt = "n", bty = "l", 
            xlab = "", ylab = "", type = "n", main = "S2")
        text(iter, obj$monitor["S2", iter], labels = iter)
    }
    for (i in which) {
        plot(iter, obj$monitor[i, iter], xaxt = "n", bty = "l", 
            xlab = "", ylab = "", type = "n", main = i)
        text(iter, obj$monitor[i, iter], labels = iter)
    }
    par(pp)
    invisible()
}
update.asreml <-
function (object, fixed., random., sparse., rcov., ..., keep.order = TRUE, 
    evaluate = TRUE) 
{
    my.update.formula <- function(old, new, keep.order = TRUE, 
        ...) {
        env <- environment(as.formula(old))
        if (R.Version()$major < 3) 
            tmp <- .Internal(update.formula(as.formula(old), 
                as.formula(new)))
        else tmp <- update.formula(as.formula(old), as.formula(new))
        out <- formula(terms.formula(tmp, simplify = TRUE, keep.order = keep.order))
        environment(out) <- env
        return(out)
    }
    if (is.null(newcall <- object$call) && is.null(newcall <- attr(object, 
        "call"))) 
        stop("need an object with call component or attribute")
    tempcall <- list(...)
    if (length(newcall$fixed)) 
        newcall$fixed <- eval(newcall$fixed)
    if (length(newcall$random)) 
        newcall$random <- eval(newcall$random)
    if (length(newcall$sparse)) 
        newcall$sparse <- eval(newcall$sparse)
    if (length(newcall$rcov)) 
        newcall$rcov <- eval(newcall$rcov)
    if (!missing(fixed.)) 
        newcall$fixed <- my.update.formula(as.formula(newcall$fixed), 
            fixed., keep.order = keep.order)
    if (!missing(random.)) 
        newcall$random <- {
            if (length(newcall$random)) 
                my.update.formula(as.formula(newcall$random), 
                  random., keep.order = keep.order)
            else random.
        }
    if (!missing(sparse.)) 
        newcall$sparse <- {
            if (length(newcall$sparse)) 
                my.update.formula(as.formula(newcall$sparse), 
                  sparse., keep.order = keep.order)
            else sparse.
        }
    if (!missing(rcov.)) 
        newcall$rcov <- {
            if (length(newcall$rcov)) 
                my.update.formula(as.formula(newcall$rcov), rcov., 
                  keep.order = keep.order)
            else rcov.
        }
    if (length(tempcall)) {
        what <- !is.na(match(names(tempcall), names(newcall)))
        for (z in names(tempcall)[what]) newcall[[z]] <- tempcall[[z]]
        if (any(!what)) {
            newcall <- c(as.list(newcall), tempcall[!what])
            newcall <- as.call(newcall)
        }
    }
    if (length(newcall$random) && length(attr(terms(as.formula(newcall$random)), 
        "factors")) == 0) 
        newcall$random <- NULL
    if (length(newcall$sparse) && length(attr(terms(as.formula(newcall$sparse)), 
        "factors")) == 0) 
        newcall$sparse <- NULL
    if (length(newcall$rcov) && length(attr(terms(as.formula(newcall$rcov)), 
        "factors")) == 0) 
        newcall$rcov <- NULL
    newcall$R.param <- call("$", as.name(deparse(substitute(object))), 
        "R.param")
    newcall$G.param <- call("$", as.name(deparse(substitute(object))), 
        "G.param")
    if (evaluate) 
        eval(newcall, sys.parent())
    else newcall
}
variogram.asreml <-
function (object, formula = ~NULL, composite = TRUE, model = c("empirical"), 
    metric = c("euclidean", "manhattan"), angle = 0, angle.tol = 180, 
    nlag = 20, maxdist = 0.5, xlag = NA, lag.tol = 0.5, grid = TRUE) 
{
    if (!inherits(object, "asreml")) 
        stop("\nObject must be of class asreml\n")
    call <- object$call
    control <- object$control
    drop.unused.levels <- control$drop.unused.levels
    fixed <- as.formula(object$fixed.formula)
    if (missingRcov <- is.null(call$rcov)) 
        rcov <- ~NULL
    else rcov <- as.formula(call$rcov)
    if (is.null(call$random)) 
        random <- ~NULL
    else random <- as.formula(object$random.formula)
    if (is.null(call$sparse)) 
        sparse <- ~NULL
    else sparse <- as.formula(object$sparse.formula)
    varNames <- names(eval(call$data, parent.frame()))
    groups <- control$group
    if (length(groups) > 0) {
        for (i in names(groups)) {
            if (is.numeric(groups[[i]])) 
                groups[[i]] <- varNames[groups[[i]]]
            else {
                if (any(is.na(match(groups[[i]], varNames)))) 
                  stop(paste("Object in group", names(groups[1]), 
                    "not in data."))
            }
        }
    }
    mbfList <- list()
    if (length(control$mbf) > 0) {
        mbfList <- asreml.getMbf(control$mbf, eval(call$data))
    }
    weights <- as.character(call$weights)
    offset <- as.character(call$offset)
    ran.order <- call$ran.order
    if (is.null(ran.order)) 
        ran.order <- "user"
    if (is.null(call$as.multivariate)) 
        as.multivariate <- NULL
    else as.multivariate <- call$as.multivariate
    if (is.null(call$na.method.Y)) 
        na.method.Y <- "include"
    else na.method.Y <- call$na.method.Y
    if (is.null(call$na.method.X)) 
        na.method.X <- "include"
    else na.method.X <- call$na.method.X
    form <- asreml.ModelFormula(fixed, random, sparse, rcov, 
        weights, offset, groups, mbfList$mbfAttr, ignore = c("mv", 
            "trait", "units"), varNames)
    model.y <- attr(form, "model.y")
    multivariate <- (length(model.y) > 1) || (!is.null(as.multivariate))
    if (is.null(call$subset)) 
        mf <- as.call(list(as.name("model.frame"), formula = form, 
            data = call$data, drop.unused.levels = drop.unused.levels, 
            na.action = "na.pass"))
    else mf <- as.call(list(as.name("model.frame"), formula = form, 
        data = call$data, subset = call$subset, drop.unused.levels = drop.unused.levels, 
        na.action = "na.pass"))
    mf[[1]] <- as.name("model.frame")
    attr(mf$formula, ".Environment") <- sys.frame(sys.parent())
    data <- eval(mf, sys.frame(sys.parent()))
    if (missingRcov) 
        rcov <- asreml.formula(parse(text = "~ units"))
    family <- asreml.getFamily(call$family)
    asr.glm <- asreml.glm(family)
    asrAttrib <- list(Control = control, GROUP = groups, MBF = mbfList)
    for (a in names(asrAttrib)) attr(data, a) <- asrAttrib[[a]]
    data <- asreml.data(data, model.y, fixed, random, sparse, 
        rcov, asr.glm, asrAttrib, na.method.X, na.method.Y, as.multivariate, 
        weights, offset, ran.order)
    thrlevels <- attr(data, "thrlevels")
    keep <- attr(data, "keep")
    ntrt <- attr(data, "ntrt")
    model.y <- attr(data, "model.y")
    mvtrait <- attr(data, "mvtrait")
    section <- asreml.sectionName(rcov, data)
    if (is.null(section)) {
        nsect <- 1
    }
    else nsect <- length(unique(data[[section]]))
    if (!is.null(call$aom)) 
        aom <- as.logical(as.character(call$aom))
    else if (!is.null(object$control$aom)) 
        aom <- object$control$aom
    else aom <- asreml.control()$aom
    if (asr.glm$id == 1) {
        if (aom) {
            type <- "stdCond"
        }
        else {
            type <- "working"
        }
    }
    else {
        type <- "deviance"
    }
    if (length(formula[[2]]) == 0) {
        vg <- NULL
        grp <- NULL
        if (nsect == 1) {
            if (multivariate) {
                dimensions <- dimensions[dimensions != mvtrait]
                for (i in unique(as.character(data[, mvtrait]))) {
                  which <- (as.character(data[, mvtrait]) == 
                    i)
                  vg <- rbind(vg, asreml.variogram(x = data[which, 
                    dimensions], z = residuals(object, type = type)[which], 
                    composite = composite, model = model, metric = metric, 
                    angle = angle, angle.tol = angle.tol, nlag = nlag, 
                    maxdist = maxdist, xlag = xlag, lag.tol = lag.tol, 
                    grid = grid))
                  grp <- c(grp, as.character(data[which, mvtrait]))
                }
                vg[["groups"]] <- grp
            }
            else {
                dimensions <- names(object$R.param[[1]])
                dimensions <- dimensions[-length(dimensions)]
                vg <- asreml.variogram(x = data[, dimensions], 
                  z = residuals(object, type = type), composite = composite, 
                  model = model, metric = metric, angle = angle, 
                  angle.tol = angle.tol, nlag = nlag, maxdist = maxdist, 
                  xlag = xlag, lag.tol = lag.tol, grid = grid)
            }
        }
        else {
            if (multivariate) {
                dimensions <- dimensions[dimensions != mvtrait]
                for (i in unique(paste(as.character(data[, section]), 
                  as.character(data[, mvtrait])))) {
                  which <- (paste(as.character(data[, section]), 
                    as.character(data[, mvtrait])) == i)
                  vg <- rbind(vg, asreml.variogram(x = data[which, 
                    dimensions], z = residuals(object, type = type)[which], 
                    composite = composite, model = model, metric = metric, 
                    angle = angle, angle.tol = angle.tol, nlag = nlag, 
                    maxdist = maxdist, xlag = xlag, lag.tol = lag.tol, 
                    grid = grid))
                  grp <- c(grp, paste(as.character(data[, section]), 
                    as.character(data[, mvtrait])))
                }
                vg[["groups"]] <- grp
            }
            else {
                for (i in unique(as.character(data[, section]))) {
                  which <- (as.character(data[, section]) == 
                    i)
                  dimensions <- names(object$R.param[[i]])
                  dimensions <- dimensions[-length(dimensions)]
                  vg <- rbind(vg, asreml.variogram(x = data[which, 
                    dimensions], z = residuals(object, type = type)[which], 
                    composite = composite, model = model, metric = metric, 
                    angle = angle, angle.tol = angle.tol, nlag = nlag, 
                    maxdist = maxdist, xlag = xlag, lag.tol = lag.tol, 
                    grid = grid))
                  grp <- c(grp, as.character(data[which, section]))
                }
                vg[["groups"]] <- grp
            }
        }
    }
    else {
        rmatch <- function(x, y) {
            nx <- nchar(x)
            ny <- nchar(y)
            for (i in seq(ny - nx + 1, 1, -1)) {
                if (x == substring(y, i, i + nx - 1)) 
                  return(i)
            }
            return(NA)
        }
        if (length(formula) != 3) 
            stop("'formula' must be of the form 'response ~ predictors'")
        rsp <- dimnames(attr(terms(formula), "factors"))[[1]][1]
        rr <- coef(object, pattern = rsp)
        ff <- asreml.getOp(formula[[3]], "|")
        if (!is.null(ff)) {
            dimensions <- attr(terms(formula(paste("~", deparse(ff[[2]])))), 
                "term.labels")
            section <- attr(terms(formula(paste("~", deparse(ff[[3]])))), 
                "term.labels")
            if (length(section) > 1) 
                stop("Only one grouping variable permitted")
        }
        else {
            dimensions <- attr(terms(formula), "term.labels")
            section <- character(0)
            nsect <- 1
        }
        df <- list()
        for (i in c(section, dimensions)) {
            v <- sapply(strsplit(dimnames(rr)[[1]], ":", fixed = TRUE), 
                function(x, i) {
                  y <- x[grep(i, x)]
                  y <- substring(y, rmatch(i, y), nchar(y))
                  y
                }, i)
            df[[i]] <- substring(v, grep("_", substring(v[1], 
                1:nchar(v[1]), 1:nchar(v[1])))[1] + 1)
        }
        names(df) <- c(section, dimensions)
        df <- data.frame(df)
        if (length(section)) 
            nsect <- length(unique(df[[section]]))
        vg <- NULL
        grp <- NULL
        if (nsect == 1) {
            vg <- asreml.variogram(x = df[, dimensions], z = rr, 
                composite = composite, model = model, metric = metric, 
                angle = angle, angle.tol = angle.tol, nlag = nlag, 
                maxdist = maxdist, xlag = xlag, lag.tol = lag.tol, 
                grid = grid)
        }
        else {
            for (i in unique(as.character(df[, section]))) {
                which <- (as.character(df[, section]) == i)
                vg <- rbind(vg, asreml.variogram(x = df[which, 
                  dimensions], z = rr[which], composite = composite, 
                  model = model, metric = metric, angle = angle, 
                  angle.tol = angle.tol, nlag = nlag, maxdist = maxdist, 
                  xlag = xlag, lag.tol = lag.tol, grid = grid))
                grp <- c(grp, as.character(df[which, section]))
            }
            vg[["groups"]] <- grp
        }
    }
    vg <- data.frame(as.list(vg))
    grid <- is.na(match("angle", names(vg)))
    if (grid) {
        if (length(dimensions) == 1) 
            vv.names <- c(dimensions, "y", "gamma", "np")
        else vv.names <- c(dimensions[1], dimensions[2], "gamma", 
            "np")
        if (length(grp) > 0) 
            vv.names <- c(vv.names, "vv.groups")
        names(vg) <- vv.names
    }
    class(vg) <- c("asrVariogram", "data.frame")
    vg
}
wald.asreml <-
function (object = NULL, Ftest = formula("~NULL"), denDF = c("none", 
    "default", "numeric", "algebraic"), ssType = c("incremental", 
    "conditional"), maxiter = 10, ...) 
{
    if (is.null(object$aovTbl)) 
        stop("The model contains no fixed terms")
    mc <- match.call(expand.dots = TRUE)
    denDF <- match.arg(denDF)
    ssType <- match.arg(ssType)
    if (is.null(object)) 
        stop("Missing asreml object\n")
    if (length(Ftest[[2]]) > 0) {
        ssType <- "conditional"
        denDF <- "default"
    }
    if (ssType == "incremental" & denDF == "none") {
        s2 <- object$sigma2
        nsect <- length(object$R.param)
        xx <- seq(1, length(object$gammas))[object$gammas.type == 
            1]
        x <- xx[length(xx)]
        ntrt <- 1
        if (ntrt + nsect == 2) 
            convrt <- object$gammas[x] == 1 & object$gammas.con[x] != 
                4
        else convrt <- FALSE
        ss <- 1
        if (convrt) 
            ss <- s2
        nfact <- length(object$yssqu)
        aod <- matrix(nrow = nfact, ncol = 4)
        aod[, 1] <- object$noeff[seq(1, nfact)]
        aod[, 2] <- object$yssqu
        aod[, 3] <- aod[, 2]/ss
        aod[, 4] <- 1 - pchisq(aod[, 3], aod[, 1])
        aod <- data.frame(aod)
        aod <- rbind(aod, c(NA, object$sigma2, NA, NA))
        heading <- c("Wald tests for fixed effects\n", paste("Response: ", 
            as.character(object$fixed.formula[2]), "\n", sep = ""), 
            "Terms added sequentially; adjusted for those above\n")
        attr(aod, "heading") <- heading
        attr(aod, "names") <- c("Df", "Sum of Sq", "Wald statistic", 
            "Pr(Chisq)")
        attr(aod, "row.names") <- c(names(object$yssqu), "residual (MS)")
        aod <- aod[c(rev(seq(1, nrow(aod) - 1)), nrow(aod)), 
            ]
        oldClass(aod) <- c("anova", "data.frame")
        return(aod)
    }
    whichDF <- c(-1, 0, 1, 6)
    names(whichDF) <- c("none", "default", "numeric", "algebraic")
    whichSS <- c(1, 3)
    names(whichSS) <- c("incremental", "conditional")
    denDF <- whichDF[denDF]
    ssType <- whichSS[ssType]
    if (object$family$family[1] == "gaussian") 
        mi <- maxiter
    else mi <- object$control$maxiter
    object <- update(object, denDF = as.numeric(eval(denDF)), 
        ssType = as.numeric(eval(ssType)), Ftest = Ftest, maxiter = mi, 
        ...)
    if (length(mc$debug) > 0 && mc$debug) 
        return(object)
    aov <- object$aovTbl
    aov[, 3] <- round(aov[, 3], 1)
    aov[, 4] <- signif(aov[, 4], 4)
    aov[, 5] <- signif(aov[, 5], 4)
    if (length(Ftest[[2]]) > 0) {
        cvec <- rep("  ", nrow(aov))
        x <- aov[, 6]
        x[x == 0] <- -64
        temp <- .C("int2char", cvec = cvec, as.integer(x + 64), 
            as.integer(length(x)))
        which <- seq(1, length(x))[as.character(temp$cvec) == 
            "O"]
        aod <- matrix(nrow = 2, ncol = 3)
        aod[, 1] <- c(aov[which, 2], aov[which, 3])
        aod[, 2] <- c(aov[which, 5], NA)
        aod[, 3] <- c(1 - pf(aod[1, 2], aod[1, 1], aod[2, 1]), 
            NA)
        aod <- data.frame(aod)
        heading <- paste(as.character(deparse(Ftest)), "\n")
        attr(aod, "heading") <- heading
        attr(aod, "names") <- c("DF", "Conditional F", "Pr(F)")
        attr(aod, "row.names") <- c(dimnames(aov)[[1]][which], 
            "Residual")
        oldClass(aod) <- c("anova", "data.frame")
        return(aod)
    }
    else if ((ssType > 1) & (denDF >= 0)) {
        aod <- data.frame(Df = aov[, 2], denDF = aov[, 3], Finc = aov[, 
            4], Fcon = aov[, 5], Margin = aov[, 6], Pr = aov[, 
            7])
        cvec <- rep("  ", nrow(aod))
        aod[, 5][aod[, 5] == 0] <- -64
        temp <- .C("int2char", cvec = cvec, as.integer(aod[, 
            5] + 64), as.integer(nrow(aod)))
        aod[, 5] <- as.character(temp$cvec)
        heading <- c("Analysis of variance for fixed effects\n", 
            paste("Response: ", as.character(object$fixed.formula[2]), 
                "\n", sep = ""))
        attr(aod, "heading") <- heading
        attr(aod, "names") <- c("Df", "denDF", "F.inc", "F.con", 
            "Margin", "Pr")
        attr(aod, "row.names") <- dimnames(aov)[[1]]
    }
    else if (ssType > 1) {
        aod <- data.frame(Df = aov[, 2], Finc = aov[, 4], Fcon = aov[, 
            5], Margin = aov[, 6])
        cvec <- rep("  ", nrow(aod))
        aod[, 4][aod[, 4] == 0] <- -64
        temp <- .C("int2char", cvec = cvec, as.integer(aod[, 
            4] + 64), as.integer(nrow(aod)))
        aod[, 4] <- as.character(temp$cvec)
        heading <- c("Analysis of variance for fixed effects\n", 
            paste("Response: ", as.character(object$fixed.formula[2]), 
                "\n", sep = ""))
        attr(aod, "heading") <- heading
        attr(aod, "names") <- c("Df", "F.inc", "F.con", "Margin")
        attr(aod, "row.names") <- dimnames(aov)[[1]]
    }
    else if (denDF >= 0) {
        aod <- data.frame(Df = aov[, 2], denDF = aov[, 3], Finc = aov[, 
            4], Pr = aov[, 7])
        heading <- c("Analysis of variance for fixed effects\n", 
            paste("Response: ", as.character(object$fixed.formula[2]), 
                "\n", sep = ""))
        attr(aod, "heading") <- heading
        attr(aod, "names") <- c("Df", "denDF", "F.inc", "Pr")
        attr(aod, "row.names") <- dimnames(aov)[[1]]
    }
    oldClass(aod) <- c("data.frame")
    list(Wald = aod, stratumVariances = object$stratumVariances)
}
wald <-
function (x, ...) 
{
    UseMethod("wald")
}
multinomial <-
function (link = "logit") 
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    okLinks <- c("logit", "probit", "cloglog")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for multinomial; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(mu > 0) && all(mu < 1)
    aic <- function(y, n, mu, wt, dev) {
        m <- if (any(n > 1)) 
            n
        else wt
        -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
            y), round(m), mu, log = TRUE))
    }
    initialize <- expression({
        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1]
            n <- rep.int(1, nobs)
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            m <- weights * y
            if (any(abs(m - round(m)) > 0.001)) warning("non-integer successes in multinomial glm!")
        } else if (NCOL(y) == 2) {
            if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in multinomial glm!")
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)
        } else stop("for the binomial family, y must be a vector of 0 and 1's\n", 
            "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })
    structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, aic = aic, 
        mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu, 
        valideta = stats$valideta), class = "family")
}
negative.binomial <-
function (link = "log", phi = 1) 
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "inverse")) 
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"inverse\"")
    }
    env <- new.env(parent = .GlobalEnv)
    assign(".Phi", phi, envir = env)
    variance <- function(mu) mu + mu^2/.Phi
    validmu <- function(mu) all(mu > 0)
    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1, 
        y)/mu) - (y + .Phi) * log((y + .Phi)/(mu + .Phi)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Phi) * log(mu + .Phi) - y * log(mu) + lgamma(y + 
            1) - .Phi * log(.Phi) + lgamma(.Phi) - lgamma(.Phi + 
            y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- env
    famname <- "negative.binomial"
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta), class = "family")
}
variogram <-
function (x, ...) 
{
    UseMethod("variogram")
}
plot.asrVariogram <-
function (obj, npanels = NA, scale = TRUE, ...) 
{
    if (!inherits(obj, "asrVariogram")) 
        stop("\nObject must be of class asrVariogram\n")
    if (scale && !is.null(obj$vv.groups)) 
        obj$gamma <- unlist(tapply(obj$gamma, obj$vv.groups, 
            function(x) x/max(x))[unique(obj$vv.groups)], use.names = FALSE)
    grid <- is.na(match("angle", names(obj)))
    onepage <- TRUE
    if (grid) {
        tmp <- names(obj)
        xlab <- tmp[1]
        ylab <- tmp[2]
        tmp[1] <- "x"
        tmp[2] <- "y"
        names(obj) <- tmp
        if (!is.null(obj$vv.groups) & is.na(npanels)) 
            npanels <- length(unique(obj$vv.groups))
        if (!is.null(obj$vv.groups) && ceiling(length(unique(obj$vv.groups))/npanels) > 
            1) 
            onepage <- FALSE
        nd <- 2 - as.numeric(all(obj$y == 0))
        if (nd == 1) {
            if (is.null(obj$vv.groups)) 
                trellis.obj <- xyplot(gamma ~ x, xlab = paste(xlab, 
                  "(lag)"), data = obj, ...)
            else trellis.obj <- xyplot(gamma ~ x | vv.groups, 
                scales = list(relation = "free"), data = obj, 
                xlab = paste(xlab, "(lag)"), layout = c(0, npanels), 
                ...)
        }
        else {
            if (is.null(obj$vv.groups)) 
                trellis.obj <- wireframe(gamma ~ x * y, data = obj, 
                  xlab = list(label = paste(xlab, "\n(lag)"), 
                    cex = 0.66), ylab = list(label = paste(ylab, 
                    "\n(lag)"), cex = 0.66), screen = list(z = 30, 
                    x = -60, y = 0), aspect = c(1, 0.66), scales = list(distance = c(1.3, 
                    1.3, 1), arrows = FALSE, cex = 0.66), zlab = "", 
                  ...)
            else trellis.obj <- wireframe(gamma ~ x * y | vv.groups, 
                data = obj, xlab = list(label = paste(xlab, "\n(lag)"), 
                  cex = 0.66), ylab = list(label = paste(ylab, 
                  "\n(lag)"), cex = 0.66), screen = list(z = 30, 
                  x = -60, y = 0), aspect = c(1, 0.66), scales = list(distance = c(1.3, 
                  1.3, 1), arrows = FALSE, cex = 0.66, relation = "free"), 
                zlab = "", layout = c(0, npanels), ...)
        }
    }
    else {
        if (is.na(npanels)) 
            npanels <- length(unique(obj$angle))
        if (!is.null(obj$vv.groups) && ceiling(length(unique(obj$angle)) * 
            length(unique(obj$vv.groups))/npanels) > 1) 
            onepage <- FALSE
        if (is.null(obj$vv.groups)) 
            trellis.obj <- xyplot(gamma ~ distance | angle, scales = list(relation = "free"), 
                xlab = "Distance", data = obj, layout = c(0, 
                  npanels), ...)
        else trellis.obj <- xyplot(gamma ~ distance | angle * 
            groups, scales = list(relation = "free"), xlab = "Distance", 
            data = obj, layout = c(0, npanels), ...)
    }
    if (onepage) 
        print(trellis.obj)
    else {
        grDevices::devAskNewPage(TRUE)
        print(trellis.obj)
        grDevices::devAskNewPage(FALSE)
    }
    invisible(trellis.obj)
}
tr <-
function (x, ...) 
{
    UseMethod("tr")
}
lrt <-
function (x, ...) 
{
    UseMethod("lrt")
}
svc <-
function (x, ...) 
{
    UseMethod("svc")
}
terms.order <-
function (form, keep.order) 
{
    xx <- asreml.getInterExp(form[[2]])
    if (is.null(xx)) 
        ffe <- attr(terms(form, keep.order = keep.order), "term.label")
    else ffe <- sapply(xx, function(x) paste(deparse(x), collapse = ""))
    trms <- attr(terms(form, keep.order = keep.order), "term.labels")
    lapply(trms, function(x, ffe) {
        ss <- strsplit(x, ":", fixed = TRUE)[[1]]
        if (length(ss) == 1) 
            terms(formula(paste("~", x)))
        else if (!is.na(which <- asreml.matchInteraction(x, ffe))) 
            terms(formula(paste("~", ffe[which])))
        else terms(formula(paste("~", x)))
    }, ffe)
}
checkPedigree <-
function (pedigree, fgen = list(character(0), 0.01), gender = character(0), 
    mv = c("NA", "0", "*", " ")) 
{
    is.missing <- function(pcol, mv) {
        isna <- is.element(as.character(pcol), mv)
        isna <- isna | is.na(pcol)
        return(isna)
    }
    if (any(dup <- duplicated(as.character(pedigree[, 1])))) {
        stop(cat("Duplicated individuals", "\n", as.character(pedigree[, 
            1])[dup], "\n"))
    }
    which <- seq(1, nrow(pedigree))[as.character(pedigree[, 1]) == 
        as.character(pedigree[, 2]) | as.character(pedigree[, 
        1]) == as.character(pedigree[, 3])]
    if (length(which <- which[!is.na(which)])) {
        cat(paste(which, pedigree[, 1][which]), sep = "\n")
        stop("Individuals appear as their own parent")
    }
    fmode <- FALSE
    xlink <- FALSE
    col4 <- -1
    col4Name <- NULL
    col4Type <- character(0)
    fg <- as.double(fgen[[2]])
    if (length(fgen[[1]]) > 0) {
        if (length(col4 <- pedigree[, fgen[[1]]]) == 0) 
            stop("Argument 'fgen' specifies a non-existent dataframe column.")
        else {
            fmode <- TRUE
            col4Name <- fgen[[1]]
            col4Type <- "fgen"
            w <- is.na(pedigree[, col4Name])
            if (sum(as.numeric(w)) > 0) {
                if (sum(as.numeric(!(is.missing(pedigree[w, 2], 
                  mv) | is.missing(pedigree[w, 3], mv)))) > 0) 
                  stop("fgen is 'NA' for non-base individuals")
            }
        }
    }
    if (length(gender) > 0) {
        if (fmode) 
            stop("Cannot specify both 'fgen' and 'gender'")
        xlink <- TRUE
        fg <- NA
        if (length(pedigree[, gender]) == 0) 
            stop("Argument 'gender' specifies non-existent dataframe column.")
        else if (!is.factor(pedigree[, gender])) 
            pedigree[, gender] <- factor(pedigree[, gender])
        sx <- levels(pedigree[, gender])
        if (length(sx) != 2) 
            stop("'gender' must only have two levels.")
        col4 <- as.numeric(pedigree[, gender])
        col4Name <- gender
        col4Type <- "gender"
    }
    lped <- nrow(pedigree)
    charPed <- character(lped * 3)
    for (i in 1:3) charPed[seq((i - 1) * lped + 1, i * lped)] <- as.character(pedigree[, 
        i])
    charPed[is.na(charPed)] <- "NA"
    lvls <- unique(charPed)
    if (any(!is.na(which <- match(mv, lvls)))) 
        lvls <- lvls[-which[!is.na(which)]]
    charPed <- NULL
    nan <- length(lvls)
    if (xlink && (nan > lped)) 
        stop("Missing individuals in sex linked pedigree")
    memumdad <- names(pedigree)[1:3]
    pedigree <- c(match(pedigree[, 1], lvls, nomatch = 0), match(pedigree[, 
        2], lvls, nomatch = 0), match(pedigree[, 3], lvls, nomatch = 0), 
        rep(0, 3 * (nan - lped)))
    if (fmode) 
        col4 <- c(col4, rep(0, (nan - lped)))
    error <- paste(rep(" ", 256), collapse = "")
    storage.mode(nan) <- "integer"
    storage.mode(lped) <- "integer"
    storage.mode(fg) <- "double"
    storage.mode(col4) <- "double"
    out <- .C("sortpedigree", lped, nan, lvls, newped = as.integer(pedigree), 
        fgsx = col4, fg, error = error, NAOK = TRUE)
    if (nchar(out$error) > 0) 
        stop(out$error)
    which <- out$newped > 0
    out$newped[which] <- lvls[out$newped[which]]
    out$newped <- data.frame(out$newped[seq(1, nan)], out$newped[seq(nan + 
        1, 2 * nan)], out$newped[seq(2 * nan + 1, 3 * nan)], 
        stringsAsFactors = FALSE)
    if (length(col4Name)) 
        out$newped <- cbind(out$newped, out$fgsx)
    names(out$newped) <- c(memumdad, col4Name)
    attr(out$newped, "rowNames") <- lvls
    if (length(col4Type)) 
        attr(out$newped, "col4Type") <- col4Type
    return(out$newped)
}
database.path <-
function (x) 
{
    x
}
na.include <-
function (x) 
{
    if (inherits(x, "data.frame")) {
        for (i in seq(along = x)) x[[i]] <- na.include(x[[i]])
        x
    }
    else {
        if (!((inherits(x, "factor") || !is.null(attr(x, "levels"))) && 
            any(pos <- is.na(x)))) 
            return(x)
        a <- attributes(x)
        a$levels <- l <- c(a$levels, "NA")
        xx <- as.vector(unclass(x))
        xx[pos] <- length(l)
        attributes(xx) <- a
        xx
    }
}
