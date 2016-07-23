build.hf.lm.split <- function(learndata, split.var, full.formula, min.formula,
    backtrans = I, rm.infl = FALSE, description = NULL, return.row.labels = FALSE,
    allow.variable.selection = TRUE, use.overall.hf = TRUE, 
    split.threshold = 100) {

    if(!is(full.formula, "formula"))
        stop("'full.formula' must be a proper formula object")

    if(allow.variable.selection && !requireNamespace("MASS"))
        stop("the 'MASS' package is needed for variable selection")
        
    if(allow.variable.selection && !is(min.formula, "formula"))
        stop("'min.formula' must be a proper formula object")

    if(!is(backtrans, "function"))
        stop("'backtrans' must be a proper function")
        
    if(!is(learndata[,split.var], "factor"))
        stop("'split.var' must be the name of a factor variable contained in `learndata'")

    learndata.build.hf.lm <- learndata

    if(use.overall.hf) {
        overall.hfc <- build.hf.lm(learndata, 
                        full.formula = full.formula, 
                        min.formula = min.formula,
                        backtrans = backtrans,
                        rm.infl = rm.infl, 
                        return.row.labels = TRUE, 
                        allow.variable.selection = allow.variable.selection)
        model.list <- list("overall" = analyse.in.hf(hf.model, overall.hfc[["hf"]]))
        row.labels <- overall.hfc[["row.labels"]]
    } else {
        model.list <- list()
        row.labels <- NULL
    }

    learndata.build.hf.lm[,split.var] <- factor(learndata.build.hf.lm[,split.var])

    for(cm in levels(learndata.build.hf.lm[,split.var]))
    {
        full.formula.m <- full.formula
    
        subset <- !is.na(learndata.build.hf.lm[,split.var]) & as.character(learndata.build.hf.lm[,split.var])==as.character(cm)
        if(any(subset)) {
            datam <- learndata[subset,all.vars(full.formula.m)]
            datam <- datam[complete.cases(datam),]

            if(dim(datam)[[1]]>split.threshold)
            {
                is.factor.variable <- sapply(dimnames(datam)[[2]], function(fv) is.factor(datam[,fv]))                
                for(factor.variable in dimnames(datam)[[2]][is.factor.variable]) {
                     if(length(levels(factor(datam[,factor.variable]))) <= 1) {
                         cat(paste("... Removing", factor.variable, "\n"))
                         full.formula.m <- update(full.formula.m, as.formula(paste(". ~ . - ",factor.variable)))
                     }
                 }
                 if(rm.infl){
                      hf.model <- lm(full.formula.m, data=datam)
                      hf.dffits <- dffits(hf.model)
                      data.noninfl <- !(abs(hf.dffits) > 2 * sqrt(summary(hf.model)$df[1]/(summary(hf.model)$df[1]+summary(hf.model)$df[2])))
                      datam <- datam[data.noninfl,]
                      rm(hf.dffits, data.noninfl, hf.model)
                      is.factor.variable <- sapply(dimnames(datam)[[2]], function(fv) is.factor(datam[,fv]))                
                      for(factor.variable in dimnames(datam)[[2]][is.factor.variable]) {
                           if(length(levels(factor(datam[,factor.variable]))) <= 1) {
                               cat(paste("... Removing", factor.variable, "\n"))
                               full.formula.m <- update(full.formula.m, as.formula(paste(". ~ . - ",factor.variable)))
                           }
                       }
                 }
                 hf.model <- try(lm(full.formula.m, data=datam), silent=TRUE)
                 cat(paste(format(cm, width=30), ":", summary(hf.model)$r.squared, "\n"))

                 if(allow.variable.selection && (class(hf.model) != "try-error")) {
                     if("datam" %in% ls(.GlobalEnv))
                         stop("an object named 'datam' must not exist in '.GlobalEnv' since it would be overwritten")

                     assign("datam", datam, envir = .GlobalEnv)
                     hf.model.n <- try(MASS::stepAIC(hf.model, scope = list(upper = full.formula.m, lower = min.formula), trace=0), silent=TRUE)
                     remove("datam", envir = .GlobalEnv)
                     if(class(hf.model.n) != "try-error") {
                         hf.model <- hf.model.n
                         cat(paste(format(" (stepAIC)", width=30), ":", summary(hf.model)$r.squared, "\n"))
                     } else {
                         cat(paste(format(" (stepAIC)", width=30), ": failed", "\n"))
                     }
                     remove(hf.model.n)
                 }

                 if(class(hf.model) != "try-error") {
                     hf.model.list <- list(hf.model)
                     names(hf.model.list) <- cm
                     model.list <- c(model.list, hf.model.list)
                     row.labels <- unique(c(row.labels, labels(datam)[[1]]))
                 }
            }
        }
    }

    purify <- function(data, model) {

        if (!is.null(data)) {
            data <- as.data.frame(data)
            unwantedentries <- rep(FALSE, nrow(data))

            for (var in names(model$xlevels)) {
                addlevels <- setdiff(levels(data[,var]), unlist(model$xlevels[var]))
                unwantedentries <- unwantedentries | is.element(data[,var], addlevels)
                data[unwantedentries, var] <- NA
                data[,var] <- factor(data[,var])
            }
        }
        data
    }

    predict.fun <- function(entries, split.var.value, interval = c("none", "parametric"), level = 0.95) {
        interval <- match.arg(interval)
        if(is.element(split.var.value,names(model.list))) {
            entries.p <- purify(entries, model.list[[split.var.value]])
            res <- backtrans(predict(model.list[[split.var.value]], newdata = entries.p,
                na.action = NULL, interval=switch(interval, none="none", parametric="prediction"), level=level))
            names(res) <- dimnames(entries.p)[[1]]
        } else {
            res <- rep(NA, dim(entries)[1])
            names(res) <- dimnames(entries.p)[[1]]
        }
        res
    }

    hf <- function (data, interval = c("none", "parametric"), level = 0.95)
    {
        interval <- match.arg(interval)
        res <- rep(NA, dim(data)[[1]])

        for(cm in names(model.list))
        {
            subset <- !is.na(data[,split.var]) & as.character(data[,split.var])==as.character(cm)
            if(any(subset)) {
                res[subset] <- predict.fun(entries = data[subset,], split.var.value = cm, interval = interval, level = level)
            }
        }

        if(!is.null(model.list[["overall"]])) {
            subset <- is.na(res)
            if(any(subset)) {
                res[subset] <- predict.fun(entries = data[subset,], split.var.value = "overall", interval = interval, level = level)
            }
        }

        res
    }

    if(is.null(description)) {
        description <- paste("Hedonic function built upon a linear model differentiating along '",split.var,"'.", sep="")
        if(rm.infl) description <- paste(description, "Influential observations removed based on DFFITS.")
        if(allow.variable.selection) description <- paste(description, "Model selection using stepAIC.")
    }

    hf <- 
        hedonic.function(hf, 
            env = list(model.list = model.list, purify = purify, predict.fun=predict.fun, backtrans = backtrans),
            characteristics.names = all.vars(full.formula),
            call = match.call(),
            description = description)

    if(return.row.labels) {
        list(hf = hf, row.labels = unique(row.labels))
    } else {
        hf
    }
}
