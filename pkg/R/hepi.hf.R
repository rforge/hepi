#Hedonic function estimators
build.hf.lm <- function(learndata, full.formula, min.formula, 
    backtrans = I, rm.infl=TRUE, description = NULL, 
    return.row.labels = FALSE, allow.variable.selection = TRUE) {
    
    if(!is(full.formula, "formula"))
        stop("'full.formula' must be a proper formula object")

    if(allow.variable.selection && !requireNamespace("MASS"))
        stop("the 'MASS' package is needed for variable selection")
        
    if(allow.variable.selection && !is(min.formula, "formula"))
        stop("'min.formula' must be a proper formula object")

    if(!is(backtrans, "function"))
        stop("'backtrans' must be a proper function")

    learndata.build.hf.lm <- learndata[complete.cases(learndata),]

    if(rm.infl){
         hf.model <- lm(full.formula, data=learndata.build.hf.lm)
         hf.dffits <- dffits(hf.model)
         data.noninfl <- !(abs(hf.dffits) > 2 * sqrt(summary(hf.model)$df[1]/(summary(hf.model)$df[1]+summary(hf.model)$df[2])))
         learndata.build.hf.lm <- learndata.build.hf.lm[data.noninfl,]
         rm(hf.dffits, data.noninfl, hf.model)
    }
    row.labels <- labels(learndata.build.hf.lm)[[1]]

    hf.model <- lm(full.formula, data = learndata.build.hf.lm)

    if(allow.variable.selection) {
        if("learndata.build.hf.lm" %in% ls(.GlobalEnv))
            stop("an object named 'learndata.build.hf.lm' must not exist in '.GlobalEnv' since it would be overwritten")
    
        assign("learndata.build.hf.lm", learndata.build.hf.lm, envir = .GlobalEnv)
        hf.model <- MASS::stepAIC(hf.model, scope = list(upper = full.formula, lower = min.formula))
        remove("learndata.build.hf.lm", envir = .GlobalEnv)
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

    hf <- function (data, interval = c("none", "parametric"), level = 0.95)
    {
        interval <- match.arg(interval)
        data <- purify(data, hf.model)
        res <- backtrans(predict(hf.model, newdata = data,
            na.action = NULL, interval=switch(interval, 
            none="none", parametric="prediction"), level=level))
        res
    }
    
    if(is.null(description)) {
        description <- "Hedonic function built upon a linear model."
        if(rm.infl) description <- paste(description, "Influential observations removed based on DFFITS.")
        if(allow.variable.selection) description <- paste(description, "Model selection using stepAIC.")
    }


    hf <- hedonic.function(hf, 
           env=list(hf.model = hf.model, purify = purify, backtrans = backtrans),
           characteristics.names = all.vars(formula(hf.model)),
           call = match.call(),
           description = description)

    if(return.row.labels) {
        list(hf = hf, row.labels = unique(row.labels))
    } else {
        hf
    }
}
