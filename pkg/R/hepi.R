.onLoad <- function(lib, pkg) require(methods)

.validHF <- function(object) {
    if(!is.element("data", names(formals(object@.Data))))
        return("The candidate hedonic function needs to have a formal argument named 'data'")
    else
        return(TRUE)
}

setClass("hedonic.function",
    representation("function",
        characteristics.names = "character",
        call = "call",
        description = "character"),
    validity = .validHF)


#Make closure (Lexical scoping, see R-FAQ, Sect. 3.3.1)
.hf.MC <- function(f, env) f

hedonic.function <- function(hf, characteristics.names, env = NULL, call = match.call(), description = "") {
    new("hedonic.function", .hf.MC(hf, env), characteristics.names = characteristics.names,
            call = call, description = description)
}

is.applicable.hf <- function(hf, data) {
    if(!is(hf, "hedonic.function"))
        stop("'hf' cannot be treated as from the class 'hedonic.function'")
    setequal(intersect(hf@characteristics.names, names(data)), hf@characteristics.names)
}

analyse.in.hf <- function(expr, hf) {
    if(!is(hf, "hedonic.function"))
        stop("'hf' cannot be treated as from the class 'hedonic.function'")
    eval.parent(substitute(eval(quote(expr), envir=environment(hf))))
}


#Hedonic elementary price indices
hepi <- function(hf0, hf1, M, type = c("jevons", "dutot", "carli", "hdutot", "hcarli"), na.rm = TRUE, debug = FALSE)
  {
    #Check whether arguments are valid
    if(!is(hf0, "hedonic.function"))
        stop("'hf0' cannot be treated as from the class 'hedonic.function'")

    if(!is(hf1, "hedonic.function"))
        stop("'hf1' cannot be treated as from the class 'hedonic.function'")
        
    if(!is.applicable.hf(hf0, M))
        stop("'hf0' cannot be applied to M")
        
    if(!is.applicable.hf(hf1, M))
        stop("'hf1' cannot be applied to M")
        
    type <- match.arg(type, several.ok = TRUE)
        
    #Predict prices for given reference characteristics
    p0hat <- hf0(M)
    p1hat <- hf1(M)

    if(length(p1hat) != length(p0hat))
        stop("'hf0(M)' and 'hf1(M)' have different lengths")

    #Calculate ratios of predicted prices
    ratios <- p1hat/p0hat

    #Return index
    calc.index <- function(type) {
      switch(type,
           carli = mean(ratios, na.rm=na.rm),
           dutot = mean(p1hat, na.rm=na.rm)/mean(p0hat, na.rm=na.rm),
           jevons = exp(mean(log(ratios), na.rm=na.rm)),
           hcarli = mean(ratios^(-1), na.rm=na.rm)^(-1),
           hdutot = (mean(p1hat^(-1), na.rm=na.rm))^(-1)/(mean(p0hat^(-1), na.rm=na.rm))^(-1),
           )
    }
    index <- sapply(type, calc.index)
    names(index) <- type
    
    if (!debug)
        index
    else
        list(index = index, p0hat = p0hat, p1hat = p1hat, ratios = ratios)
  }

#Wrapper functions for individual index types
hepi.jevons <- function(hf0, hf1, M) as.numeric(hepi(type="jevons", hf0=hf0, hf1=hf1, M=M))
hepi.carli <- function(hf0, hf1, M) as.numeric(hepi(type="carli", hf0=hf0, hf1=hf1, M=M))
hepi.dutot <- function(hf0, hf1, M) as.numeric(hepi(type="dutot", hf0=hf0, hf1=hf1, M=M))
hepi.hcarli <- function(hf0, hf1, M) as.numeric(hepi(type="hcarli", hf0=hf0, hf1=hf1, M=M))
hepi.hdutot <- function(hf0, hf1, M) as.numeric(hepi(type="hdutot", hf0=hf0, hf1=hf1, M=M))
