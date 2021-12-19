#' Kling-Gupta Efficiency
#' 
#' @param obs numeric 'data.frame', 'matrix' or 'vector' with observed values
#' @param sim numeric 'data.frame', 'matrix' or 'vector' with simulated values
#' @param s scaling factors
#' 
#' @return 
#' Kling-Gupta Efficiency between 'sim' and 'obs'
#' 
#' @references 
#' 1. Hoshin V. Gupta, Harald Kling, Koray K. Yilmaz, Guillermo F. Martinez,
#' Decomposition of the mean squared error and NSE performance criteria:
#' Implications for improving hydrological modelling,
#' Journal of Hydrology, Volume 377, Issues 1-2, 20 October 2009, Pages 80-91,
#' DOI: 10.1016/j.jhydrol.2009.08.003. ISSN 0022-1694
#' 
#' 2. Kling, H., M. Fuchs, and M. Paulin (2012), Runoff conditions in the upper
#' Danube basin under an ensemble of climate change scenarios,
#' Journal of Hydrology, Volumes 424-425, 6 March 2012, Pages 264-277,
#' DOI:10.1016/j.jhydrol.2012.01.011.
#' 
#' @export 
KGE.default <- function(obs, sim, s=c(1,1,1), na.rm=TRUE,
                        method=c("2009", "2012"), out.type=c("single", "full"), ...) {
  # If the user provided a value for 's'
  if (!identical(s, c(1,1,1)) )  {
    if ( length(s) != 3 ) stop("Invalid argument: lenght(s) must be equal to 3 !")
    if ( sum(s) != 1 )    stop("Invalid argument: sum(s) must be equal to 1.0 !")
  } # IF end

  method   <- match.arg(method)
  out.type <- match.arg(out.type)

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

  vi <- valindex(sim, obs)

  if (length(vi) > 0) {

    obs <- as.numeric(obs[vi])
    sim <- as.numeric(sim[vi])

    # Mean values
    mean.sim <- mean(sim, na.rm=na.rm)
    mean.obs <- mean(obs, na.rm=na.rm)

    # Standard deviations
    sigma.sim <- sd(sim, na.rm=na.rm)
    sigma.obs <- sd(obs, na.rm=na.rm)

    # Pearson product-moment correlation coefficient
    r     <- rPearson(sim, obs)

    # Alpha is a measure of relative variability between simulated and observed values (See Ref1)
    Alpha <- sigma.sim / sigma.obs

    # Beta is the ratio between the mean of the simulated values to the mean of observations
    Beta <- mean.sim / mean.obs

    # CV.sim is the coefficient of variation of the simulated values [dimensionless]
    # CV.obs is the coefficient of variation of the observations [dimensionless]
    CV.sim <- sigma.sim / mean.sim
    CV.obs <- sigma.obs / mean.obs

    # Gamma is the variability ratio, which is used instead of Alpha (See Ref2)
    Gamma <- CV.sim / CV.obs

    # Variability ratio depending on 'method'
    if(method=="2012") {
      vr     <- Gamma
      vr.stg <- "Gamma"
    } else {
      vr     <- Alpha
      vr.stg <- "Alpha"
    } # ELSE end

    # KGE Computation
    if ( (mean.obs != 0) | (sigma.obs != 0) ) {
      KGE <- 1 - sqrt( (s[1]*(r-1))^2 + (s[2]*(vr-1))^2 + (s[3]*(Beta-1))^2 )
    } else {
      if ( mean.obs != 0)  warning("Warning: 'mean(obs)==0'. Beta = Inf")
      if ( sigma.obs != 0) warning("Warning: 'sd(obs)==0'. ", vr.stg, " = Inf")
      KGE <- NA
    } # ELSE end
  } else {
    r    <- NA
    Beta <- NA
    vr   <- NA
    if(method=="2012") {
      vr.stg <- "Gamma"
    } else vr.stg <- "Alpha"
    KGE <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end

  if (out.type=="single") {
    out <- KGE
  } else {
    out <- list(KGE.value=KGE, KGE.elements=c(r, Beta, vr))
    names(out[[2]]) <- c("r", "Beta", vr.stg)
  } # ELSE end

  return(out)
} # 'KGE.default' end

rPearson <- function (sim, obs, ...) {
    if (is.na(match(class(sim), c("integer", "numeric", 
        "ts", "zoo"))) | is.na(match(class(obs), 
        c("integer", "numeric", "ts", "zoo")))) 
        stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
    rPearson <- cor(sim, obs, method = "pearson", use = "pairwise.complete.obs")
    return(rPearson)
}
