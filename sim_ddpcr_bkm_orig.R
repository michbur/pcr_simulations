sim_ddpcr_bkm_orig <- function(m, n=20000, mexp=T, n_exp = 8, pos_sums=F, fluo = NULL,sddropc=0,mudropr=1,sddropr=0,Pvar=T,piperr=0,seed=runif(1),dropsd=0,falpos=0,falneg=0,rain=0){
  
  ##############
  ### checks ###
  ##############
  
  if(!is.numeric(seed)) stop("seed must have a numeric argument", call. = TRUE, domain = NA)
  set.seed(seed)
  
  testInteger <- function(x){
    test <- all.equal(x, as.integer(x), check.attributes = FALSE)
    isTRUE(test)
  }
  
  if(!is.logical(mexp)) stop("mexp must be a logical argument (TRUE or FALSE)", call. = TRUE, domain = NA)
  if(max(!is.finite(m))) stop("Concentrations should all be numeric", call. = TRUE, domain = NA)
  if(min(m)<0) stop("Concentrations cannot be negative", call. = TRUE, domain = NA)
  if(mexp){lambda <- m}else{lambda <- m*0.89/1000}
  
  if(!is.numeric(n)) stop("number of droplets must have a numeric argument", call. = TRUE, domain = NA)
  if(n < 10) stop("number of droplets must be large", call. = TRUE, domain = NA)
  if(!testInteger(n)) {warning("number of droplets will be rounded up to the next integer", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                       n <- ceiling(n)}
  
  if(!is.numeric(n_exp)) stop("number of replicates must have a numeric argument", call. = TRUE, domain = NA)
  if(n_exp < 1) stop("number of replicates must be at least 1", call. = TRUE, domain = NA)
  if(!testInteger(n_exp)) {warning("number of replicates will be rounded up to the next integer", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                           n_exp <- ceiling(n_exp)}
  
  if(!is.logical(pos_sums)) stop("pos_sums must be a logical argument (TRUE or FALSE)", call. = TRUE, domain = NA)
  
  if(!is.numeric(sddropc)) stop("sddropc must have a numeric argument", call. = TRUE, domain = NA)
  if(sddropc < 0) {warning("sddropc will be set to 0", call. = TRUE, domain = NA)
                   sddropc <- 0}
  if(sddropc > n/5) {warning("sddropc will be set to n/5", call. = TRUE, domain = NA)
                     sddropc <- n/5}
  
  if(!is.numeric(mudropr)) stop("mudropr must have a numeric argument", call. = TRUE, domain = NA)
  if(mudropr*n < 10) stop("mudropr too small, too few droplets will be returned", call. = TRUE, domain = NA)
  if(mudropr > 1) warning("mudropr will be set to 1", call. = TRUE, domain = NA) # happens in code
  
  if(mudropr < 1){if(!is.numeric(sddropr)) stop("sddropr must have a numeric argument", call. = TRUE, domain = NA)
                  if(sddropr < 0) {warning("sddropr will be set to 0", call. = TRUE, domain = NA)
                                   sddropr <- 0}
                  if((sddropr >= mudropr) | ((sddropr+mudropr) >= 1)) stop("sddropr too large", call. = TRUE, domain = NA)
  }
  
  if(!is.numeric(piperr)) stop("pipette error must have a numeric argument", call. = TRUE, domain = NA)
  if(piperr < 0) stop("pipette error should be positive or 0", call. = TRUE, domain = NA)
  
  if(!is.numeric(dropsd)) stop("dropsd must have a numeric argument", call. = TRUE, domain = NA)
  if(dropsd < 0) {warning("dropsd will be set to 0", call. = TRUE, domain = NA)
                  dropsd <- 0}
  
  if(!is.logical(Pvar)) stop("Pvar must be a logical argument (TRUE or FALSE)", call. = TRUE, domain = NA)
  
  if(is.null(fluo)){fluoselect <- 1
                    if(!is.numeric(falpos)) stop("falpos must have a numeric argument", call. = TRUE, domain = NA)
                    if(falpos < 0) {warning("falpos will be set to 0", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                                    falpos <- 0}
                    if(!is.numeric(falneg)) stop("falneg must have a numeric argument", call. = TRUE, domain = NA)
                    if(falneg < 0) {warning("falneg will be set to 0", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                                    falneg <- 0}
                    if(falpos >= 1) stop("falpos too large, set a number between 0 and 1", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                    if(falneg >= 1) stop("falneg too large, set a number between 0 and 1", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
  }
  else{
    if(!is.numeric(rain)) stop("rain must have a numeric argument", call. = TRUE, domain = NA)
    if(rain < 0) {warning("rain will be set to 0", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                  rain <- 0}
    if(rain >= 1) {warning("rain will be set to 1", call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NA)
                   rain <- 1}
    if(is.logical(fluo)){fluoselect <- ifelse(fluo,2,1)}
    else{if(is.numeric(fluo)){if(fluo>0){fluoselect <- 3}else{stop("fluo does not have a valid argument", call. = TRUE, domain = NA)}}
         else{stop("fluo does not have a valid argument", call. = TRUE, domain = NA)}
    }
  }
  
  
  ###############
  ### repfunc ###
  ###############
  
  # Same procedure for all replicates
  # repfunc is called internally in samfunc
  repfunc <- function(repdat){
    dropmem <- sample(repdat[1],repdat[2],replace=T,prob=rlnorm(repdat[1],0,dropsd))
    # droplet membership, probability proportional to size, size following a lognormal distribution
    dropn <- ifelse(mudropr>=1,repdat[1],round(repdat[1]*plogis(rnorm(1,log(mudropr/(1-mudropr)),log((mudropr+sddropr)/(mudropr-sddropr)*(1-mudropr+sddropr)/(1-mudropr-sddropr))/2))))
    # number of droplets retained
    dropmem <- dropmem[dropmem<=dropn]
    # only retain copies of which the droplet is retained (lower rank)
    if(fluoselect==1){
      if(!pos_sums){dropno <- dropn - length(as.vector(table(dropmem)))
                    # number of droplets without copy (total - number with copies)
                    droppos <- dropn-(rbinom(1,dropno,1-falpos)+rbinom(1,dropn-dropno,falneg))
                    return.drops <- c(droppos,dropn)}
      # number of droplets with a negative signal (true neg + false neg)
      else{dropvec <- sapply(1:dropn,function(i){any(dropmem==i)})
           # vector with TRUE for positive droplets and FALSE for negative
           dropfin <- (sapply(dropvec,function(x){ifelse(x,rbinom(1,1,falneg),rbinom(1,1,falpos))}))%%2
           return.drops <- dropfin}
      # vector TRUE for positive signal and FALSE for negative signal
      return.fluo <- NULL
    }
    else{dropvec <- sapply(1:dropn,function(i){any(dropmem==i)})
         # vector with TRUE for positive droplets and FALSE for negative
         fluopeaks <- rnorm(dropn,1000,100)+8000*dropvec*(1-runif(dropn)^(1/rain-1))*(1-rain^2)+2000*(1-dropvec)*(1-runif(dropn)^rain)*(1-rain^2)
         # random variation+downward rain+upward rain
         dropfin <- (fluopeaks>2500)
         # hard threshold as in most software these days
         # vector TRUE for positive signal and FALSE for negative signal
         if(pos_sums){return.drops <- dropfin}
         else{return.drops <- c(sum(dropfin),dropn)}
         if(fluoselect==2){return.fluo <- fluopeaks}
         else{fluopos <- (1:dropn)*fluo+9+runif(dropn)*2+rnorm(dropn)
              # vector of positions where the peak was found
              # peaks on average 10 positions away from each other.
              fluox <- 1:(round(dropn*fluo+90,-2))
              fluoy <- rnorm(length(fluox),50,10)
              # define fluo vectors with random background
              j <- 1
              for(i in fluox){
                if(j<=dropn){
                  if(fluopos[j]<(i-20)){j <- j+1}
                  for(k in j:round(j+30/fluo)){
                    # move j such that only the influence of peaks close by are counted
                    # influence of peaks further away would be marginally small anyway
                    if(k <= dropn){
                      dist <- fluopos[k]-i
                      # distance between peak and current location
                      fluoy[i] <- fluoy[i]+dnorm(dist/2+rnorm(1,0,0.1))/dnorm(0)*fluopeaks[k]*0.95
                      # add fluorescence signal stemming from this specific droplet
                    }
                  }
                }
              }
              return.fluo <- fluoy
         }
    }
    return(list(drop=return.drops,fluo=return.fluo))
    # returns list with first element vector of droplets 1/0
    # or pair of number of pos droplets and total number of droplets
    # and second element either NULL, peak fluorescence of droplets
    # or vector with continuous fluorescence output 
  }
  
  
  ###############
  ### samfunc ###
  ###############
  
  # Same procedure for all simulations
  # samfunc is called in sim_ddpcr
  samfunc <- function(lambdan){
    dropstart <- round(rnorm(n_exp,n,sddropc))
    # number of droplets
    copyvar <- lambdan*rnorm(n_exp,1,piperr)*dropstart
    copyvar[copyvar<0] <- 0
    # expected number of copies after pipette variation
    copyn <- ifelse(rep(Pvar,n_exp),rpois(n_exp,copyvar),round(copyvar))
    # number of copies
    lamdummy <- rep(lambdan,n_exp)
    repdat <- as.list(data.frame(rbind(dropstart,copyn,lamdummy)))
    # number of droplets and copies in a list with n_exp elements, all pairs
    # droplets is the first element, copies the second
    repres <- sapply(repdat,repfunc)
    # returns a list with n_exp*2 elements with
    # first element vector of droplets 1/0 or
    #   pair of number of pos droplets and total number of droplets
    # second element either NULL, peak fluorescence of droplets or
    #   vector with continuous fluorescence output 
    # and so on for each replicate
  }
  
  
  ###############
  ### execute ###
  ###############
  
  out <- sapply(lambda,samfunc)
  # out returns a list with entries for each lambda
  # each entry is a list with n_exp*2 elements from samfunc
  # first element vector of droplets 1/0 or
  #   pair of number of pos droplets and total number of droplets
  # second element either NULL, peak fluorescence of droplets or
  #   vector with continuous fluorescence output
  # and so on
  
  # To do: create ddpcr object?
  # create_ddpcr(res, rep(n, n_exp), threshold = 2500, type = type)
}