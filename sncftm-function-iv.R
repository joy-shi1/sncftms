# ------ G-ESTIMATION OF STRUCTUARL NESTED CUMULATIVE FAILURE TIME MODELS ------
# -------------------- USING INSTRUMENTAL VARIABLE ANALYSIS --------------------
#
# -------------------------------- BY: Joy Shi ---------------------------------
# -------------------------- LAST MODIFIED: 2021-12-06 -------------------------
#
# NOTES: 
# Please see simulated data as an example of how the data needs to be
# set up for the analysis
# 
# REQUIREMENTS:
# The function relies on packages 'optimx' and 'parallel' (for paralleling 
# the code)
#
# BASED ON THE SNCFTM SAS MACRO BY SALLY PICCIOTTO:
# For more information, refer to
# https://causalab.sph.harvard.edu/software/ and 
# https://pubmed.ncbi.nlm.nih.gov/24347749/

# ARGUMENTS:
#  - data: name of the data frame containing the variables in the model
#  - id: name of the variable (as a string) corresponding to participant index
#  - time: name of the variable (as a string) corresponding to time index 
#          (minimum must be 1)
#  - z: name of the variable (as a string) corresponding to the instrument
#  - z.modelvars: formula for the model for the instrument
#  - z.family: family for the model for the instrument (options are "gaussian" 
#    for linear regression; "binomial" for logistic regression)
#  - x: name of the variable (as a string) corresponding to the treatment
#  - y: name of the variable (as a string) corresponding to the outcome
#  - clost: name of the variable (as a string) corresponding to censoring due
#    to lost to follow-up
#  - clost.modelvars: formula for the model for censoring due to lost to 
#    follow-up
#  - cdeath: name of the variable (as a string) corresponding to censoring due
#    to death
#  - z.indicator: (optional) name of the variable (as a string) corresponding
#    to the indicator for which observations should be used in the model for 
#    for the treatment (e.g. only among controls in a nested case-control sample)
#  - z.timefixed: set to T if the instrument is time fixed (i.e. the model for
#    the instrument will be estimated only among observations at baseline;
#    otherwise will use all observation time points)
#  - death.modelvars: formula for the model for censoring due to death
#  - blipfunction: options are 1 (for 1+[exp(psi*Am)-1]/(k-m)) or 2 
#    (for psi*Am)
#  - start.value: starting value for grid search
#  - grid: set to T to obtain output from the estimating equation across
#    a range of psi value
#  - grid.range: range of psi values to calculate the estimating equation; a
#    single value c is given and the range is (+c, -c)
#  - grid.increment: increments of psi used to calculate the estimating equation
#  - blip.data: name of the data frame for blipping up and down; if NULL then
#    will use the data frame specified under the "data" argument
#  - blipupdown: set to T to obtain marginal cumulative risks under the "never
#    treat" and intervention regimes by blipping down and blipping up
#  - intervention.regime: specification of the treatment regime that is 
#    targeted when blipping up; should be a vector of the same length as 
#    the number of time points; if not specified, then the intervention.regime
#    will be set to "always treat"
#  - boot: set to T to obtain 95% CI by bootstrapping
#  - R: number of bootstraps
#  - parallel: set to T to parallelize
#  - seed: seed used for bootstrapping

# Installing and loading required packages
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('optimx')) install.packages('optimx'); library('optimx')

# SNCFTM function 
sncftm.iv <- function(data, id, time, z, z.modelvars=~1, z.family="gaussian", x, y,
                      clost=NULL, clost.modelvars=NULL, cdeath=NULL, cdeath.modelvars=NULL,
                      z.indicator=NULL, z.timefixed=T,
                      blipfunction, start.value=0,
                      grid=F, grid.range=1.5, grid.increment=0.01,
                      blip.data=NULL, blipupdown=T, intervention.regime=NULL,
                      boot=T, R=1000, parallel=T, seed=549274){
  
  # Data Prep and calculations not required for psi 
  if (parallel==T){numCores <- max(detectCores()-1, 1)}
  
  estf.dataprep <- function(data){
    d <- data.frame(data)
    
    # Predicted values from Z and censoring models
    if (z.timefixed==T){
      if (is.null(z.indicator)==T){
        if (z.family=="gaussian"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         data=d[which(d[[time]]==min(d[[time]])),])
        } else if (z.family=="binomial"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         family=binomial, data=d[which(d[[time]]==min(d[[time]])),])
        }
      } else if (is.null(z.indicator)==F){
        if (z.family=="gaussian"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         data=d[which(d[[time]]==min(d[[time]]) & d[[z.indicator]]==1),])
        } else if (z.family=="binomial"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         family=binomial, data=d[which(d[[time]]==min(d[[time]]) & d[[z.indicator]]==1),])
        }
      }
    } else if (z.timefixed==F){
      if (is.null(z.indicator)==T){
        if (z.family=="gaussian"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         data=d)
        } else if (z.family=="binomial"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         family=binomial, data=d)
        }
      } else if (is.null(z.indicator)==F){
        if (z.family=="gaussian"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         data=d[which(d[[z.indicator]]==1),])
        } else if (z.family=="binomial"){
          z.model <- glm(as.formula(paste(z, "~", paste(z.modelvars)[2], sep="")), 
                         family=binomial, data=d[which(d[[z.indicator]]==1),])
        }
      }      
    }
    d$z.pred <- predict(z.model, newdata=d, type="response")
    if (is.null(clost)==F & is.null(clost.modelvars)==F){
      clost.model <- glm(as.formula(paste(clost, "==0~", paste(clost.modelvars)[2], sep="")), family=binomial(), data=d)
      d$clost.pred <- predict(clost.model, d, type="response")
      d$clost <- d[,clost]
    } else{
      d$clost.pred <- 1
      d$clost <- 0
    }
    if (is.null(cdeath)==F & is.null(cdeath.modelvars)==F){
      cdeath.model <- glm(as.formula(paste(cdeath, "==0~", paste(cdeath.modelvars)[2], sep="")), family=binomial(), data=d)
      d$cdeath.pred <- predict(cdeath.model, d, type="response")
      d$cdeath <- d[,cdeath]
    } else{
      d$cdeath.pred <- 1
      d$cdeath <- 0
    }
    d <- d[order(d[,id], d[,time]),]
    n.followup <- length(unique(d[[time]]))
    
    # Creating dataset restricted to participants who ever had an event
    ever.y.id <- d[which(d[,y]==1),][[id]]
    ever.treat.id <- d[which(d[,x]!=0),][[id]]
    cero <- d[which(d[[id]] %in% ever.y.id),] # Restrict to ever had an event
    cero$ever_treat <- ifelse((cero[[id]] %in% ever.treat.id)==T, 1, 0) # Indicate if ever treated
    cero$count <- ave(rep(1, nrow(cero)), cero[[id]], FUN = sum) # Count for each ID
    cero <- cero[which(!is.na(cero[,z])),]  # Restrict to non-missing z  
    
    # Calculating contributions to estimating equation among untreated
    cero.untreated <- cero[which(cero$ever_treat==0),]
    if (nrow(cero.untreated)==0){
      newcov.untreated <- 0
      newu1.untreated <- 0
    } else{
      # Calculating censoring weights
      ## Note:
      ## Censoring weight for time=1 is cumulative product from time=1 to time=k for ID i
      ## Censoring weight for time=2 is cumulative product from time=2 to time=k for ID i
      ## ...etc.
      tpw.untreated <- 1/(cero.untreated[["clost.pred"]]*cero.untreated[["cdeath.pred"]]) # Weight for each time point
      cumw.untreated <- unname(ave(tpw.untreated, cero.untreated[[id]], FUN=prod)) # Cumulative product of weights per ID
      lagtpw.untreated <- suppressWarnings(unname(ave(tpw.untreated, cero.untreated[[id]], FUN=function(j) c(1, j[1:(length(j)-1)]), 1))) # Lagged weights by ID
      cumlagw.untreated <- unname(ave(lagtpw.untreated, cero.untreated[[id]], FUN=cumprod)) # Cumulative product of lagged weights
      w.untreated <- cumw.untreated/cumlagw.untreated
      
      # Contribution to estimating equation from H(psi)
      ## Note:
      ## Only contribution is when Y = 1
      ## Note that if, say, Y = 1 at k = 3, then also contribute to estimating equation at k = 4, 5, K (i.e. until end of follow-up)
      ## Contribution is exp(psi*A)=1 when A=0 (untreated)
      totalfu.untreated <- unname(ave(cero.untreated[[id]], cero.untreated[[id]], FUN=length)) # Duration of follow-up per ID
      tpcontributed.untreated <- n.followup+1-totalfu.untreated # Number of time points contributed to est eq
      hm.untreated <- tpcontributed.untreated*w.untreated # Contribution, weighted by censoring weights
      
      # Multiply H(psi) by Z-E[Z]
      u0i.untreated <- hm.untreated * (cero.untreated[[z]]-cero.untreated[["z.pred"]])
      
      # Calculating covariace matrix
      newu1.untreated <- aggregate(u0i.untreated, list(cero.untreated[[id]]), FUN=sum)[,2]
      newcov.untreated <- newu1.untreated %*% newu1.untreated
    }
  
    # Calculating contributions to estimating equation among treated
    # Note: only include calculations that aren't dependent on psi here
    cero.treated <- cero[which(cero$ever_treat==1),]
    
      # Calculating censoring weights
      ## Note:
      ## Censoring weight for time=1 is cumulative product from time=1 to time=k for ID i
      ## Censoring weight for time=2 is cumulative product from time=2 to time=k for ID i
      ## ...etc.
      tpw.treated <- 1/(cero.treated[["clost.pred"]]*cero.treated[["cdeath.pred"]]) # Weight for each time point
      cumw.treated <- unname(ave(tpw.treated, cero.treated[[id]], FUN=prod)) # Cumulative product of weights per ID
      lagtpw.treated <- suppressWarnings(unname(ave(tpw.treated, cero.treated[[id]], FUN=function(j) c(1, j[1:(length(j)-1)]), 1))) # Lagged weights by ID
      cumlagw.treated <- unname(ave(lagtpw.treated, cero.treated[[id]], FUN=cumprod)) # Cumulative product of lagged weights
      w.treated <- cumw.treated/cumlagw.treated
    
      # Contribution to estimating equation from H(psi)
      a.treated <- cero.treated[[x]]
      tcount.treated <- as.integer(table(cero.treated[[id]]))
      y.treated <- cero.treated[[y]]

    # Return environment
    dataprep.env <- list(
      d = d,
      n.followup = n.followup,
      newcov.untreated = newcov.untreated,
      newu1.untreated = newu1.untreated,
      cero.treated = cero.treated,
      w.treated = w.treated,
      a.treated = a.treated,
      tcount.treated = tcount.treated,
      y.treated = y.treated
    )
    
    return(dataprep.env)
  }
  dataprep.env <- estf.dataprep(data)

  # Calculations for estimating equation dependent on psi
  estf.iv <- function(psi, dataprep.results, blipfunction){
    hm <- rep(0, nrow(dataprep.results$cero.treated))
    last <- 0
    for (i in 1:length(dataprep.results$tcount.treated)){
      count <- dataprep.results$tcount.treated[i]
      start <- last+1
      last <- start+count-1
      atmp <- dataprep.results$a.treated[start:last]
      ytmp <- dataprep.results$y.treated[start:last]
      hm.id <- rep(0, count) 
      for (m in 1:count){
        for (k in count:dataprep.results$n.followup){
          sumblip <- 0
          for (j in m:k){
            if (j <= count){
              if (blipfunction==1){
                numtmp <- (k+1)-j
                denomtmp <- k-j+exp(psi*atmp[j])
              }
              if (blipfunction==2){
                numtmp <- 1
                denomtmp <- exp(psi*atmp[j])
              }
              blip.tmp2 <- log(numtmp)-log(denomtmp)
              sumblip <- sumblip + blip.tmp2
            }
          }
          exp.blip <- exp(sumblip)
          if (k>=count){
            hm.id[m] <- hm.id[m] + exp.blip
          }
        }
      }
      hm[start:last] <- hm.id
    }
    # Multiply Hm with weights and Z-E[Z]    
    u0i.treated <- hm*dataprep.results$w.treated*(dataprep.results$cero.treated[[z]]-dataprep.results$cero.treated[["z.pred"]])
    
    # Covariance matrix
    newu1.treated <- aggregate(u0i.treated, list(dataprep.results$cero.treated[[id]]), FUN=sum)[,2]
    newcov.treated <- newu1.treated %*% newu1.treated
    
    # Estimating Equation
    newcov <- newcov.treated + dataprep.results$newcov.untreated
    # newcov <- ifelse(newcov==0, 1e-12, newcov)
    newu <- sum(newu1.treated) + sum(dataprep.results$newu1.untreated)
    # return(newu %*% solve(newcov) %*% newu)
    return(as.numeric(newu*newu/newcov))
  }

  # Finding minimum of estimating equation
  psi <- NULL
  psi.esteq <- NULL
  psi.converge <- NULL
  estf.results <- suppressWarnings(optimx::optimx(start.value, estf.iv, dataprep.results=dataprep.env, blipfunction=blipfunction, method=c("nlminb")))
  if (estf.results$convcode!=0|estf.results$value>0.0001){
    estf.results1 <- suppressWarnings(optimx::optimx(start.value, estf.iv, dataprep.results=dataprep.env, blipfunction=blipfunction, method=c("nlm")))
    if (estf.results1$value<estf.results$value){
      psi <- estf.results1$p1
      psi.esteq <- estf.results1$value
      psi.converge <- estf.results1$convcode
    }else{
      psi <- estf.results$p1
      psi.esteq <- estf.results$value
      psi.converge <- estf.results$convcode
    }
  } else{
    psi <- estf.results$p1
    psi.esteq <- estf.results$value
    psi.converge <- estf.results$convcode
  }  
  
  results <- list(psi=psi,
                  psi.esteq=psi.esteq,
                  psi.converge=psi.converge)

  #Estimating equation across range of psi values
  if (grid==T){
    if (parallel==T){
      cl <- makeCluster(numCores)
      clusterExport(cl, ls(), envir=environment())
      est.eq.results <- parLapply(cl, seq(-grid.range, grid.range, by=grid.increment), function(i) {estf.iv(i, dataprep.env, blipfunction)})
      stopCluster(cl)
      psi.grid <- data.frame(cbind(psi=seq(-grid.range, grid.range, by=grid.increment), est.eq=do.call(rbind, est.eq.results)))
      }
    if (parallel==F){
      psi.grid <- data.frame(cbind(psi=seq(-grid.range, grid.range, by=grid.increment),
        est.eq=sapply(seq(-grid.range, grid.range, by=grid.increment), function(i){estf.iv(i, dataprep.env, blipfunction)})))
    }
  results[["psi.grid"]] <- psi.grid
  }
  
  # Function for blipping down/up
  blipf <- function(dataprep.results, psi.estimate){
    if (is.null(blip.data)==T){
      blip.d <- dataprep.results$d
    } else{
      blip.d <- blip.data
      if (is.null(clost)==F & is.null(clost.modelvars)==F){
        clost.model <- glm(as.formula(paste(clost, "==0~", paste(clost.modelvars)[2], sep="")), family=binomial(), data=blip.d)
        blip.d$clost.pred <- predict(clost.model, blip.d, type="response")
        blip.d$clost <- blip.d[,clost]
      } else{
        blip.d$clost.pred <- 1
        blip.d$clost <- 0
      }
      if (is.null(cdeath)==F & is.null(cdeath.modelvars)==F){
        cdeath.model <- glm(as.formula(paste(cdeath, "==0~", paste(cdeath.modelvars)[2], sep="")), family=binomial(), data=blip.d)
        blip.d$cdeath.pred <- predict(cdeath.model, blip.d, type="response")
        blip.d$cdeath <- blip.d[,cdeath]
      } else{
        blip.d$cdeath.pred <- 1
        blip.d$cdeath <- 0
      }
      blip.d <- blip.d[order(blip.d[,id], blip.d[,time]),]    
    }    
    
    blipdown <- expand.grid(unique(blip.d[[id]]), unique(blip.d[[time]]))
    colnames(blipdown) <- c(id, time)
    blipdown <- merge(x=as.data.frame(blip.d), y=blipdown, by=c(id, time), all=T)
    blipdown$ever_y <- ifelse(blipdown[[id]] %in% unique(blipdown[which(blipdown[[y]]==1),][[id]]), 1, 0)
    blipdown$ever_clost <- ifelse(blipdown[[id]] %in% unique(blipdown[which(blipdown$clost==1),][[id]]), 1, 0)
    blipdown$ever_cdeath <- ifelse(blipdown[[id]] %in% unique(blipdown[which(blipdown$cdeath==1),][[id]]), 1, 0)
    blipdown$ever_cens <- ifelse(blipdown$ever_clost==1|blipdown$ever_cdeath==1, 1, 0)
    blipdown[,y] <- ifelse(is.na(blipdown[[y]]) & blipdown$ever_y==1, 1, blipdown[,y])
    blipdown$ipcw <- ave(1/(blipdown$cdeath.pred*blipdown$clost.pred), blipdown[[id]], FUN=cumprod)
    blipdown$ipcw_max <- ave(blipdown$ipcw, blipdown[[id]], FUN=function(i) max(i, na.rm=T))
    blipdown$ipcw <- ifelse(is.na(blipdown$ipcw), blipdown$ipcw_max, blipdown$ipcw)
    blipdown <- blipdown[!is.na(blipdown[[y]]),]
    blipdown[,x] <- ifelse(is.na(blipdown[[x]]), 0, blipdown[[x]])
    y.0 <- paste(y, 0, sep = ".")
    blipdown[, y.0] <- NA
    
    for (t in unique(blip.d[[time]])){
      var <- paste("blip0", t, sep="_")
      varcum <- paste("blip0cum", t, sep="_")
      if (blipfunction==1){
        blipdown[,var] <- ifelse(blipdown[,time]>t, 1, ifelse(blipdown$ever_cens==0, (t+1-blipdown[[time]])/(t-blipdown[[time]]+exp(psi.estimate*blipdown[[x]])),0))
      }
      if (blipfunction==2){
        blipdown[,var] <- ifelse(blipdown[,time]>t, 1, ifelse(blipdown$ever_cens==0, 1/exp(psi.estimate*blipdown[,x]), 0))
      }
      blipdown[,varcum] <- ave(blipdown[[var]], blipdown[[id]], FUN=cumprod)
      blipdown[,varcum] <- ifelse(blipdown[[time]]>t, 1, blipdown[[varcum]])
      blipdown[, y.0] <- ifelse(blipdown[[time]]==t, blipdown[[y]]*blipdown[[varcum]], blipdown[[y.0]])
    }
    
    meanY <- data.frame(1:dataprep.results$n.followup)
    colnames(meanY) <- time
    meanY[,y] <- sapply(split(blipdown, blipdown[[time]]), function(i) weighted.mean(i[[y]], i[["ipcw"]]))
    meanY[,y.0] <- sapply(split(blipdown, blipdown[[time]]), function(i) weighted.mean(i[[y.0]], i[["ipcw"]]))
    
    Y0avgs <- rev(meanY[[y.0]])
    
    # Blipping Up
    if (is.null(intervention.regime)==T){
      regime <- rep(1, dataprep.results$n.followup)
    } else{
      regime <- intervention.regime
    }
    for (k in 1:dataprep.results$n.followup){
      matrix.tmp <- matrix(0, nrow=k, ncol=k)
      km1 <- k-1
      km3 <- k-3
      if (blipfunction==1|blipfunction==2){
        matrix.tmp[1,k] <- exp(psi.estimate*regime[k])
      }
      if (km3>=-1){
        for (c in seq(km3, -1, by=-1)){
          j <- c+2
          jp1 <- j+1
          cp1 <- c+1
          for (rowin in 0:km1){
            m <- rowin+1
            if (rowin<k-j){
              if (blipfunction==1){
                matrix.tmp[m,j] <- matrix.tmp[m,jp1]*((exp(psi.estimate*regime[jp1])-1)/(k-rowin-1-c)+1)
              }
              if (blipfunction==2){
                matrix.tmp[m,j] <- matrix.tmp[m,jp1]*exp(psi.estimate*regime[jp1])
              }
            } 
            if (rowin==(k-j)){
              sumcol <- sum(matrix.tmp[,jp1])
              if (blipfunction==1|blipfunction==2){
                matrix.tmp[m,j] <- (1-sumcol)*exp(psi.estimate*regime[jp1])
              }
            }
          }
        }
      }
      assign(paste("tg", k, sep="_"), matrix.tmp)
    }
    
    EYgs <- rep(0, dataprep.results$n.followup)
    for (matin in 1:dataprep.results$n.followup){
      EYgs[matin] <- Y0avgs[(dataprep.results$n.followup-matin+1):dataprep.results$n.followup] %*% get(paste("tg", matin, sep="_"))[,1]
    }
    meanY[,paste(y, "g", sep=".")] <- EYgs
    return(meanY)
  }
  
  # Checking that we found a minimum that is equal to zero and algorithm converged
  if (blipupdown==T){
    if ((blipfunction==1|blipfunction==2) & (psi.converge!=0|psi.esteq>0.0001)){
      print("Algorithm did not converge, and/or psi estimated on boundary. Cumulative risks under interventions will not be calculated.")
    } else{
      blip.results <- blipf(dataprep.env, psi)
      results[["blip.results"]] <- blip.results
    }
  }
  
  # Bootstrapping
  sncftm.boot <- function(i){
    # Bootstrapping sample
    set.seed(seed.vector[i])
    d.boot <- data.frame(original.id=sample(unique(data[[id]]), replace=T),
                         new.id=1:length(unique(data[[id]])))
    d.boot <- merge(d.boot, data, by.x="original.id", by.y=id, all.x=T)
    colnames(d.boot)[2] <- id
    # Data prep
    dataprep.boot <- estf.dataprep(d.boot)
    # Finding minimum
    psiboot <- NULL
    psiboot.esteq <- NULL
    psiboot.converge <- NULL
    estf.bootresults <- suppressWarnings(optimx::optimx(start.value, estf.iv, dataprep.results=dataprep.boot, blipfunction=blipfunction, method=c("nlminb")))
    if (estf.bootresults$convcode!=0|estf.bootresults$value>0.0001){
      estf.bootresults1 <- suppressWarnings(optimx::optimx(start.value, estf.iv, dataprep.results=dataprep.boot, blipfunction=blipfunction, method=c("nlm")))
      if (estf.bootresults1$value<estf.bootresults$value){
        psiboot <- estf.bootresults1$p1
        psiboot.esteq <- estf.bootresults1$value
        psiboot.converge <- estf.bootresults1$convcode
      }else{
        psiboot <- estf.bootresults$p1
        psiboot.esteq <- estf.bootresults$value
        psiboot.converge <- estf.bootresults$convcode
      }
    } else{
      psiboot <- estf.bootresults$p1
      psiboot.esteq <- estf.bootresults$value
      psiboot.converge <- estf.bootresults$convcode
    }
    # Blipping down/up
    if (blipupdown==T & psiboot.converge==0 & psiboot.esteq<0.0001){
      blip.results <- blipf(dataprep.boot, psiboot)
      results <- c(psiboot, psiboot.esteq, psiboot.converge, 
                   blip.results[,2], blip.results[,3], blip.results[,4])
      names(results) <- c("psi", "psi.esteq", "psi.converge",
                          paste("Y.t", blip.results[,1], sep=""),
                          paste("Y0.t", blip.results[,1], sep=""),
                          paste("Yg.t", blip.results[,1], sep=""))
      return(results)
    } else{
      results <- c(psi=psiboot, psi.esteq=psiboot.esteq, psi.converge=psiboot.converge)
      return(results)
    }
  }
  if (boot==T){
    set.seed(seed)
    seed.vector <- round(runif(R, min=0, max=1)*10000)
    if (parallel==T){
      cl <- makeCluster(numCores)
      clusterEvalQ(cl, library(optimx))
      clusterExport(cl, ls(), envir=environment())
      boot.results <- parLapply(cl, 1:R, function(i){sncftm.boot(i)})
      stopCluster(cl)
      boot.results <- do.call(rbind, boot.results)
    }
    if (parallel==F){
      boot.results <- lapply(1:R, function(i){sncftm.boot(i)})
      boot.results <- do.call(rbind, boot.results)
    }
  results[["boot.results"]] <- boot.results  
  }
  
  # Returning results
  return(results)
}
