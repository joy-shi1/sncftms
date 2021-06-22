# ------ G-ESTIMATION OF STRUCTUARL NESTED CUMULATIVE FAILURE TIME MODELS ------
# --------------------------------- SIMULATIONS --------------------------------
#
# -------------------------------- BY: Joy Shi ---------------------------------
# -------------------------- LAST MODIFIED: 2021-06-21 -------------------------
#
# Functions for simulating time-to-event data and analyzing the data using 
# structural nested cumulative failure time models. 
# 
# For all functions, arguments include:
#  - n: number of participants
#  - b: number of simulated data sets
#  - hazard: baseline hazard of the outcome (when treatment=0, confounder=0)
#  - seed: set seed value for simulations
#  - analysis: options include 
#         "unadj" for unadjusted SNCFTM
#         "adj" for confounding adjusted SNCFTM
#         "iv" for instrumental variable SNCFTM
#         "iv-cc" to simulate selecting a nested case-control sample in the 
#                 full simulated cohort and running an instrumental variable 
#                 SNCFTM in the case-control sample; marginal risks are
#                 calculated in the full simulated cohort
#         "g-formula" for confounding adjustment using the g-formula
#                     (note: only marginal risks are calculated; psi is not
#                     calculated)
#  - beta_ay: magnitude of the effect of A_k on Y_(k+1) at all time points
# ------------------------------------------------------------------------------

# ------------------ Loading and installing relevant packages ------------------

if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('optimx')) install.packages('optimx'); library('optimx')
if (!require('gfoRmula')) install.packages('gfoRmula'); library('gfoRmula')

# --------------------------- Simulation parameters ----------------------------

b <- 1000

# -------------------------- Loading SNCFTM functions --------------------------

source('./SNCFTM Function/SNCFTM Confounding Function.R')
source('./SNCFTM Function/SNCFTM IV Function.R')

# --------------------- Function for simulations for DAG 1 ---------------------
# - Time-fixed instrument Z
# - Time-fixed baseline confounder U
# - Time-fixed exposure A
# - Time-fixed outcome Y
# ------------------------------------------------------------------------------

sncftm.sim1 <- function(n=25000, b, hazard, seed=76245, analysis, beta_ay=0){
  
  # Creating list of seeds
  set.seed(seed)
  seed.list <- round(runif(b)*100000)
  
  # Looping through b iterations
  results <- do.call(rbind, lapply(1:b, function(i){  
    
    set.seed(seed.list[i])
    
    # Generating data
    id <- seq(1:n)
    l <- rbinom(n, size=1, prob=0.25)
    z <- rbinom(n, size=1, prob=0.5)
    a <- rbinom(n, size=1, prob=(0.25*z+0.5*l))
    logodds.y <- log(hazard/(1-hazard))+0.5*l+beta_ay*a
    prob.y <- 1/(1+exp(-logodds.y))
    y <- rbinom(n, size=1, p=prob.y)
    time <- rep(1, n)
    sim.data <- data.frame(id, time, l, z, a, y)
    
    # Generating case-control sample
    if (analysis=="iv-cc"){
      
      # Select controls
      id.list <- sample(sim.data[which(sim.data$y==0),"id"], size=sum(sim.data$y)*2, replace=T)
      sim.controls <- data.frame(originalid=id.list, newid=1:length(id.list))
      sim.controls <- merge(sim.controls, sim.data, by.x="originalid", by.y="id")
      sim.controls$z.include <- 1
      
      # Adding cases
      sim.cases <- sim.data[which(sim.data$y==1),]
      names(sim.cases)[1] <- "originalid"
      sim.cases$newid <- (nrow(sim.controls)+1):(nrow(sim.controls)+nrow(sim.cases))
      sim.cases$z.include <- 0
      
      # Combining data
      sim.ccdata <- rbind(sim.controls, sim.cases)
      
      # Renaming ID variable of original simulated data
      colnames(sim.data)[1] <- "newid"
      
    }

    # Analysis
    if (analysis=="unadj"){
      results.unadj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
        x.modelvars=~1, x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
        parallel=F)
      
      results2 <- c(results.unadj$psi, results.unadj$psi.esteq, results.unadj$blip.results[,3], results.unadj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))
      
    } else if (analysis=="adj"){
      results.adj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
        x.modelvars=~l, x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
        parallel=F)
      
      results2 <- c(results.adj$psi, results.adj$psi.esteq, results.adj$blip.results[,3], results.adj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""),paste("yg.t", paste(1:max(sim.data$time)), sep=""))  
      
    } else if (analysis=="iv"){
      results.iv <- sncftm.iv(data=sim.data, id="id", time="time", z="z", z.modelvars=~1,
         z.family="binomial", x="a", y="y", blipfunction=1, boot=F, parallel=F)  
      
      results2 <- c(results.iv$psi, results.iv$psi.esteq, results.iv$blip.results[,3], results.iv$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))  
      
    } else if (analysis=="iv-cc"){
      results.ivcc <- sncftm.iv(data=sim.ccdata, id="newid", time="time", z="z", z.modelvars=~1,
        z.family="binomial", x="a", y="y", blipfunction=1, blip.data=sim.data, boot=F, parallel=F,
        z.indicator="z.include")
      
      results.pre <- c(results.ivcc$psi, results.ivcc$psi.esteq, results.ivcc$blip.results[,3], results.ivcc$blip.results[,4])
      if (length(results.pre)<4){
        results2 <- rep(NA,4)
        results2[1:length(results.pre)] <- results.pre
      } else{
        results2 <- results.pre
      }
      
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))  
      
    } else if (analysis=="g-formula"){
      y.model <- glm(y~a+l+a:l, data=sim.data)
      sim.data0 <- sim.data 
      sim.data0$a <- 0
      sim.data0$pred.y <- predict(y.model, newdata=sim.data0)
      sim.data1 <- sim.data 
      sim.data1$a <- 1
      sim.data1$pred.y <- predict(y.model, newdata=sim.data1) 
      results2 <- c(y0.t1=mean(sim.data0$pred.y), yg.t1=mean(sim.data1$pred.y))
    }
    return(results2)
  }))
  
  # Returning results
  return(results)
}    

# --------------------- Function for simulations for DAG 2 ---------------------
# - Time-fixed instrument Z
# - Time-fixed baseline confounder U
# - Time-fixed exposure A
# - Time-varying outcome Y (with two time points)
# ------------------------------------------------------------------------------

sncftm.sim2 <- function(n=25000, b, hazard, seed=84756, analysis, beta_ay=0){
  
  # Creating list of seeds
  set.seed(seed)
  seed.list <- round(runif(b)*100000)
  
  # Looping through b iterations
  results <- do.call(rbind, lapply(1:b, function(i){  

    set.seed(seed.list[i])
        
    # Generating data
    id <- seq(1:n)
    l <- rbinom(n, size=1, prob=0.25)
    z <- rbinom(n, size=1, prob=0.5)
    a <- rbinom(n, size=1, prob=(0.25*z+0.5*l))
    logodds.y1 <- log(hazard/(1-hazard))+0.5*l+beta_ay*a
    logodds.y2 <- log(hazard/(1-hazard))+0.5*l+beta_ay*a
    prob.y1 <- 1/(1+exp(-logodds.y1))
    prob.y2 <- 1/(1+exp(-logodds.y2))
    y1 <- rbinom(n, size=1, p=prob.y1)
    y2 <- rbinom(n, size=1, p=prob.y2)
    y2 <- ifelse(y1==1, NA, y2)
    sim.data <- rbind(data.frame(id, time=rep(1,n), l, z, a, y=y1),
                      data.frame(id, time=rep(2,n), l, z, a, y=y2))
    sim.data <- sim.data[which(!is.na(sim.data$y)),]
    sim.data <- sim.data[order(sim.data$id, sim.data$time),]
    
    
    # Generating case-control sample
    if (analysis=="iv-cc"){
      
      # Selecting controls
      id.list1 <- sample(id[which(y1==0)], size=sum(y1)*2, replace=T)
      id.list2 <- sample(id[which(y1==0 & y2==0)], size=sum(y1==0 & y2==1)*2, replace=T)
      sim.controls <- data.frame(originalid=c(id.list1, id.list2), newid=1:(length(id.list1)+length(id.list2)))
      sim.controls <- merge(sim.controls, sim.data, by.x="originalid", by.y="id")
      sim.controls$z.include <- 1
      
      # Adding cases
      sim.cases <- sim.data[which(sim.data$id %in% id[which(y1==1|y2==1)]),]
      names(sim.cases)[1] <- "originalid"
      sim.cases$newid <- (max(sim.controls$newid)+1):(max(sim.controls$newid)+nrow(sim.cases))
      sim.cases$newid <- ave(sim.cases[['newid']], sim.cases[['originalid']], FUN=min)
      sim.cases$z.include <- 0
      
      # Combining data
      sim.ccdata <- rbind(sim.controls, sim.cases)
      
      # Renaming ID variable of original simulated data
      colnames(sim.data)[1] <- "newid"
             
    }
    
    # Analysis
    if (analysis=="unadj"){
      results.unadj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
                                   x.modelvars=~as.factor(time), x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
                                   parallel=F)
      
      results2 <- c(results.unadj$psi, results.unadj$psi.esteq, results.unadj$blip.results[,3], results.unadj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""),paste("yg.t", paste(1:max(sim.data$time)), sep=""))
      
    } else if (analysis=="adj"){
      results.adj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
                                 x.modelvars=~l*as.factor(time), x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
                                 parallel=F)
      
      results2 <- c(results.adj$psi, results.adj$psi.esteq, results.adj$blip.results[,3], results.adj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))    
      
    } else if (analysis=="iv"){
      results.iv <- sncftm.iv(data=sim.data, id="id", time="time", z="z", z.modelvars=~1,
                              z.family="binomial", x="a", y="y", blipfunction=1, boot=F, parallel=F)  
      
      results2 <- c(results.iv$psi, results.iv$psi.esteq, results.iv$blip.results[,3], results.iv$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep="")) 
      
    } else if (analysis=="iv-cc"){
      results.ivcc <- sncftm.iv(data=sim.ccdata, id="newid", time="time", z="z", z.modelvars=~1,
                                z.family="binomial", x="a", y="y", blipfunction=1, blip.data=sim.data, boot=F, parallel=F,
                                z.indicator="z.include")
      
      results.pre <- c(results.ivcc$psi, results.ivcc$psi.esteq, results.ivcc$blip.results[,3], results.ivcc$blip.results[,4])
      if (length(results.pre)<6){
        results2 <- rep(NA,6)
        results2[1:length(results.pre)] <- results.pre
      } else{
        results2 <- results.pre
      }
      
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))  
   
    } else if (analysis=="g-formula"){

      y.model <- glm(y~a*l*as.factor(time), data=sim.data)
      
      # Creating expanded dataset
      sim.expanded <- rbind(data.frame(id, time=rep(1,n), l, z, a),
                            data.frame(id, time=rep(2,n), l, z, a))     
      
      # Estimates for y0
      sim.data0 <- sim.expanded 
      sim.data0$a <- 0
      sim.data0$pred.y <- predict(y.model, newdata=sim.data0)
      sim.data0$pred.survival <- ave(1-sim.data0$pred.y, sim.data0$id, FUN=cumprod)
      sim.data0$cum.y <- 1-sim.data0$pred.survival
      
      # Estimates for yg
      sim.data1 <- sim.expanded 
      sim.data1$a <- 1
      sim.data1$pred.y <- predict(y.model, newdata=sim.data1)
      sim.data1$pred.survival <- ave(1-sim.data1$pred.y, sim.data1$id, FUN=cumprod)
      sim.data1$cum.y <- 1-sim.data1$pred.survival
      
      results2 <- c(aggregate(sim.data0$cum.y, list(sim.data0$time), FUN=mean)[,2],
                    aggregate(sim.data1$cum.y, list(sim.data1$time), FUN=mean)[,2])
      
      names(results2) <- c(paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep="")) 
    }
    return(results2)
  }))
  
  # Returning results
  return(results)
}    

# --------------------- Function for simulations for DAG 3 ---------------------
# - Time-fixed instrument Z
# - Time-varying confounder U (with two time points)
# - Time-varying exposure A
# - Time-varying outcome Y 
# ------------------------------------------------------------------------------

sncftm.sim3 <- function(n=25000, b=1000, hazard, seed=84756, analysis, beta_ay=0){
  
  # Creating list of seeds
  set.seed(seed)
  seed.list <- round(runif(b)*100000)
  
  # Looping through b iterations
  results <- do.call(rbind, lapply(1:b, function(i){  
    
    set.seed(seed.list[i])
    
    # Generating data
    id <- seq(1:n)
    l1 <- rbinom(n, size=1, prob=0.25)
    z <- rbinom(n, size=1, prob=0.5)
    a1 <- rbinom(n, size=1, prob=(0.25*z+0.5*l1))
    u <- rbinom(n, size=1, prob=0.25)
    l2 <- rbinom(n, size=1, prob=(0.1+0.25*a1+0.25*l1+0.25*u))
    a2 <- rbinom(n, size=1, prob=(0.25*z+0.5*l2))
    logodds.y1 <- log(hazard/(1-hazard))+0.5*l1+beta_ay*a1
    logodds.y2 <- log(hazard/(1-hazard))+0.5*u+beta_ay*a2
    prob.y1 <- 1/(1+exp(-logodds.y1))
    prob.y2 <- 1/(1+exp(-logodds.y2))
    y1 <- rbinom(n, size=1, p=prob.y1)
    y2 <- rbinom(n, size=1, p=prob.y2)
    y2 <- ifelse(y1==1, NA, y2)
    sim.data <- rbind(data.frame(id, time=rep(1,n), l=l1, z, a=a1, y=y1),
                      data.frame(id, time=rep(2,n), l=l2, z, a=a2, y=y2))
    sim.data <- sim.data[which(!is.na(sim.data$y)),]    
    
    # Generating case-control sample
    if (analysis=="iv-cc"){
      
      # Selecting controls
      id.list1 <- sample(id[which(y1==0)], size=sum(y1)*2, replace=T)
      id.list2 <- sample(id[which(y1==0 & y2==0)], size=sum(y1==0 & y2==1)*2, replace=T)
      sim.controls <- data.frame(originalid=c(id.list1, id.list2), newid=1:(length(id.list1)+length(id.list2)))
      sim.controls <- merge(sim.controls, sim.data, by.x="originalid", by.y="id")
      sim.controls$z.include <- 1
      
      # Adding cases
      sim.cases <- sim.data[which(sim.data$id %in% id[which(y1==1|y2==1)]),]
      names(sim.cases)[1] <- "originalid"
      sim.cases$newid <- (max(sim.controls$newid)+1):(max(sim.controls$newid)+nrow(sim.cases))
      sim.cases$newid <- ave(sim.cases[['newid']], sim.cases[['originalid']], FUN=min)
      sim.cases$z.include <- 0
      
      # Combining data
      sim.ccdata <- rbind(sim.controls, sim.cases)
      
      # Renaming ID variable of original simulated data
      colnames(sim.data)[1] <- "newid"      
    }
    
    # Analysis
    if (analysis=="unadj"){
      results.unadj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
                                   x.modelvars=~as.factor(time), x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
                                   parallel=F)
      
      results2 <- c(results.unadj$psi, results.unadj$psi.esteq, results.unadj$blip.results[,3], results.unadj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""),paste("yg.t", paste(1:max(sim.data$time)), sep=""))
      
    } else if (analysis=="adj"){
      results.adj <- sncftm.conf(data=sim.data, id="id", time="time", x="a", 
                                 x.modelvars=~l*as.factor(time), x.linkfunction="logit", y="y", blipfunction=1, boot=F, 
                                 parallel=F)
      
      results2 <- c(results.adj$psi, results.adj$psi.esteq, results.adj$blip.results[,3], results.adj$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))    
      
    } else if (analysis=="iv"){
      results.iv <- sncftm.iv(data=sim.data, id="id", time="time", z="z", z.modelvars=~1,
                              z.family="binomial", x="a", y="y", blipfunction=1, boot=F, parallel=F)  
      
      results2 <- c(results.iv$psi, results.iv$psi.esteq, results.iv$blip.results[,3], results.iv$blip.results[,4])
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep="")) 
      
    } else if (analysis=="iv-cc"){
      results.ivcc <- sncftm.iv(data=sim.ccdata, id="newid", time="time", z="z", z.modelvars=~1,
                                z.family="binomial", x="a", y="y", blipfunction=1, blip.data=sim.data, boot=F, parallel=F,
                                z.indicator="z.include")
    
      results.pre <- c(results.ivcc$psi, results.ivcc$psi.esteq, results.ivcc$blip.results[,3], results.ivcc$blip.results[,4])
      if (length(results.pre)<6){
        results2 <- rep(NA,6)
        results2[1:length(results.pre)] <- results.pre
      } else{
        results2 <- results.pre
      }      
      
      names(results2) <- c("psi", "esteq", paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep=""))  
      
    } else if (analysis=="g-formula"){
      
      gformula.data <- sim.data
      gformula.data$time <- gformula.data$time-1

      interventions <- list(list(c(static, rep(0, 2))), list(c(static, rep(1, 2))))
      
      results.gformula <- gformula(gformula.data, id="id", time_points=2, time_name="time",
                       covnames=c("l","a"), outcome_name="y", outcome_type="survival",
                       covtypes=c("binary", "binary"), histories=c(lagged, lagged),
                       histvars=list(c("a", "l"), c("l")),
                       covparams=list(covmodels=c(l~lag1_a + lag1_l, a~lag1_l)),
                       ymodel=y~a*l*time,
                       intvars=list('a', 'a'),
                       interventions=interventions,
                       int_descript <- c('Never treat', 'Always treat'),
                       seed=123,
                       model_fits=F)
      
      interventions <- list(list(c(static, rep(0, 2))),
                            list(c(static, rep(1, 2))))

      results2 <- as.numeric(c(results.gformula$result[2,4], results.gformula$result[5,4],
                               results.gformula$result[3,4], results.gformula$result[6,4]))
  
      names(results2) <- c(paste("y0.t", paste(1:max(sim.data$time)), sep=""), paste("yg.t", paste(1:max(sim.data$time)), sep="")) 
    }
    return(results2)
  }))
  
  # Returning results
  return(results)
}    

# ---------------------- Example execution of simulations ----------------------
sim.results1 <- sncftm.sim1(hazard=1/20, b=b, analysis="adj", beta_ay=0)
sim.results2 <- sncftm.sim2(hazard=1/20, b=b, analysis="iv", beta_ay=0.5)
sim.results3 <- sncftm.sim3(hazard=1/10, b=b, analysis="iv-cc", beta_ay=0)