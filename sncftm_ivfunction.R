# data=d
# data$alcohol <- -data$alcohol
# id="id"
# time="time"
# z="z"
# z.model=z.model
# x="alcohol"
# y="event"
# center=0
# blipfunction=1
# blipupdown=T
# grid=T
# clost <- NULL
# clost.model <- NULL
# cdeath <- NULL
# cdeath.model <- NULL

sncftm.iv <- function(data, id, time, z, z.model, x, y, 
                      clost=NULL, clost.model=NULL, cdeath=NULL, cdeath.model=NULL, 
                      grid=F, center=0, iteration=6, blipfunction, blipupdown=T){
  
  # Data Prep
  z.predicted <- predict(z.model, newdata=data, type="response")
  if (is.null(clost)==F & is.null(clost.model)==F){
    clost.predicted <- predict(clost.model, data, type="response")
    data <- data %>% mutate(clost=!!as.name(clost))
  } else{
    clost.predicted <- rep(1, nrow(data))
    data <- data %>% mutate(clost=0)
  }
  if (is.null(cdeath)==F & is.null(cdeath.model)==F){
    cdeath.predicted <- predict(cdeath.model, data, type="response")
    data <- data %>% mutate(cdeath=!!as.name(cdeath))
  } else{
    cdeath.predicted <- rep(1, nrow(data))
    data <- data %>% mutate(cdeath=0)
  }
  
  data <- data %>% 
    mutate(z.pred=z.predicted) %>%
    mutate(clost.pred=clost.predicted) %>%
    mutate(cdeath.pred=cdeath.predicted) %>%
    mutate(c.pred=clost.pred*cdeath.pred) %>%
    arrange(!!as.name(id), !!as.name(time)) 
  
  n.followup <- data %>% distinct(!!as.name(time)) %>% nrow()
  min.time <- min(data[,time], na.rm=T)
  max.time <- max(data[,time], na.rm=T)
  
  # Estimating Function for SNCFTM
  estf.iv <- function(psi, blipfunction){
    cero <- data %>% 
      group_by(!!as.name(id)) %>% 
      mutate(ever_y=suppressWarnings(max(!!as.name(y), na.rm=T))) %>% 
      filter(ever_y==1) %>% 
      dplyr::select(-ever_y) %>% 
      mutate(count=n()) %>% 
      mutate(ever_treat=max(!!as.name(x), na.rm=T)) %>%
      filter(!is.na(!!as.name(z)))
    
    # Untreated
    cero.untreated <- cero %>% filter(ever_treat==0)
    if (nrow(cero.untreated)==0){
      newcov.untreated <- 0
      newu1.untreated <- 0
    } else{
      predc.untreated <- cbind(cero.untreated[["clost.pred"]], cero.untreated[["cdeath.pred"]])
      residz.untreated <- cero.untreated[[z]]-cero.untreated[["z.pred"]]
      tcount.untreated <- cero.untreated %>% distinct(!!as.name(id), count) %>% pull(count)
      y.untreated <- cero.untreated[[y]]
      time.untreated <- cero.untreated[[time]]
      
      hm.untreated <- rep(0, nrow=nrow(cero.untreated))
      last <- 0
      for (i in 1:length(tcount.untreated)){
        count <- tcount.untreated[i]
        start <- last+1;
        last <- start+count-1;
        ytmp <- rep(1, n.followup)
        ytmp[1:length(y.untreated[start:last])] <- y.untreated[start:last]
        ctmp <- matrix(predc.untreated[start:last,], ncol=2)
        ipcwnew <- matrix(0, nrow=count, ncol=n.followup)
        for (m in 1:count){
          for (k in m:n.followup){
            for (j in m:k){
              if (j<=count){
                if (j==m){
                  tmpw <- ctmp[j,]
                } else{
                  tmpw <- tmpw*ctmp[j,]
                }            
              }
              ipcwsep <- as.matrix(1/tmpw)
              ipcwnew[m,k] <- apply(ipcwsep, 2, prod) 
            }
          }
          hm.id <- ipcwnew %*% ytmp
        }
        hm.untreated[start:last] <- hm.id
      }
      
      u0i.untreated <- hm.untreated * residz.untreated
      newu1.untreated <- cbind(cero.untreated[[id]], u0i.untreated) %>% data.frame()
      colnames(newu1.untreated) <- c("id", "u0i.untreated")
      newu1.untreated <- newu1.untreated %>% group_by(id) %>%summarise_all(sum) %>% pull(u0i.untreated)
      newcov.untreated <- newu1.untreated %*% newu1.untreated
    }
    
    # Treated
    cero.treated <- cero %>% filter(ever_treat!=0) %>% filter(!is.na(!!as.name(z)))
    predc.treated <- cbind(cero.treated[["clost.pred"]], cero.treated[["cdeath.pred"]])
    a.treated <- cero.treated[[x]]
    residz.treated <- cero.treated[[z]]-cero.treated[["z.pred"]]
    tcount.treated <- cero.treated %>% distinct(!!as.name(id), count) %>% pull(count)
    y.treated <- cero.treated[[y]]
    time.treated <- cero.treated[[time]]
    
    hm <- rep(0, nrow(cero.treated))
    last <- 0
    for (i in 1:length(tcount.treated)){
      count <- tcount.treated[i]
      start <- last+1
      last <- start+count-1
      ctmp <- matrix(predc.treated[start:last,], ncol=2)
      atmp <- a.treated[start:last]
      ytmp <- y.treated[start:last]
      hm.id <- rep(0, count) 
      ipcwnew <- matrix(0, nrow=count, ncol=n.followup)
      for (m in 1:count){
        for (k in m:n.followup){
          sumblip <- 0
          for (j in m:k){
            if (j <= count){
              if (j==m){
                tmpw=ctmp[j,]
              } else{
                tmpw=tmpw*ctmp[j,]
              }
              ipcwsep <- as.matrix(1/tmpw)
              ipcwnew <- apply(ipcwsep, 2, prod) 
              if (blipfunction==1){
                numtmp <- (k+1)-j
                denomtmp <- k-j+exp(psi*atmp[j])
              }
              if (blipfunction==2){
                numtmp <- 1
                denomtmp <- exp(psi*atmp[j])
              }
              if (blipfunction==3){
                numtmp <- 1
                denomtmp <- exp(psi[m]*atmp[j])
              }
              blip.tmp2 <- log(numtmp)-log(denomtmp)
              sumblip <- sumblip + blip.tmp2
            }
          }
          exp.blip <- exp(sumblip)
          if (k>=count){
            hm.id[m] <- hm.id[m] + ipcwnew*exp.blip
          }
        }
      }
      hm[start:last] <- hm.id
    }
    
    u0i.treated <- hm*residz.treated
    newu1.treated <- cbind(cero.treated[[id]], u0i.treated) %>% data.frame()
    colnames(newu1.treated) <- c("id", "u0i.treated")
    newu1.treated <- newu1.treated %>% group_by(id) %>%summarise_all(sum) %>% pull(u0i.treated)
    newcov.treated <- newu1.treated %*% newu1.treated
    
    # Estimating Equation
    newcov <- newcov.treated + newcov.untreated
    newu <- sum(newu1.treated) + sum(newu1.untreated)
    return(newu %*% solve(newcov) %*% newu)
  }
  
  # Finding minimum of estimating equation
  if(blipfunction==1|blipfunction==2){
    psi <- NULL
    psi.esteq <- NULL
    psi.converge <- NULL
    bound <- center
    i <- 1
    while(is.null(psi)){
      bound <- bound+0.25
      estf.results <- optimize(estf.iv, interval=c(-bound, bound), blipfunction=blipfunction)
      if ((bound-sum(abs(estf.results$minimum))>0.05 & estf.results$objective<0.0001)|i==iteration){
        psi <- estf.results$minimum
        psi.esteq <- estf.results$objective  
        psi.converge <- 0
      } else{
        i <- i+1
      }
    }
    
    # Estimating equation across range of psi values
    if (grid==T){
      psi.grid <- cbind(psi=seq(-(center+0.25*iteration), (center+0.25*iteration), by=0.03),
                        est.eq=sapply(seq(-(center+0.25*iteration), (center+0.25*iteration), by=0.03), function(i){estf.iv(i, blipfunction)})) %>%
        data.frame()
    } else{
      psi.grid <- NULL
    }
    
  }
  
  if (blipfunction==3){
    estf.results <- optim(rep(0, n.followup), estf.iv, blipfunction=blipfunction, )
    psi <- estf.results$par
    psi.esteq <- estf.results$value
    psi.converge <- estf.results$convergence
    psi.grid <- NULL
  }
  
  # Checking psi is on boundary; if not, then blip down and up
  if (blipupdown==F){
    return(list(psi.grid=psi.grid, psi=psi, minimum=psi.esteq, counterf.risks=NULL))
  }
  if (blipupdown==T){
    if ((blipfunction==1|blipfunction==2) & (bound-sum(abs(psi))<0.05)){
      print("Algorithm did not converge, and/or psi estimated on boundary. Cumulative risks under interventions will not be calculated.")
      list(psi.grid=psi.grid, psi=psi, minimum=psi.esteq, counterf.risks=NULL)
    } else if (blipfunction==3 & (psi.converge!=0)){
      print("Algorithm did not converge, and/or psi estimated on boundary. Cumulative risks under interventions will not be calculated.")
      list(psi.grid=psi.grid, psi=psi, minimum=psi.esteq, counterf.risks=NULL)
    } else{
      blipdown <- expand.grid(unique(data[[id]]), unique(data[[time]])) %>%
        mutate(!!id:=Var1) %>%
        mutate(!!time:=Var2) %>%
        full_join(data, by=c(id, time)) %>%
        group_by(!!as.name(id)) %>%
        mutate(ever_y=suppressWarnings(max(!!as.name(y), na.rm=T))) %>%
        mutate(ever_clost=max(clost, na.rm=T)) %>%
        mutate(ever_cdeath=max(cdeath, na.rm=T)) %>%
        mutate(ever_cens=ifelse(ever_clost==1|ever_cdeath==1, 1, 0)) %>%
        mutate(!!y:=ifelse(is.na(!!as.name(y)) & ever_y==1, 1, !!as.name(y))) %>%
        mutate(ipcw=cumprod(1/c.pred)) %>%
        fill(ipcw) %>%
        filter(!is.na(!!as.name(y))) %>%
        mutate(!!x:=ifelse(is.na(!!as.name(x)), 0, !!as.name(x))) %>%
        arrange(!!as.name(id), !!as.name(time)) %>%
        mutate(!!(paste(y, 0, sep=".")):=NA)
      
      for (t in unique(data[[time]])){
        var <- paste("blip0", t, sep="_")
        varcum <- paste("blip0cum", t, sep="_")
        if (blipfunction==1){
          blipdown <- blipdown %>% mutate(!!var:=if_else(!!as.name(time)>t, 1, if_else(ever_cens==0,
                                                                                       (t+1-!!as.name(time))/(t-!!as.name(time)+exp(psi*!!as.name(x))), 0)))
        }
        if (blipfunction==2){
          blipdown <- blipdown %>% mutate(!!var:=if_else(!!as.name(time)>t, 1, if_else(ever_cens==0, (1/exp(psi*!!as.name(x))), 0)))
        }
        if (blipfunction==3){
          blipdown <- blipdown %>% mutate(!!var:=if_else(!!as.name(time)>t, 1, if_else(ever_cens==0, (1/exp(psi[time]*!!as.name(x))), 0)))
        }      
        blipdown <- blipdown %>% mutate(!!varcum:=!!as.name(var))
        blipdown <- blipdown %>% group_by(id) %>% mutate(!!varcum:=cumprod(!!as.name(var)))
        blipdown <- blipdown %>% mutate(!!varcum:=ifelse(!!as.name(time)>t, 1, !!as.name(varcum)))
        blipdown <- blipdown %>% mutate(!!(paste(y, 0, sep=".")):=ifelse(!!as.name(time)==t, !!as.name(y)*!!as.name(varcum),
                                                                         !!as.name(paste(y, 0, sep="."))))
      }
      
      meanY <- blipdown %>% group_by(!!as.name(time)) %>% summarize_at(vars(!!as.name(y), !!as.name(paste(y, 0, sep="."))), 
                                                                       ~weighted.mean(., w=ipcw))
      Y0avgs <- meanY[[paste(y, 0, sep=".")]] %>% rev()
      
      # Blipping Up
      regime <- rep(1, n.followup)
      for (k in 1:n.followup){
        matrix.tmp <- matrix(0, nrow=k, ncol=k)
        km1 <- k-1
        km3 <- k-3
        if (blipfunction==1|blipfunction==2){
          matrix.tmp[1,k] <- exp(psi*regime[k])
        }
        if (blipfunction==3){
          matrix.tmp[1,k] <- exp(psi[k]*regime[k])
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
                  matrix.tmp[m,j] <- matrix.tmp[m,jp1]*((exp(psi*regime[jp1])-1)/(k-rowin-1-c)+1)
                }
                if (blipfunction==2){
                  matrix.tmp[m,j] <- matrix.tmp[m,jp1]*exp(psi*regime[jp1])
                }
                if (blipfunction==3){
                  matrix.tmp[m,j] <- matrix.tmp[m,jp1]*exp(psi[jp1]*regime[jp1])
                }
              } 
              if (rowin==(k-j)){
                sumcol <- sum(matrix.tmp[,jp1])
                if (blipfunction==1|blipfunction==2){
                  matrix.tmp[m,j] <- (1-sumcol)*exp(psi*regime[jp1])
                }
                if (blipfunction==3){
                  matrix.tmp[m,j] <- (1-sumcol)*exp(psi[jp1]*regime[jp1])
                }
              }
            }
          }
        }
        assign(paste("tg", k, sep="_"), matrix.tmp)
      }
      
      EYgs <- rep(0, n.followup)
      for (matin in 1:n.followup){
        EYgs[matin] <- Y0avgs[(n.followup-matin+1):n.followup] %*% get(paste("tg", matin, sep="_"))[,1]
      }
      meanY <- cbind(meanY, EYgs) %>% rename(!!paste(y, "g", sep="."):=EYgs)
      
      # Returning Values
      list(psi.grid=psi.grid, psi=psi, minimum=psi.esteq, counterf.risks=as.matrix(meanY))
    }
  }
}