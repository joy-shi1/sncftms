# ------ G-ESTIMATION OF STRUCTUARL NESTED CUMULATIVE FAILURE TIME MODELS: -----
# -------------------------------- DATA EXAMPLE --------------------------------
#
# -------------------------------- BY: Joy Shi ---------------------------------
# -------------------------- LAST MODIFIED: 2021-06-21 -------------------------
#

# ------------------------------- Data parameters ------------------------------
n <- 10000
set.seed(54432)
time.points <- 5


# ---------------------------- Generating variables ----------------------------

  # id
  id <- seq(1:n)
  
  # baseline u (note: can be time-varying)
  u <- rbinom(n, size=1, prob=0.25)
  
  # baseline z (note: can be time-varying)
  z <- rbinom(n, size=1, prob=0.5)
  
  # time-varying treatment a
  a <- rbinom(n*time.points, size=1, prob=(0.25*rep(z, time.points) + 0.5*rep(u, time.points)))
  
  # time-varying outcome y 
  y <- log(1/39)+0.3*rep(u, time.points)+0.5*rep(a, time.points)
  y <- rbinom(n*5, size=1, p=1/(1+exp(-y)))
  
  # censoring due to lost to follow-up
  censor.lost <- log(1/19) + 0.3*rep(u, time.points) + 0.3*a
  censor.lost <- rbinom(n*5, size=1, p=1/(1+exp(-censor.lost)))
  
  # censoring due to death
  censor.death <- log(1/19) + 0.4*rep(u, time.points) + 0.4*a
  censor.death <- rbinom(n*5, size=1, p=1/(1+exp(-censor.death)))
    
  # Combining into one dataset
  sim.data <- data.frame(id=rep(id, time.points),
                         time=rep(1:time.points, each=n),
                         u=rep(u, time.points),
                         z=rep(z, time.points),
                         a=a,
                         y=y,
                         censor.lost=censor.lost,
                         censor.death=censor.death)

  # Note: since this is simulated data, the number of rows for each participant
  # is equal to the number of time points simulated. In reality, there may be 
  # fewer rows per participant because follow-up ends earlier (i.e. due to 
  # developing the outcome or censoring). As such, we will remove observations
  # after the first instance of developing the outcome or censoring.

  # Create variable for time point where outcome/censoring occurred
  y.time <- ave(ifelse(y==0, 5, y*sim.data$time), id, FUN=min)
  censor.lost.time <- ave(ifelse(censor.lost==0, 5, censor.lost*sim.data$time), id, FUN=min)
  censor.death.time <- ave(ifelse(censor.death==0, 5, censor.death*sim.data$time), id, FUN=min)
  
  # Drop observations if after event or censoring
  sim.data <- sim.data[which(sim.data$time<=y.time & 
                             sim.data$time<=censor.lost.time &
                             sim.data$time<=censor.death.time),]
  
  # Only permit event or censoring (not both) to occur at any given time point
  # Arbitrarily prioritize y>censor.death>censor.lost
  sim.data$censor.lost <- ifelse(sim.data$y==1|sim.data$censor.death==1, 0, censor.lost)
  sim.data$censor.death <- ifelse(sim.data$y==1, 0, censor.death)
  
  # outcome should be missing if censored
  sim.data$y <- ifelse(sim.data$censor.lost==1|sim.data$censor.death==1, NA, sim.data$y)
    
  # Arrange by id and time
  sim.data <- sim.data[order(sim.data$id, sim.data$time),]

  # Final set-up of data:
  rownames(sim.data) <- NULL
  head(sim.data, n=20)
  
  
# ---------------------- Analysis of simulated data using ----------------------
# ------------------- SNCFTM with adjustment for confounding -------------------
  
sncftm.confresults <- sncftm.conf(data = sim.data,
                  id = "id",
                  time = "time",
                  x = "a",
                  x.modelvars = ~u+as.factor(time),
                  x.linkfunction="logit",
                  y = "y",
                  clost = "censor.lost",
                  clost.modelvars = ~a + u + as.factor(time),
                  cdeath = "censor.death",
                  cdeath.modelvars = ~a + u + as.factor(time),
                  blipfunction = 1,
                  start.value = 0,
                  grid=T,
                  grid.range=1.5,
                  grid.increment=0.01,
                  blipupdown=T,
                  boot=T,
                  R=10,
                  parallel=T)

# ---------------------- Analysis of simulated data using ----------------------
# ----------------------- instrumental variable analysis -----------------------

sncftm.ivresults <- sncftm.iv(data = sim.data,
                  id = "id",
                  time = "time",
                  z = "z",
                  z.modelvars = ~1,
                  z.family="gaussian",
                  x = "a",
                  y = "y",
                  clost = "censor.lost",
                  clost.modelvars = ~a + u + as.factor(time),
                  cdeath = "censor.death",
                  cdeath.modelvars = ~a + u + as.factor(time),
                  blipfunction = 1,
                  start.value = 0,
                  grid=T,
                  grid.range=1.5,
                  grid.increment=0.01,
                  blipupdown=T,
                  boot=T,
                  R=10,
                  parallel=T)
