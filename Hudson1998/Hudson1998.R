## Unit 1 - Applied Ecology of Infectious Diseases
## ZOO4926/ZOO6927 Spring 2014, Department of Biology
## University of Florida, Gainesville, Florida, USA
##
## Juliet R.C. Pulliam, 2016
##
## This file contains the model presented in:
## Hudson, PJ, AP Dobson, and D Newborn (1998) Prevention
## of population cycles by parasite removal. Science 
## 282: 2256-2258.
##

## Model 1: Equations 1 - 6
HudsonModel <- function(beta=0.11, # transmission coefficient
                       a = 1.8, # per capita birth rate (grouse)
                       b = 1.05, # per capita death rate (background,grouse)
                       c = 52, # 
                       p = 0, # proportion treated
                       k = 1.0, # aggregation of parasites in hosts (variance / mean)
                       alpha = 3*10^-4, # disease induced mortality rate
                       delta = 5*10^-4, # per capita reduction in birth rate for infected grouse
                       gamma = 7, # per capita mortality rate of free-living stages
                       mu = 1, # per capita mortality rate of adult worms
                       theta = 1, # maturation rate of chicks into reproductive class - not specified in paper
                       Delta = 4*10^-3, # density-dependent reduction in grouse reproduction and survival
                       plot=T, # whether or not to produce plots
                       logTraj=T, # whether to plot phase space trajectory on log scale
                       silent=T, # if F, returns model results as data frame
                       pop.init = 10^3, # initial population - not specified in paper
                       inf.init = .1, # initial infected - not specified in paper
                       years = 20,timestep=1/52){ # timeframe for output
  
  # Anderson and May 1980 - model 1
  model <- function(t,y,parms){
    with(c(as.list(y),parms),{
      dXdt <- aa*(X+Y)-bb*X-BB*X*Y+CC*Y
      dYdt <- BB*X*Y - (AA+bb+CC)*Y
      dNdt <- (aa-bb)*N - AA*Y
      list(c(dXdt,dYdt,dNdt))
    })
  }
  # Weekly rates
  parameters <- c(BB=beta,
                  aa = a,
                  bb = b,
                  AA = alpha,
                  CC = gamma,
                  N0 = pop.init,
                  I0 = inf.init)
  
  ts <- data.frame(lsoda(
    y = c(X = pop.init-inf.init,
          Y = inf.init,
          N = pop.init),# Initial conditions for population
    times = seq(0,years*52,timestep),   # Timepoints for evaluation
    func = model,                   		# Function to evaluate
    parms = parameters                	# Vector of parameters
  ))
  
  if(plot==T){
    par(mfcol=c(2,2),bty='L')
    plot(ts$time*timestep,log10(ts$N),xlab='Time (years)',ylab=expression(paste(log[10],' N')),main='Total population',type='l')
    plot(ts$time*timestep,log10(ts$Y),xlab='Time (years)',ylab=expression(paste(log[10],' Y')),main='Infected population',type='l')
    plot(ts$time*timestep,ts$Y/ts$N,xlab='Time (years)',ylab='Y/N',main='Prevalence',ylim=c(0,1),type='l')
    if(logTraj){
      plot(log10(ts$X),log10(ts$Y),xlab='Susceptible populaiton',ylab='Infected Population',main='Phase space trajectory',type='l')
      mtext('(log-log)',3,cex=.8,line=.5)
    }else{
      plot(ts$X,ts$Y,xlab='Susceptible populaiton',ylab='Infected Population',main='Phase space trajectory',type='l')
    }
  }
  if(silent==F){
    return(ts)
  }
}

## Model 2 - Equations 3, 5, 6
FreeLivingStageModel <- function(nu=10^-10, # transmission coefficient
                                 a = 4.3, # per capita birth rate
                                 b = 3.3, # per capita death rate (background)
                                 alpha = 14, # disease-induced mortality rate
                                 gamma = 0, # recovery rate
                                 mu = 3, # per capita death rate of free-living stage
                                 lambda = 10^6, # per capita shedding rate of infected hosts
                                 plot=T, # whether or not to produce plots
                                 logTraj=T, # whether to plot phase space trajectory on log scale
                                 silent=T, # if F, returns model results as data frame
                                 pop.init = 10^3, # inital population 
                                 inf.init = 1, # initial infected popuation
                                 free.init = 100, # intial populaiton of free-living stage
                                 years = 32,timestep=1/52){ # timeframe for output
  
  # Anderson and May 1980 - model 2
  model <- function(t,y,parms){
    with(c(as.list(y),parms),{
      dNdt <- (aa-bb)*N - AA*Y
      dYdt <- BB*W*(N-Y) - (AA+bb+CC)*Y
      dWdt <- LL*Y-(MM+BB*N)*W
      list(c(dNdt,dYdt,dWdt))
    })
  }
  
  # Annual rates
  parameters <- c(BB=nu,
                  aa = a,
                  bb = b,
                  AA = alpha,
                  CC = gamma,
                  MM = mu,
                  LL = lambda)	
  
  ts <- data.frame(lsoda(
    y = c(N = pop.init,
          Y = inf.init,
          W = free.init),# Initial conditions for population
    times = seq(0,years,timestep),   # Timepoints for evaluation
    func = model,                   		# Function to evaluate
    parms = parameters                	# Vector of parameters
  ))
  
  if(plot==T){
    par(mfcol=c(2,2),bty='L')
    plot(ts$time,log10(ts$N),xlab='Time (years)',ylab=expression(paste(log[10],' N')),main='Total population',type='l')
    plot(ts$time,log10(ts$Y),xlab='Time (years)',ylab=expression(paste(log[10],' Y')),main='Infected population',type='l')
    plot(ts$time,ts$Y/ts$N,xlab='Time (years)',ylab='Y/N',main='Prevalence',ylim=c(0,1),type='l')
    if(logTraj){
      plot(log10(ts$N-ts$Y),log10(ts$Y),xlab='Susceptible populaiton',ylab='Infected Population',main='Phase space trajectory',type='l')
      mtext('(log-log)',3,cex=.8,line=.5)
    }else{
      plot(ts$N-ts$Y,ts$Y,xlab='Susceptible populaiton',ylab='Infected Population',main='Phase space trajectory',type='l')
    }
  }
  if(silent==F){
    return(ts)
  }
}

BasicModel()
FreeLivingStageModel()