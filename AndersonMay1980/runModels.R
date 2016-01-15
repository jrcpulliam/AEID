rm(list=ls()) # Clears the workspace

## The following line needs to be run once, but after installing the package won't
## need to be run again. To run the line, remove the preceeding "#":
# install.packages("deSolve") # Installs the package used for numerical integration of ODE's

require(deSolve) # Loads the installed package

load("AndersonMay1980.Rdata") # Loads the model function definitions

BasicModel() # Runs basic model with default parameters
BasicModel(beta=0) # Runs the basic model, setting beta to 0
BasicModel(alpha=10,years=30) # etc...

FreeLivingStageModel() # Runs free-living stage model with default parameters
