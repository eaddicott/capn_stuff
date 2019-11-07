######################################################################################
# {capn} data: Kansas groundwater (gW.data) data
# Date: 12/30/2017 Updated: 11/06/2019
# Reference: Fenichel et al. (2016) PNAS, Addicott & Fenichel (2019) JEEM
#####################################################################################

# Load Required Packages -------------------------------------------------------------
if (!require("pacman")) install.packages("pacman") #pacman allows you to use the p_load function
p_load(capn, R.oo, repmis, ggplot2,devtools) #the p_load function checks if a library is installed. 
#                                            #If not, it installs it, then it attaches the called library

#capn documentation: https://cran.r-project.org/web/packages/capn/capn.pdf 
#repmis documentation: https://cran.r-project.org/web/packages/repmis/index.html

# Clear Workspace -------------------------------------------------------------------
rm(list=ls()) #clear workspace

# Get data for the problem set from Github ------------------------------------------
#KSwater_data <- readRDS(gzcon(url("https://github.com/eaddicott/capn_stuff/blob/master/KSwater_data.RDS?raw=true"))) # RDS file upload
source_data("https://github.com/eaddicott/capn_stuff/raw/master/KSwater_data.RData") #Rdata file upload
#str(KSwater_data) # this line will show you the structure of the data

## STRUCTURE INFO
#The object will be a list of 7 lists. Each of the 7 lists corresponds to a groundwater management district (1-5), 
# the outgroup(6), or the entire state of Kansas(7).
# Each list will have 11 named elements:
# [1] $gmdnum: int   1:7
# [2] $mlogitcoeff: df with dim = (#cropcodes x 23) containing the coefficients and intercept terms from the multinomial logit model
# [3] $mlogitmeans: df containing the mean for each of the variables run in the mlogit model
# [4] $cropamts: df containing the mean acres planted to each of the 5 crop cover types for each cropcode
# [5] $watercoeff: df containing the water withdrawal regression coefficients and intercept
# [6] $wwdmeans: df containing the means for each of the variables in the water withdrawal regression
# [7] $costcropacre: df dim=(1x5) containing the cost per acre of planting each of the 5 crops
# [8] $cropprices: df dim=(5x1) containing the per unit prices of each of the crops
# [9] $meanwater: num mean AF water in the region of interest
# [10]$recharge: num recharge rate for the region of interest
# [11]$watermax: num upper bound of domain for node space, max water observed in region.


# Modeling parameters to control   --------------------------------------------------
###
region <- 6 # Select region,  1:5 are GMD, 6 is outgroup, 7 is state
###
# After setting the region, create capn data structure
region_data <- KSwater_data[[region]] #double-brackets here important. Load in region specific data
gw.data <- datasetup(region) # if running line by line, be sure to run this function at the bottom of this script
#Economic parameters
dr <- 0.03 #discount rate

#System parameters
recharge <- region_data[[10]] #units are inches per year constant rate 
#recharge <- 1.25 uncomment to input your own for sensitivity analysis

#capN parameters
order <- 10 # approximaton order
NumNodes <- 100 #number of nodes
wmax <- region_data[[11]]
## More parameters for GMD level analysis in RDS file



# SYSTEM MODELS -------------------------------------------------------------
#source a separate script from github that contains the functions for gw system
if (!exists("cropFwater", mode = "function")) { # check if one of the groundwater functions exists
  require(devtools) # if not, load up package devtools so we can source an R file directly from github
  source_url("https://github.com/eaddicott/capn_stuff/raw/master/system_fns.R") 
  # Note: You may need to provide y/n input on creating a temporary directory.
  #       SHA-1 hash file info that is read out is extraneous
}

# Prepare capn input ----------------------------------------------------------

#     TEST capn               -----------------------------
######################################################################

#Prepare {capn}

#prepare capN-------------------------
Aspace <- aproxdef(order,0,wmax,dr) #defines the approximation space
nodes <- chebnodegen(NumNodes,0,wmax) #define the nodes

#prepare for simulation-------------------
simuData <- matrix(0,nrow = NumNodes, ncol = 5)

#simulate at nodes------------------
for(j in 1:NumNodes){
  simuData[j,1]<-nodes[j] #water depth nodes
  simuData[j,2]<-sdot(nodes[j],recharge,gw.data) #change in stock over change in time
  simuData[j,3]<- 0-WwdDs1(nodes[j],gw.data) # d(sdot)/ds, of the rate of change in the change of stock
  simuData[j,4]<-ProfDs1(nodes[j],gw.data) #Change in profit with the change in stock
  simuData[j,5]<-profit(nodes[j],gw.data) #profit
}


#recover approximating coefficents---------------------------
pC<-paprox(Aspace, simuData[,1],simuData[,2],simuData[,3],simuData[,4])  #the approximated coefficent vector for prices

#project shadow prices, value function, and inclusive wealth----------------
waterSim <- psim(pcoeff = pC,
                 stock = simuData[,1],
                 wval = simuData[,5],
                 sdot = simuData[,2])
#convert to data frame ---------------------------
waterSim <- as.data.frame(waterSim)

#use ggplot to plot the water shadow price function -------------------------  
ggplot() + 
  geom_line(data = waterSim, aes(x = stock, y = shadowp),
            color = 'blue') +
            labs( 
            x= "Stored groundwater",
            y = "Shadow price")  +
            theme(  #http://ggplot2.tidyverse.org/reference/theme.html
              axis.line = element_line(color = "black"), 
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA)
            )
### REVISE HERE AT END ------------------------------------------------------------
cat("if everything runs well the next line should say 17.44581", "\n")  # this is only for the state runs, NEED TO REVISE
cat("At 21.5 acre feet of water, the shadow price is" , psim(pC,21.5)$shadowp, "\n")

# PLOT VALUE --------------------------------------------------------------
#use ggplot plot the value function
ggplot() + 
  geom_line(data = waterSim[5:100,], aes(x = stock[5:100], y = vfun[5:100]),
            # four rows are removed because they are two far outside the data
            color = 'blue') +
  labs( 
    x= "Stored groundwater",
    y = "Intertemporial Welfare")  +
  theme(  #http://ggplot2.tidyverse.org/reference/theme.html
    axis.line = element_line(color = "black"), 
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

cat("you may get an warning that some number or rows were removed from the plot. This is not important.", "\n")

testme<-psim(pcoeff = pC, 
             stock = c(18.5,21.5), 
             wval = c(profit(18.5,gw.data),profit(21.5, gw.data)),
             sdot = c(sdot(18.5, recharge, gw.data),sdot(21.5, recharge, gw.data))
             )
testme

rm(j, testme)

# Convert regression output and summary statistics into capn parameters -------------------------------------
datasetup <- function(gmdnum) {
  # Region Specific Parameters ---------------------------------------------------
  #CASE 1: GMD 1-5
  if (gmdnum < 6) {
  # Crop Choice Parameters Setup ------------------------------------------------
  `mlogit.gmd.coeff` <- region_data[[2]]
  `mlogit.gmd.means` <- region_data[[3]]
  
  mlogit.gmd.coeff   <- data.frame(mlogit.gmd.coeff[-c(1,24,25),], 
                                   stringsAsFactors = FALSE)
  mlogit.gmd.means   <- data.frame(mlogit.gmd.means[-c(22,23), gmdnum + 1], 
                                   stringsAsFactors = FALSE)
  
  # Transpose betas, chop off variable names to create crop.coeff column
  crop.gmd.betas     <- t(mlogit.gmd.coeff[1,-c(1,ncol(mlogit.gmd.coeff))])
  crop.gmd.alphas    <- rep(NA, ncol(mlogit.gmd.coeff) - 2)
  crop.mean.vector <- append(as.numeric(unlist(mlogit.gmd.means[-c(1),])), 1)
  max.j <- ncol(mlogit.gmd.coeff) - 1
  for (j in 2:max.j) {
    crop.gmd.alphas[j - 1] <- as.numeric(mlogit.gmd.coeff[-c(1), j]) %*% 
      as.numeric(crop.mean.vector)
  }
  crop.coeff.gmd <- data.frame(crop.gmd.alphas,crop.gmd.betas, 
                               stringsAsFactors = FALSE)
  colnames(crop.coeff.gmd) <- c("alpha","beta")
  
  # Crop Amounts Parameter Setup ------------------------------------------------
  `crop.amts.gmd` <- region_data[[4]]
  crop.amts.gmd <- data.matrix(crop.amts.gmd[-c(1,8),-c(1)])
  crop.amts.gmd <- t(crop.amts.gmd)
  crop.amts.gmd <- transform(crop.amts.gmd, numeric)
  #Check that crop amounts are one more than 
  if (nrow(crop.amts.gmd) != max.j) {
    stop("BAD LENGTH")
  }
  # Water Withdrawal Parameter Setup --------------------------------------------
  `water.coeff` <- region_data[[5]]
  rmse <- as.numeric(water.coeff[36,2])
  water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
  
  # Beta ------------------------------------------------------------------------
  beta <- as.numeric(water.coeff[21,2])
  
  # Gamma -----------------------------------------------------------------------
  gamma <- as.numeric(water.coeff[c(1:10),2])
  gamma1 <- gamma[c(1,3,5,7,9)]
  gamma2 <- gamma[c(2,4,6,8,10)]
  
  # Alpha -----------------------------------------------------------------------
  alpha.water.coeff <- water.coeff[-c(1:11,21),2]
  `wwd.means` <- region_data[[6]]
  wwd.means <- data.frame(wwd.means[-c(1,11,22,23),gmdnum + 2], stringsAsFactors = FALSE)
  wwd.mean.vector <- append(as.numeric(unlist(wwd.means)),1)
  alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
  alpha <- alpha + (rmse**2)/2 #correction for log estimation
  alpha <- as.vector(alpha, mode = "numeric")
  }
  #CASE 2: GMD 6 Outgroup
  else if (gmdnum == 6){
    # Crop Choice Parameters Setup ------------------------------------------------
    `mlogit.state.coeff` <- KSwater_data[[7]][[2]]
    `mlogit.state.means` <- KSwater_data[[7]][[3]]
    mlogit.state.coeff   <- data.frame(mlogit.state.coeff[-c(1,24,25),], 
                                       stringsAsFactors = FALSE)
    mlogit.state.means   <- data.frame(mlogit.state.means[-c(22,23),], 
                                       stringsAsFactors = FALSE)
    # Transpose betas, chop off variable names to create crop.coeff column
    crop.state.betas     <- t(mlogit.state.coeff[1,-c(1,ncol(mlogit.state.coeff))])
    crop.state.alphas    <- rep(NA, ncol(mlogit.state.coeff) - 2)
    crop.mean.vector <- append(as.numeric(mlogit.state.means$mean[-c(1)]),1)
    max.j <- ncol(mlogit.state.coeff) - 1
    for (j in 2:max.j) {
      crop.state.alphas[j - 1] <- as.numeric(mlogit.state.coeff[-c(1), j]) %*% 
        as.numeric(crop.mean.vector)
    }
    crop.coeff.state <- data.frame(crop.state.alphas,crop.state.betas, 
                                   stringsAsFactors = FALSE)
    colnames(crop.coeff.state) <- c("alpha","beta")
    
    # Crop Amounts Parameter Setup ------------------------------------------------
    `crop.amts.state` <- KSwater_data[[7]][[4]]
    crop.amts.state <- data.matrix(crop.amts.state[-c(1,8),-c(1)])
    crop.amts.state <- t(crop.amts.state)
    
    #Check that crop amounts are one more than 
    if (nrow(crop.amts.state) != max.j) {
      stop("BAD LENGTH")
    }
    
    # Water Withdrawal Parameter Setup --------------------------------------------
    `water.coeff` <- region_data[[5]]
    rmse <- as.numeric(water.coeff[36,2])
    water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
    
    # Beta ------------------------------------------------------------------------
    beta <- as.numeric(water.coeff[21,2])
    
    # Gamma -----------------------------------------------------------------------
    gamma <- as.numeric(water.coeff[c(1:10),2])
    gamma1 <- gamma[c(1,3,5,7,9)]
    gamma2 <- gamma[c(2,4,6,8,10)]
    
    # Alpha -----------------------------------------------------------------------
    alpha.water.coeff <- water.coeff[-c(1:11,21),2]
    `wwd.means` <- region_data[[6]]
    wwd.means <- data.frame(wwd.means[-c(1,11,22,23),gmdnum + 2], stringsAsFactors = FALSE)
    wwd.mean.vector <- append(as.numeric(unlist(wwd.means)),1)
    alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
    alpha <- alpha + (rmse**2)/2 #correction for log estimation
    alpha <- as.vector(alpha, mode = "numeric")
  }
  #CASE 3: STATE
  else if (gmdnum == 7){
    # Crop Choice Parameters Setup ------------------------------------------------
    `mlogit.state.coeff` <- region_data[[2]]
    `mlogit.state.means` <- region_data[[3]]
    mlogit.state.coeff   <- data.frame(mlogit.state.coeff[-c(1,24,25),], 
                                       stringsAsFactors = FALSE)
    mlogit.state.means   <- data.frame(mlogit.state.means[-c(22,23),], 
                                       stringsAsFactors = FALSE)
    # Transpose betas, chop off variable names to create crop.coeff column
    crop.state.betas     <- t(mlogit.state.coeff[1,-c(1,ncol(mlogit.state.coeff))])
    crop.state.alphas    <- rep(NA, ncol(mlogit.state.coeff) - 2)
    crop.mean.vector <- append(as.numeric(mlogit.state.means$mean[-c(1)]),1)
    max.j <- ncol(mlogit.state.coeff) - 1
    for (j in 2:max.j) {
      crop.state.alphas[j - 1] <- as.numeric(mlogit.state.coeff[-c(1), j]) %*% 
        as.numeric(crop.mean.vector)
    }
    crop.coeff.state <- data.frame(crop.state.alphas,crop.state.betas, 
                                   stringsAsFactors = FALSE)
    colnames(crop.coeff.state) <- c("alpha","beta")
    
    # Crop Amounts Parameter Setup ------------------------------------------------
    `crop.amts.state` <- region_data[[4]]
    crop.amts.state <- data.matrix(crop.amts.state[-c(1,8),-c(1)])
    crop.amts.state <- t(crop.amts.state)
    
    #Check that crop amounts are one more than 
    if (nrow(crop.amts.state) != max.j) {
      stop("BAD LENGTH")
    }
    
    # Water Withdrawal Parameter Setup --------------------------------------------
    `water.coeff` <- region_data[[5]]
    rmse <- as.numeric(water.coeff[36,2])
    water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
    
    # Beta ------------------------------------------------------------------------
    beta <- as.numeric(water.coeff[21,2])
    
    # Gamma -----------------------------------------------------------------------
    gamma <- as.numeric(water.coeff[c(1:10),2])
    gamma1 <- gamma[c(1,3,5,7,9)]
    gamma2 <- gamma[c(2,4,6,8,10)]
    
    # Alpha -----------------------------------------------------------------------
    alpha.water.coeff <- water.coeff[-c(1:11,21),2]
    `wwd.means` <- region_data[[6]]
    wwd.means <- data.frame(wwd.means[-c(1,11,22,23),], stringsAsFactors = FALSE)
    wwd.mean.vector <- append(as.numeric(wwd.means$all),1)
    alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
    alpha <- alpha + (rmse**2)/2 #correction for log estimation
    alpha <- as.vector(alpha, mode = "numeric")
  }
  
  # CONSISTENT PARAMETERS ------------------------------------------------------
  # Cost Crop Acre 
  `cost.crop.acre` <- region_data[[7]]
  cost.crop.acre <- as.vector(unlist(cost.crop.acre), mode = "numeric")
  
  # Crop Prices 
  `crop.prices` <- region_data[[8]]
  crop.prices <- as.vector(unlist(crop.prices), mode = "numeric")

  
  # Store Parameters in Data Structure ------------------------------------------
  gw.data <- rep(NA,9) #preallocate
  # GMD 1-5
  if (gmdnum < 6) {
    gw.data <- list(crop.coeff.gmd, crop.amts.gmd, alpha, beta, gamma, 
                    gamma1, gamma2, crop.prices, cost.crop.acre)
  }
  # GMD 6 (Outgroup) or State
  if (gmdnum >= 6) {
    gw.data <- list(crop.coeff.state, crop.amts.state, alpha, beta, gamma, 
                        gamma1, gamma2, crop.prices, cost.crop.acre)
  }

  return(gw.data)
} 
