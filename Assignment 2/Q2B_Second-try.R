## Set parameters
## Function to run for 1 year and check max(T) at each depth
## 
rm(list=ls())

setwd("C:\\Users\\Nick\\Git\\PermafrostTools\\Assignment 2")
source("RequiredFunctions.R")


##  Set Global Parameters

Ti <- -2
#make material properties and discretization
dzmin <- .01 # 10 cm grid spacing 
zmax <- 5 # lowermost node
base <- 1.1 #  grid coarsening with depth
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

soil.initial <- data.frame( zj = -soild$z, # depth of node [m]
                            dz = soild$dz, # width of node [m]
                            Tj = rep(Ti,nz), # initial temperature [C]
                            soi = rep(0.1, nz), # volumetric proportion of soil matrix [0-1]
                            wat = rep(0.9, nz), # volumetric proportion of total water [0-1]
                            liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                            ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                            k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                            c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
                            
)

type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization (van genuchten)
# [1]: van Genuchten alpha [mm-1], 
#'                     unfrozen.par[2]: van Genuchten n [-], unfrozen.par[3]: 
#'                     residual water content [m3/m3].
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

## Function to spin up experiment for a number of years subject to oscillating boundary conditions
Spinup <- function(soil,years,MAAT){  #could add parameters variable
  Ti <- MAAT
  mat <- soil
  years.to.run <- years  # run for number of years 
  
  
    
    #Set time parameters
    dt <- toSeconds(1,"d") # 1 day time steps
    ns <- 15  # save every 15 days
    st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
    t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
    nt    <- length(t)
    
    #== Boundary Conditions ==============
    bc.upper <- rep(yearly.T(Ti,6)$T.year,years.to.run) # oscillate about mean Ti
    bc.lower <- rep(0, nt)  #  Zero heat flux at depth
    
    #== Target arrays to hold results =========
    out.Tj  <- rbind(NULL,t(mat$Tj))
    out.liq <- rbind(NULL,t(mat$liq))
    out.t   <- t[0]
    out.mat  <- mat
   
    #== Loop over time and solve =======
    for (i in 1:nt) {
      mat <- HcGroundImplicit(dt, mat, bc.lower[i], bc.upper[i], unfrozen.type = unfrozen.type,
                                  unfrozen.par = unfrozen.par, bcutype = bcutype, bcltype = bcltype,
                                  layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
      
      
      #output only in specified interval
      check <- i/ns - floor(i/ns)
      if (check == 0) {
        out.Tj  <- rbind( out.Tj,t(mat$Tj))
        out.liq <- rbind(out.liq,t(mat$liq))
        out.t   <- t[i]
        out.mat <- mat
        print(paste(i,nt))
      }
    }
  return(list(mat=out.mat,out.Tj=out.Tj))
  } # END Spinup

##  Function to warm up soil discretization, runs until no more permafrost or until stopyears is reached
Warm <- function(soil,years,warming){  #could add parameters variable
  ## Warming is annual warming
  mat <- soil$mat[,c("zj","dz","Tj","soi","wat","liq","ice","k.soi","c.soi")]
  years.to.run <- years  # run for number of years before stopping
  Ti <- soil$mat$Tj[1]
  #make material properties and discretization
  
  
  #Set time parameters
  dt <- toSeconds(1,"d") # 1 day time steps
  ns <- 5  # save every 5 days
  st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
  t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
  nt    <- length(t)
  
  #== Boundary Conditions ==============
  bc.upper <- yearly.T(Ti,6,years=years.to.run,warming=warming)$T.year # oscillate about mean Ti
  bc.lower <- rep(0, nt)  #  Zero heat flux at depth
  
  #== Target arrays to hold results =========
  out.Tj  <- rbind(NULL,t(mat$Tj))
  out.liq <- rbind(NULL,t(mat$liq))
  out.t   <- t[0]
  out.dz  <- mat$dz
  out.mat  <- mat
  out.PF <- array(dim=c(nz,years.to.run))  # empty array dim # years and num soil disc
  
  #== Loop over time and solve =======
  for (i in 1:nt) {
    mat <- HcGroundImplicit(dt, mat, bc.lower[i], bc.upper[i], unfrozen.type = unfrozen.type,
                            unfrozen.par = unfrozen.par, bcutype = bcutype, bcltype = bcltype,
                            layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
    
    
    #output only in specified interval
    check <- i/ns - floor(i/ns)
    if (check == 0) {
      out.Tj  <- rbind( out.Tj,t(mat$Tj))
      out.liq <- rbind(out.liq,t(mat$liq))
      out.t   <- t[i]
      print(paste(i,nt))
    }
    
    #Test for permafrost,
    if (i%%365 == 0 & i>365*2) {
      snapshot_2yr <- apply(out.Tj[((i/5)-146):(i/5),],2,max) # True if <0 for last 2 yr (note 146 = 2yr*(365/5), where 5=ns)
      out.PF[,i/365]  <- snapshot_2yr < 0 # build true false if PF or no
      
    }
  }
  return(list(mat=out.mat,out.Tj=out.Tj,PF=out.PF,Bc=bc.upper))
}

## Visualize!!  # change FPS so it only prints every nth row
vis <- function(PF.output,fps=0.1){
  for (i in seq(1,length(PF.output$out.Tj[,1]),1)) {
    plot(PF.output$out.Tj[i,],PF.output$mat$zj, type="l", xlim=c(-10,3))
    Sys.sleep(fps)
  }
}

## Spin up for 10 years
B <- Spinup(soil.initial,3,MAAT=-2)

##  Watch spinup
vis(B)

##  Warm Soil
W<-Warm(B,10,3)

## Watch Warming
vis(W,fps=.05)

#W$out.Tj[last year:now,]
#

