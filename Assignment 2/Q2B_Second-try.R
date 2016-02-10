## Set parameters
## Function to run for 1 year and check max(T) at each depth
## 
rm(list=ls())

setwd("C:\\Users\\Nick\\Code\\PermafrostTools\\Assignment 2")
source("RequiredFunctions.R")


##  Set Global Parameters

Ti <- 0
dzmin <- .1 # 10 cm grid spacing 
zmax <- 5 # lowermost node
base <- 1.1 #  grid coarsening with depth
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

soil.initial <- data.frame( zj = -soild$z, # depth of node [m]
                dz = soild$dz, # width of node [m]
                Tj = rep(Ti,nz), # initial temperature [C]
                soi = rep(1, nz), # volumetric proportion of soil matrix [0-1]
                wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
            
)

type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.0)  # parameter(s) for unfrozen water content parameterization (van genuchten)
# [1]: van Genuchten alpha [mm-1], 
#'                     unfrozen.par[2]: van Genuchten n [-], unfrozen.par[3]: 
#'                     residual water content [m3/m3].
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

## Function to input mat and output mat + years
Spinup <- function(soil,years){  #could add parameters variable
  mat <- soil.initial
  years.to.run <- 2  # run for number of years 
  #make material properties and discretization
  
    
    #Set time parameters
    dt <- toSeconds(1,"d") # 1 day time steps
    ns <- 365  # save every year 
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
    out.dz  <- mat$dz
   
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
    }
  return(list(mat=mat,out.Tj=out.Tj))
  } # END Spinup

##  Function to input mat, runs until no more permafrost or until stopyears is reached
##  After every year
Run <- function(soil, stopyears) {
  print("in progress")
}
  
  
