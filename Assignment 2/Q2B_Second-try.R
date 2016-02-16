rm(list=ls()) # Clean up workspace

library(PermafrostTools)
source("RequiredFunctions.R")


##  Set Global Parameters

Ti <- -1
#make material properties and discretization
dzmin <- .01 # 10 cm grid spacing 
zmax <- 5 # lowermost node
base <- 1.1 #  grid coarsening with depth
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

soil.initial <- data.frame( zj = -soild$z, # depth of node [m]
                            dz = soild$dz, # width of node [m]
                            Tj = rep(Ti,nz), # initial temperature [C]
                            soi = rep(0.5, nz), # volumetric proportion of soil matrix [0-1]
                            wat = rep(0.5, nz), # volumetric proportion of total water [0-1]
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
    ns <- 5  # save every 5 days
    st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
    t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
    nt    <- length(t)
    
    #== Boundary Conditions ==============
    bc.upper <- yearly.T(Ti,6,years=years.to.run)$T.year # oscillate about mean Ti
    bc.lower <- rep(0.2, nt)  #  heat flux with depth
    
    #== Target arrays to hold results =========
    out.Tj  <- rbind(NULL,t(mat$Tj))
    out.liq <- rbind(NULL,t(mat$liq))
    out.t   <- 0
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
        out.t   <- c(out.t,t[i])
        out.mat <- mat
        print(paste(i,nt))
      }
    }
  return(list(mat=out.mat,out.Tj=out.Tj,bc=bc.upper,t=out.t))
  } # END Spinup

##  Function to warm up soil discretization, runs until no more permafrost or until stopyears is reached
Warm <- function(soil,years,warming){  #could add parameters variable
  ## Warming is annual warming
  mat <- soil$mat[,c("zj","dz","Tj","soi","wat","liq","ice","k.soi","c.soi")]
  years.to.run <- years  # run for number of years before stopping
  Ti <- -0.5
  #make material properties and discretization
  
  
  #Set time parameters
  dt <- toSeconds(1,"d") # 1 day time steps
  ns <- 5  # save every 5 days
  st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
  t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
  nt    <- length(t)
  
  #== Boundary Conditions ==============
  bc.upper <- yearly.T(Ti,6,years=years.to.run,warming=warming)$T.year # oscillate about mean Ti
  bc.lower <- rep(0.2, nt)  #   heat flux at depth
  
  #== Target arrays to hold results =========
  out.Tj  <- rbind(NULL,t(mat$Tj))
  out.liq <- rbind(NULL,t(mat$liq))
  out.t   <- 0
  out.dz  <- mat$dz
  out.mat  <- mat
  out.PF <- array(dim=c(nz,years.to.run),
                  dimnames = list(Depth=as.character(round(mat$zj,3)), 
                                  Time=paste("y",c(1:years.to.run),sep="")))  # empty array dim # years and num soil disc
  
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
      out.t   <- c(out.t,t[i])
      print(paste(i,nt))
    }
    
    #Test for permafrost,
    if (i%%365 == 0 & i>365*2) {
      snapshot_2yr <- apply(out.Tj[((i/5)-146):(i/5),],2,max) # get max temp for last 2 years (note 146 = 2yr*(365/5), where 5=ns)
      out.PF[,i/365]  <- snapshot_2yr < 0 # build true false if PF or no
    }
    
    if (sum(out.PF[,i/365]) == 0 & i>365*5) break  #stop running if no more PF, but run for at least 5yr
  }
  return(list(mat=out.mat,out.Tj=out.Tj,PF=out.PF,Bc=bc.upper,t=out.t))
} # End Warm

## Visualize!!  # Watch a soil profile evolve
vis <- function(PF.output,skip=1){
  for (i in seq(1,length(PF.output$out.Tj[,1]),skip)) {
    plot(PF.output$out.Tj[i,],PF.output$mat$zj, type="l", xlim=c(-10,3),
         xlab = "Temperature (C)", ylab = "Depth (m)")
    abline(v=0,col="red",lty=2)
    Sys.sleep(0.1)
  }
}

## Spin up for 20 years
B <- Spinup(soil.initial,20,MAAT=-0.5)

##  Watch spinup
#vis(B)

##  Warm Soil
W<-Warm(B,50,0.1)  #50 year warming @ 0.1 deg/yr

## Watch Warming
#vis(W,skip=15)

# Make Plots
dev.new()
plot(seq(1,10,length.out=365),B$bc[0:365*10],xlab="Years Elapsed", ylab="Temperature (C)",type="l",
     main="Spin-up temperatures (First 10 Years)")

dev.new()
plot(B$t/toSeconds(1,'y'),B$out.Tj[,41], xlab="Years Elapsed", ylab="Temperature (C)",type="l",
     main="Spin-up Temperature at 4.71 m")

dev.new()
plot(seq(1,10,length.out=365),W$Bc[0:365*10],xlab="Years Elapsed", ylab="Temperature (C)",type="l",
     main="Simulated Warming (First 10 Years)")

dev.new()
plot(W$t/toSeconds(1,'y'),W$out.Tj[,39], xlab="Years Elapsed", ylab="Temperature (C)",type="l",
     main="Ground Temperature at 3.83 m")

dev.new()
plot(apply(W$out.Tj[(1388-74*2):(1388-74),],2,mean),W$mat$zj, xlab="Temperature (C)", ylab="Depth (m)",type="l",
     main="Mean Ground Temperature for year 18")
abline(v=0,col="red",lty=2)

