rm(list=ls())  # Clean up workspace
  
library(PermafrostTools)
source("RequiredFunctions.R")
  
  
  ##  Set Global Parameters
  
  Ti <- -1.5
  #make material properties and discretization
  dzmin <- 0.01 # 10 cm grid spacing 
  zmax <- 150 # lowermost node
  base <- 1.3 #  grid coarsening with depth
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
  soil.initial$Tj=Ti-soil.initial$zj*(0.02/soil.initial$k.soi)
  
  
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
    dt <- toSeconds(73,"d") # 73 day time steps
    ns <- 3  # save every 3*5=15 days
    st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
    t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
    nt    <- length(t)
    
    #== Boundary Conditions ==============
    bc.upper <- yearly.T(Ti,0,years=years.to.run)$T.year # flat temp at Ti
    bc.lower <- rep(0.02, nt)  #  heat flux with depth
    
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
    return(list(mat=out.mat,out.Tj=out.Tj,bc=bc.upper))
  } # END Spinup
  
  ##  Function to warm up soil discretization, runs until no more permafrost or until stopyears is reached
  Warm <- function(soil,years,warming){  #could add parameters variable
    ## Warming is annual warming
    mat <- soil$mat[,c("zj","dz","Tj","soi","wat","liq","ice","k.soi","c.soi")]
    years.to.run <- years  # run for number of years before stopping
    Ti <- -1.5
    #make material properties and discretization
    
    
    #Set time parameters
    dt <- toSeconds(73,"d") # 5 day time steps
    ns <- 5  # save every year
    st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
    t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
    nt    <- length(t)
    
    #== Boundary Conditions ==============
    bc.upper <- yearly.T(Ti,0,years=years.to.run,warming=warming)$T.year # linear ramp from Ti
    bc.upper[bc.upper>10]=10
    bc.lower <- rep(0.02, nt)  #   heat flux at depth
    
    #== Target arrays to hold results =========
    out.Tj  <- rbind(NULL,t(mat$Tj))
    out.liq <- rbind(NULL,t(mat$liq))
    out.t   <- t[0]
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
        out.t   <- t[i]
        print(paste(i,nt))
      }
      
      #Test for permafrost,
      if (i%%5 == 0 & i>5*2) {
        snapshot <- out.Tj[(i/5),] # get max temp for last 2 years (note 146 = 2yr*(365/5), where 5=ns)
        out.PF[,(i/5)]  <- snapshot < 0 # build true false if PF or no
      }
      
      if (sum(out.PF[,i/5]) == 0 & i>5*5) break  #stop running if no more PF, but run for at least 5yr
    }
    return(list(mat=out.mat,out.Tj=out.Tj,PF=out.PF,Bc=bc.upper))
  } # END Warm
 
  vis <- function(PF.output,skip=1){
    for (i in seq(1,length(PF.output$out.Tj[,1]),skip)) {
      plot(PF.output$out.Tj[i,],PF.output$mat$zj, type="l", xlim=c(-10,3),
           xlab = "Temperature (C)", ylab = "Depth (m)")
      abline(v=0,col="red",lty=2)
      Sys.sleep(0.1)
    }
  }
  
  ## Spin up for 500 years
  B500 <- Spinup(soil.initial,500,MAAT=Ti)
  
  ##  Warm Soil
  W500<-Warm(B500,5000,0.1)  #5000 year warming @ 0.1 deg/yr


## Watch Warming
dev.new()
par(mfrow=c(1,1))
vis(W500,skip=50)


# Make Plots

dev.new() 
par(mfrow=c(1,2))
## Initial soil profile
plot(soil.initial$Tj,soil.initial$zj,type="l",
     xlab="Temperature (C)",
     ylab="Depth (m)",xlim=c(-1.5,0.2),
     main="Initial Soil Temperatures")
abline(v=0,col="red",lty=2)

## Spin-up soil profile
plot(B500$out.Tj[200,],B500$mat$zj,type="l",
     xlab="Temperature (C)",
     ylab="Depth (m)",xlim=c(-1.5,0.2),
     main="Soil Temperatures after 500yr Spin-up")
abline(v=0,col="red",lty=2)


dev.new()
PFloss<-sapply(c(1:length(W500$mat$zj)),function(x) min(which(W500$PF[x,]==F))) # which year did PF thaw completely?
plot(PFloss,W500$mat$zj,type="l",xlab="Years",ylab="Depth (m)",lty=2,
     main="Location of 0 Degree Isotherm")

dev.new()
plot(seq(1,max(PFloss),length.out=dim(W500$out.Tj)[1]),W500$out.Tj[,36], xlab="Years Elapsed", ylab="Temperature (C)",type="l",
     main="Ground Temperature at 110.9 m")
abline(v=2925,lty=2,col="grey")
abline(h=0,lty=3,col="red")

dev.new()
par(mfrow=c(1,2))
plot(W500$out.Tj[2850,],W500$mat$zj, xlab="Temperature (C)", ylab="Depth (m)",type="l",
     main="Year 2850")
abline(h=-110.9,lty=2,col="grey")
abline(v=0,lty=3,col="red")

plot(W500$out.Tj[3100,],W500$mat$zj, xlab="Temperature (C)", ylab="Depth (m)",type="l",
     main="Year 3100")
abline(h=-110.9,lty=2,col="grey")
abline(v=0,lty=3,col="red")

##################
# dev.new()
# PFloss<-sapply(c(1:length(W500$mat$zj)),function(x) min(which(W500$PF[x,]==F))) # which year did PF thaw completely?
# plot(PFloss,W500$mat$zj,type="l",xlab="Years",ylab="Depth (m)",lty=2,
#      main="Location of 0 Degree Isotherm")
# 
# dev.new()
# plot(seq(1,max(PFloss),length.out=dim(W500$out.Tj)[1]),W500$out.Tj[,40], xlab="Years Elapsed", ylab="Temperature (C)",type="l",
#      main="Ground Temperature at 117.4 m \n (2nd Discretization)")
# text(3000,-.3,"dzmin = 1.0 \n zmax = 150 \n base = 1.05"
#          )
# 
# 
# dev.new()
# par(mfrow=c(1,2))
# plot(W500$out.Tj[2850,],W500$mat$zj, xlab="Temperature (C)", ylab="Depth (m)",type="l",
#      main="Year 2850")
# abline(h=-110.9,lty=2,col="grey")
# abline(v=0,lty=3,col="red")
# 
# plot(W500$out.Tj[3100,],W500$mat$zj, xlab="Temperature (C)", ylab="Depth (m)",type="l",
#      main="Year 3100")
# abline(h=-110.9,lty=2,col="grey")
# abline(v=0,lty=3,col="red")
# 
# 
# 


