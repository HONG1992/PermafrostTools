rm(list=ls())
library(PermafrostTools)
source("RequiredFunctions.R")

### Comparison with analytical solution

######################################
## First, build analytical solution ##
######################################
Ottawa <- PF_site(heatflux = 0.2, heatcap=1.8e6, th_cond=2.5, MST=-2)
OtClim <- Climate(Amplitude=4, MAGST = -4)

Ott<-Stepchange(Ottawa,OtClim,+3)

f=c(1,2,4,8,15)
lapply(f,function(x) Ott$add.step(x,"y"))
dat1<-Ott$get.values()


########################################
## Plot analytical solution as a base ##
########################################

dev.new()

plot(dat1$Ti,dat1$depths,type="l",xlim=c(-4,2),axes=F,ann=F)
axis(3);mtext("Temperature (C)", side=3, line=2)
axis(2);mtext("Depth (m)", side=2, line=3)
box()
lines(dat1[,3], dat1$depths); text(dat1[10,3], dat1$depths[8], "1y")
lines(dat1[,4], dat1$depths); text(dat1[10,4], dat1$depths[8], "2y")
lines(dat1[,5], dat1$depths); text(dat1[10,5], dat1$depths[8], "4y")
lines(dat1[,6], dat1$depths); text(dat1[10,6]-.05, dat1$depths[8], "8y")
lines(dat1[,7], dat1$depths); text(dat1[10,7]+.05, dat1$depths[8], "15y")
lines(x=rep(0,111),y=seq(-100,10,1),col="red",lty=4)

legend(0,-20,legend=c("Analytical", "Numeric"),col=c("black","blue"),lty=c(1,3),
       lwd=c(1,2.5),title ="Solution Type")

##############################################
##  Next we'll build the numerical solution ##
##############################################

#follow the example from the help file, with a few exceptions

#Set parameters
dt <- toSeconds(1,"d") # 1 day time steps
ns <- 365  # save every year
st <- toSeconds(15,"y")  # [s] run it for 15 year

dzmin <- .5 # 10 cm grid spacing 
zmax <- 75 # lowermost node
base <- 1.25 #  grid coarsening with depth
Ti <--4
step <- +3

t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations


type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization (van genuchten)
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

#Boundary Conditions

bc.upper <- c(Ti,rep(Ti+step,nt-1))
#bc.upper <- seq(from = -5, to = 0, length.out = nt) #linear warming
bc.lower <- rep(0.2, nt)

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # width of node [m]
                   Tj = Ti+soild$z*(Ottawa$HeatFlux/Ottawa$ThermalConduct), # initial temperature [C]
                   soi = rep(1, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

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

#################################
#  Annotate with numeric data  ##
#################################

#mtext("Figure 1", side=1, line=1,cex=2)

text(-3,-40,paste(
"Time Step = 1 day \n
     Grid spacing = ",dzmin,"m \n
Grid Coarsening = ",base,sep=""))

lines(out.Tj[1,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[2,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[3,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[5,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[9,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[16,],mat$zj,lty=3,lwd=2.5, col="blue")


##############################################
##  Make new plotting window for new figure ##
##############################################
dev.new()
plot(dat1$Ti,dat1$depths,type="l",xlim=c(-4,2),axes=F,ann=F)
axis(3);mtext("Temperature (C)", side=3, line=2)
axis(2);mtext("Depth (m)", side=2, line=3)
box()
lines(dat1[,3], dat1$depths); text(dat1[10,3], dat1$depths[8], "1y")
lines(dat1[,4], dat1$depths); text(dat1[10,4], dat1$depths[8], "2y")
lines(dat1[,5], dat1$depths); text(dat1[10,5], dat1$depths[8], "4y")
lines(dat1[,6], dat1$depths); text(dat1[10,6]-.05, dat1$depths[8], "8y")
lines(dat1[,7], dat1$depths); text(dat1[10,7]+.05, dat1$depths[8], "15y")
lines(x=rep(0,111),y=seq(-100,10,1),col="red",lty=4)

legend(0,-20,legend=c("Analytical", "Numeric"),col=c("black","blue"),lty=c(1,3),
       lwd=c(1,2.5),title ="Solution Type")


# Run numerical model with coarse step size -------------------------------

#follow the example from the help file, with a few exceptions

#Set parameters
dt <- toSeconds(73,"d") # 73 day time steps
ns <- 5  # save every year
st <- toSeconds(15,"y")  # [s] run it for 15 year

dzmin <- .5 # 10 cm grid spacing 
zmax <- 75 # lowermost node
base <- 1.25 #  grid coarsening with depth
Ti <--4
step <- +3

t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations


type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization (van genuchten)
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

#Boundary Conditions

bc.upper <- c(Ti,rep(Ti+step,nt-1))
#bc.upper <- seq(from = -5, to = 0, length.out = nt) #linear warming
bc.lower <- rep(0.2, nt)

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # width of node [m]
                   Tj = Ti+soild$z*(Ottawa$HeatFlux/Ottawa$ThermalConduct), # initial temperature [C]
                   soi = rep(1, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

#Target arrays to hold results =========
out.Tj  <- rbind(NULL,t(mat$Tj))
out.liq <- rbind(NULL,t(mat$liq))
out.t   <- t[0]
out.dz  <- mat$dz

#Loop over time and solve =======
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

# Annotate Figure 2 -------------------------------------------------------

#mtext("Figure 2", side=1, line=1,cex=2)

text(-3,-40,paste(
  "Time Step = 73 day \n
     Grid spacing = ",dzmin,"m \n
Grid Coarsening = ",base,sep=""))

lines(out.Tj[1,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[2,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[3,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[5,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[9,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[16,],mat$zj,lty=3,lwd=2.5, col="blue")

##############################################
##  Make new plotting window for new figure 3##
##############################################
dev.new()
plot(dat1$Ti,dat1$depths,type="l",xlim=c(-4,2),axes=F,ann=F)
axis(3);mtext("Temperature (C)", side=3, line=2)
axis(2);mtext("Depth (m)", side=2, line=3)
box()
lines(dat1[,3], dat1$depths); text(dat1[10,3], dat1$depths[8], "1y")
lines(dat1[,4], dat1$depths); text(dat1[10,4], dat1$depths[8], "2y")
lines(dat1[,5], dat1$depths); text(dat1[10,5], dat1$depths[8], "4y")
lines(dat1[,6], dat1$depths); text(dat1[10,6]-.05, dat1$depths[8], "8y")
lines(dat1[,7], dat1$depths); text(dat1[10,7]+.05, dat1$depths[8], "15y")
lines(x=rep(0,111),y=seq(-100,10,1),col="red",lty=4)

legend(0,-20,legend=c("Analytical", "Numeric"),col=c("black","blue"),lty=c(1,3),
       lwd=c(1,2.5),title ="Solution Type")


# Run numerical model with coarse soil profile -------------------------------

#Set parameters
dt <- toSeconds(1,"d") # 1 day time steps
ns <- 365  # save every year
st <- toSeconds(15,"y")  # [s] run it for 15 year

dzmin <- .5 # 50 cm grid spacing 
zmax <- 75 # lowermost node
base <- 2#  grid coarsening with depth
Ti <--4
step <- +3

t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations


type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization (van genuchten)
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

#Boundary Conditions

bc.upper <- c(Ti,rep(Ti+step,nt-1))
#bc.upper <- seq(from = -5, to = 0, length.out = nt) #linear warming
bc.lower <- rep(0.2, nt)

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # width of node [m]
                   Tj = Ti+soild$z*(Ottawa$HeatFlux/Ottawa$ThermalConduct), # initial temperature [C]
                   soi = rep(1, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

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

# Annotate Figure 3 -------------------------------------------------------

#mtext("Figure 3", side=1, line=1,cex=2)

text(-3,-40,paste(
  "Time Step = 1 day \n
     Grid spacing = ",dzmin,"m \n
Grid Coarsening = ",base,sep=""))

lines(out.Tj[1,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[2,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[3,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[5,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[9,],mat$zj,lty=3,lwd=2.5, col="blue")
lines(out.Tj[16,],mat$zj,lty=3,lwd=2.5, col="blue")