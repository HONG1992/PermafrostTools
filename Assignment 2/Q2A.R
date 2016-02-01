## Build surface temperature for 1 year
yearly.T <- function(mean, amplitude) {
  jul.day <- c(1:365)
  T.year <- mean+amplitude*sin(jul.day*2*pi/365)
  return(data.frame(jul.day,T.year))
}

years.to.run <- 10  # run for 10 years
k.factor <- 2.5
#follow the example from the help file, with a few exceptions

#Set parameters
dt <- toSeconds(1,"d") # 1 day time steps
ns <- 5  # save every 5 days
st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
k.thawed <- 1
k.frozen <- k.thawed*k.factor



dzmin <- .1 # 10 cm grid spacing 
zmax <- 5 # lowermost node
base <- 1.1 #  grid coarsening with depth
Ti <- 0


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

bc.upper <- c(Ti,rep(yearly.T(Ti,6)$T.year,years.to.run))
#bc.upper <- c(Ti,rep(Ti+step,nt-1))
#bc.upper <- seq(from = -5, to = 0.1, length.out = nt) #linear warming
bc.lower <- rep(0, nt)  #  Zero heat flux at depth

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # width of node [m]
                   Tj = rep(Ti,nz), # initial temperature [C]
                   soi = rep(1, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(k.frozen, nz), # thermal conductivity of soil matrix [W m-1 K-1]
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
  mat$k.soi[mat$Tj<0] <- k.frozen
  mat$k.soi[mat$Tj>=0] <- k.thawed
  
  #output only in specified interval
  check <- i/ns - floor(i/ns)
  if (check == 0) {
    out.Tj  <- rbind( out.Tj,t(mat$Tj))
    out.liq <- rbind(out.liq,t(mat$liq))
    out.t   <- t[i]
    print(paste(i,nt))
  }
}

for (r in c(1:740)) {
  plot(out.Tj[r,],mat$zj,type="l",xlim=c(-10,10))
  Sys.sleep(.1)
}

#Examine convergence of annual temperature plots
# 1, 74, 147 (+n*73) will be annual
dev.new()
plot(1,1,type="n",xlim=c(-10,10),ylim=c(-5,0),
     xlab="End-of-winter Temperature",ylab="Depth (m)",
     main="Figure 2-1")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
lines(out.Tj[1,],mat$zj,type="l")
lines(out.Tj[74,],mat$zj,type="l",col="blue")
lines(out.Tj[146,],mat$zj,type="l",col="green")
lines(out.Tj[219,],mat$zj,type="l",col="red")

text(1,-3,"Year 1")
text(-1.5,-2,"Year 2")
text(-4,-1,"year 3, 4")

## Examine convergence of annual mean temperatures
dev.new()
plot(1,1,type="n",xlim=c(-10,10),ylim=c(-5,0),
     xlab="Mean Annual Temperature",ylab="Depth (m)",
     main="Figure 2-2")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
lines(sapply(as.data.frame(out.Tj[2:74,]),mean),mat$zj,type="l")
       # mean for year 1
lines(sapply(as.data.frame(out.Tj[75:147,]),mean),mat$zj,col="green")  # mean for year 2
lines(sapply(as.data.frame(out.Tj[148:220,]),mean),mat$zj,col="blue") # mean for year 3
lines(sapply(as.data.frame(out.Tj[221:294,]),mean),mat$zj,col="red") # mean for year 4
text(1,-2,"Year 1")
text(-3,-2,"Year 2, 3, 4")


lines(out.Tj[74+73,],mat$zj,type="l")
lines(out.Tj[74+2*73,],mat$zj,type="l")
lines(out.Tj[74+3*73,],mat$zj,type="l")

#### Run experiment for more 
