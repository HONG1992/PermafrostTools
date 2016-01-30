library(PermafrostTools)
library(ggplot2)


time = c(1:365)
dep=seq(1,15,2.5)
b<-TempAnalyticHarmonic(-1, As=10,z=dep,t=time,3e8/pi,kappa=.5)

dim(b)

### Comparison with analytical solution

## First we'll build the numerical solution
#Set parameters
dt <- 3600  # 1 hour time steps
ns <- 240  # save every 10 days
st <- 365  # run it for 1 year
dzmin <- 1 # 1m grid spacing 
zmax <- 25 # lowermost node
base <- 1 # no grid coarsening with depth
Ti <--5

type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEUMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

st    <- st * 3600 * 24 # convert to [s]
t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

bc.upper <- seq(from = -1, to = 3, length.out = nt) #linear warming
bc.lower <- rep(0.05, nt)

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # with of node [m]
                   Tj = rep( Ti, nz), # initial temperature [C]
                   soi = rep(0.5, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(2e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

#== Target arrays to hold results =========
out.Tj  <- rbind(NULL,t(mat$Tj))
out.liq <- rbind(NULL,t(mat$liq))
out.t   <- t[0]
out.dz  <- mat$dz

HcGroundImplicit(dt, mat, bc.lower, bc.upper, unfrozen.type = "DALLAMICO",            # use dall'amico
                 unfrozen.par = c(1), bcutype = "NEUMANN", bcltype = "NEUMANN",
                 layers.sno = c(0), type.sno = "CALONNE", type.gnd = "COSENZA")

