
setwd("C:\\Users\\Nick\\Code\\PermafrostTools\\Assignment 2")
source("RequiredFunctions.R")

runQ2 <- function(years,kthaw,kfactor){
years.to.run <- years  # run for number of years
k.factor <- kfactor
#follow the example from the help file, with a few exceptions

#Set parameters
dt <- toSeconds(1,"d") # 1 day time steps
ns <- 1  # save every 1 days
st <- toSeconds(years.to.run,"y")  # [s] run it for number of years specified above
k.thawed <- kthaw
k.frozen <- k.thawed*k.factor

dzmin <- .001 # 10 cm grid spacing 
zmax <- 5 # lowermost node
base <- 1.1 #  grid coarsening with depth
Ti <- 0

t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations


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

#Boundary Conditions

bc.upper <- yearly.T(mean=Ti,6,years=years.to.run)$T.year
bc.lower <- rep(0, nt)  #  Zero heat flux at depth

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # width of node [m]
                   Tj = rep(Ti,nz), # initial temperature [C]
                   soi = rep(.9, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0.1, nz), # volumetric proportion of total water [0-1]
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
  
  ##  At the end of each time-step, adjust the soil properties accordingly
  mat$k.soi[mat$Tj<0] <- k.frozen   
  mat$k.soi[mat$Tj>=0] <- k.thawed
#   mat$kj <- mat$k.soi
#   mat$liq <- 0
#   mat$ice <- 0
#   mat$cj <- mat$c.soi
  
  #output only in specified interval
  check <- i/ns - floor(i/ns)
  if (check == 0) {
    out.Tj  <- rbind( out.Tj,t(mat$Tj))
    out.liq <- rbind(out.liq,t(mat$liq))
    out.t   <- t[i]
    print(paste(i,nt))
  }
}
return(list(mat=mat,out.Tj=out.Tj, bc=c(bc.upper,bc.lower)))
}


##  Run the spin-up experiments

firstTest <- runQ2(10,1,2.5)

#Examine convergence of annual temperature plots
dev.new()
plot(1,1,type="n",xlim=c(-5,2),ylim=c(-5,0),
     xlab="End-of-winter Temperature (C)",ylab="Depth (m)"
#     main="Figure 2-1"
     )
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
lines(firstTest$out.Tj[1,],firstTest$mat$zj,type="l") #1
lines(firstTest$out.Tj[365,],firstTest$mat$zj,type="l",col="blue") #2
lines(firstTest$out.Tj[365*2,],firstTest$mat$zj,type="l",col="green")
lines(firstTest$out.Tj[365*3,],firstTest$mat$zj,type="l",col="red",lty=1)
lines(firstTest$out.Tj[365*4,],firstTest$mat$zj,type="l",col="orange",lty=1)
lines(firstTest$out.Tj[365*5,],firstTest$mat$zj,type="l",col="yellow",lty=1)
lines(firstTest$out.Tj[365*6,],firstTest$mat$zj,type="l",col="black",lty=2)
lines(firstTest$out.Tj[365*7,],firstTest$mat$zj,type="l",col="blue",lty=2)
lines(firstTest$out.Tj[365*8,],firstTest$mat$zj,type="l",col="green",lty=2)
lines(firstTest$out.Tj[365*9,],firstTest$mat$zj,type="l",col="red",lty=2)
lines(firstTest$out.Tj[365*10,],firstTest$mat$zj,type="l",col="orange",lty=2)

cols21 = c("black","blue","green","red","orange","yellow","black","blue","green","red","orange")
years =c("Year 1","Year 2","Year 3","Year 4","Year 5","Year 6","Year 7","Year 8","Year 9","Year 10")
ty21 = c(1,1,1,1,1,1,2,2,2,2,2)

legend(0.25,-2,legend=years,col=cols21,lty=ty21, ## add legend
       title ="Year")


## Examine convergence of annual mean temperatures
dev.new()
plot(1,1,type="n",xlim=c(-3,2),
     ylim=c(-5,0),
     xlab="Mean Annual Temperature (C)",ylab="Depth (m)"
#     main="Figure 2-2"
     )
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
lines(sapply(as.data.frame(firstTest$out.Tj[0:365,]),mean),firstTest$mat$zj,type="l")
lines(sapply(as.data.frame(firstTest$out.Tj[366:(365*2),]),mean),firstTest$mat$zj,col="blue")  # mean for year 2
lines(sapply(as.data.frame(firstTest$out.Tj[(365*2):(365*3),]),mean),firstTest$mat$zj,col="green") # mean for year 3
lines(sapply(as.data.frame(firstTest$out.Tj[(365*3):(365*4),]),mean),firstTest$mat$zj,col="red") # mean for year 4
lines(sapply(as.data.frame(firstTest$out.Tj[(365*4):(365*5),]),mean),firstTest$mat$zj,col="orange")
lines(sapply(as.data.frame(firstTest$out.Tj[(365*5):(365*6),]),mean),firstTest$mat$zj,col="yellow")
lines(sapply(as.data.frame(firstTest$out.Tj[(365*6):(365*7),]),mean),firstTest$mat$zj,col="black",lty=2)
lines(sapply(as.data.frame(firstTest$out.Tj[(365*7):(365*8),]),mean),firstTest$mat$zj,col="blue",lty=2)
lines(sapply(as.data.frame(firstTest$out.Tj[(365*8):(365*9),]),mean),firstTest$mat$zj,col="green",lty=2)
lines(sapply(as.data.frame(firstTest$out.Tj[(365*9):(365*10),]),mean),firstTest$mat$zj,col="red",lty=2)



cols21 = c("black","blue","green","red","orange","yellow","black","blue","green","red")
years =c("Year 1","Year 2","Year 3","Year 4","Year 5","Year 6","Year 7","Year 8","Year 9")
ty21 = c(1,1,1,1,1,1,2,2,2,2)

legend(1,-2,legend=years,col=cols21,lty=ty21, ## add legend
       title ="Year")


## look at temperature with depth
dev.new()
plot(seq(1,20,length.out = dim(firstTest$out.Tj)[1]),
     firstTest$out.Tj[,65],type="l",xlab="Years",ylab="Temperature")

#### Run experiment for multiple factors
##  runQ2(years,k_thawed,thawed:frozen factor)

factor1 <- runQ2(10,1,1)
factor2 <- runQ2(10,1,2)
factor3 <- runQ2(10,1,3)
factor4 <- runQ2(10,1,4)
factor5 <- runQ2(10,1,5)

L <- length(factor1$out.Tj[,1])  # number of time steps in model
Lm <- 365 # number of steps in a year

## Calcualte mean temperatures for the final year of the model
F1temp <-sapply(as.data.frame(factor1$out.Tj[(L-Lm):L,]),mean)
F2temp <-sapply(as.data.frame(factor2$out.Tj[(L-Lm):L,]),mean)
F3temp <-sapply(as.data.frame(factor3$out.Tj[(L-Lm):L,]),mean)
F4temp <-sapply(as.data.frame(factor4$out.Tj[(L-Lm):L,]),mean)
F5temp <-sapply(as.data.frame(factor5$out.Tj[(L-Lm):L,]),mean)

##  Plot mean annual temperature plots 
Fig4.15 <- as.array(cbind(F1temp,F2temp,F3temp,F4temp,F5temp))
depths <- array(rep(factor1$mat$zj,5),dim=c(length(factor1$mat$zj),5))

##  Plot
dev.new()
plot(Fig4.15[,1],depths[,1],type="n",xlim=c(-1.5,0),ylim=c(-5,0),
     xlab="Mean Annual Temperature (C)", ylab="Depth (m)") # make plot
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")  # add background
ratio = c(1,2,3,4,5)
plotcols = c("black","red","blue","yellow","pink")
plottyp = c(1,2,3,4,5)

lapply(c(1:5),function(x) lines(Fig4.15[,x],depths[,x],col=plotcols[x],lty=plottyp[x]))  ## add lines

legend(0.25,-2,legend=as.character(ratio),col=plotcols,lty=plottyp, ## add legend
      title ="Frozen:Thawed Conductivity Ratio")


