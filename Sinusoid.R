
###Function which creates climate variables
Climate = function(Amplitude=4,MAGST=-0.5) {
  clim=list(Amplitude=Amplitude,
            MAGST=MAGST,
            units = list(Amplitude="C",MAGST="C"))  
  return(clim)
}

## Function which returns list describing field site
PF_site = function(heatflux = 0.2, heatcap=1.8e6,
                   th_cond=2.5, MST=-2) {
  sitechars=list(HeatFlux=c(heatflux),
                 HeatCapacity=c(heatcap),
                 ThermalConduct=c(th_cond),
                 MST=c(MST),
                 ThermDiff=th_cond/heatcap,
                 Units=list(HeatFlux="W/m2",
                            HeatCapacity="J/m3K",
                            ThermalConduct="W/mK",
                            MST="C",
                            ThermDiff="[unitless]") )
  return(sitechars)
}

## Function to generate temperatures for a given depth and time
SinudoidalVariation <-function(site,climate){
  #produces a closure for a site and climate

T_analytical = function(depth,time) {  
  if (length(depth)>1&length(time)>1){
    #build arrays
    z=array(rep(depth),length(time),dim=c(length(depth),length(time)))
    t=t(array(rep(time),length(depth),dim=c(length(time),length(depth))))
    
  }
  else{
    z=depth
    t=time
  }
  #creates a function to yield a temperature for a given depth and time
  t=t*60*60*24 # time to days
  z=-abs(z)
  T0=climate$MAGST - (site$HeatFlux/site$ThermalConduct)*z
  w=(2*pi)/(31536000)  #angular velocity
  constant = sqrt((w)/(2*site$ThermDiff)) # value of sqrt(w/2k)
  
  DelT = climate$Amplitude*exp(z*constant)*cos(w*t+z*constant)
  Tfin=DelT+T0
  rownames(Tfin) = depth
  colnames(Tfin) = paste("t",time,sep="")
  return(Tfin)
  #return(list(DelT=DelT,Amp=climate$Amplitude,z=z,cons=constant,om=w,time=t,cosfun=cos(w*t+z*constant)))
}
}


## Build a data frame -- obsolete after inclusion of array capabilities in SinusoidalVariation()
buildDF = function(TempFunc,depths,timestep,numsteps){
  output <- NULL
  Cnames <- NULL
  times <- seq(0,numsteps)*timestep
  for (i in seq(0,numsteps)){
    colm<-sapply(depths,function(x) TempFunc(z=x,t=i*timestep))
    Cnames <- c(Cnames,paste("t",i*timestep,sep=""))
    output <- cbind(output,colm)
  }
  colnames(output)=Cnames
  return(output)
}
