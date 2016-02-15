
## Download necessary packages
list.of.packages <- c("pracma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Load Libraries
library(PermafrostTools)
library(pracma)

## Functions

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

## Analytic Heat step change


# Convert time to seconds
toSeconds <- function(s,units="s",leap="FALSE") {
  lookup <- c(1,60,3600,3600*24,2628000,3600*24*365)
  names(lookup)<-c("s","m","h","d","m","y")
  if (leap==TRUE) { lookup[which(names(lookup)=="y")]<-3600*24*366 }
  if (leap=="average") { lookup[which(names(lookup)=="y")]<-3600*24*365.25 }
  output <- lookup[which(names(lookup)==units)]*s
  return(output)
} 


Stepchange <- function(site,climate,Tstep,zbound=c(0,60), zstep=1) {   
  depths <- -seq(zbound[1],zbound[2],zstep)
  Ti <- climate$MAGST-depths*(site$HeatFlux/site$ThermalConduct)
  output<-cbind(depths, Ti)
  
  addstep <- function(time=NULL, unit="s"){
    if (is.numeric(time)) {
      time.s = toSeconds(time,unit)
      paren <- sqrt(site$ThermDiff*time.s)
      DelT <- Tstep*erfc(-depths/(2*paren))
      add <-(Ti + DelT)
      output <<- as.data.frame(cbind(output, add))
      names(output)[which(names(output)=="add")]<<-paste(time,unit,sep="")
      #output<-list(output,add)
    }
  }
  getvalues <- function() {
    return(output)
  }
  return(list(add.step=addstep,get.values=getvalues,site.params=site,climate.params=climate))
}

## Build surface temperature for 1 year
yearly.T <- function(mean, amplitude, warming = 0, years = 1) {
  jul.day <- c(0:(365*years))
  T.year <- mean+amplitude*sin(jul.day*2*pi/365) + (jul.day*warming/365) #(1-erfc(2*B$jul.day/(365*years)))*warming
  return(data.frame(jul.day,T.year))
}