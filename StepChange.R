## Analytic Heat step change
library(pracma)

# Convert time to seconds
toSeconds <- function(s,units="s",leap="FALSE") {
  lookup <- c(1,60*60,3600*24,3600*24*365)
  names(lookup)<-c("s","h","m","y")
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

Ott<-Stepchange(Ottawa,OtClim,+3)


f=c(1,2,4,8,15)
lapply(f,function(x) Ott$add.step(x,"y"))
dat1<-Ott$getValues

plot(dat1$Ti,dat1$depths,type="l",xlim=c(-10,10))
lines(dat1[,3],dat1$depths)
lines(dat1[,4],dat1$depths)
lines(dat1[,5],dat1$depths)
lines(dat1[,6],dat1$depths)
lines(dat1[,7],dat1$depths)


