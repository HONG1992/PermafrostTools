###  Borehole Plotting.R
###  Nick Brown, 2016
###
### Functions for plotting borehole temperature profiles 
### Input data must have first column as time and the subsequent columns as data
### Depth values are found by searching the text string for a number, e.g. "X3.42" would return 3.42
### Time must be specified as a posix string



options(stringsAsFactors = FALSE)
try(library(xts),library(dygraphs))
#source("PF_DataCleaning.R")


boreholeTimeseries = function(dataset, dateformat="%d.%m.%Y",Ptz="UTC",
                              depths="all",daterange=c("",""),DYG=FALSE) {
  ### Function to plot time-series 
  ###  Must be 1st column time, all other columns depth
  
  # Begin #
  
  #Convert time to POSIX
  dataset[,1]=as.POSIXct(dataset[,1],tz=Ptz,format=dateformat)
  
  #Subset data based on date range.  Must use same POSIX format as dataset
  if (daterange[1]==""&daterange[2]==""){  #Default use all data
    dataset<-dataset
  }
  else {  #subset max/min
    daterange=as.POSIXct(daterange,tz=Ptz,format=dateformat)
    dataset<-dataset[dataset[,1]>=daterange[1]&dataset[,1]<=daterange[2],]
  } 
  
  #Split dataset into dates and data for convenience.  This should be fixed later 
  #  so as not to create extra memory load
  PFdates=as.POSIXct(dataset[,1],tz=Ptz,format=dateformat)
  PFdata=dataset[,-1]
  
  
  #Instantiate Generator for PF lines
  addPFdata = addlines_x(dataset)
  
  #Make list of depths to plot if specified, otherwise plot all depths.  Might be unneccesary
  if (depths[1] == "all") {
    Plot_depths <- rep(TRUE,length(PFdata))
    depths = seq(1,length(PFdata))
  }
  else {
    Plot_depths = rep(FALSE,length(PFdata))
    Plot_depths[depths]=TRUE
  }
  
  if (DYG==F) {
  #Make Plotting window
  plot(dataset[1,1],dataset[1,2],type="n",
       xlim=range(PFdates),ylim=c(min(PFdata[,Plot_depths],na.rm=T),max(PFdata[,Plot_depths],na.rm=T)),    # Dimensions
       main="Temperature Time Series", xlab="Date", ylab="Temperature"   # Labels
       )
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")  #gray background
  
  #Generate plotting function and run for specified depths
  DoPlot = addlines_x(dataset)
  invisible(lapply(depths,DoPlot))
  }
  else if (DYG==T){  #Option to use dynamic graphs package
    DYG_PF<-xts(x=dataset[,-1][,Plot_depths],order.by = dataset[,1])
    dygraph(DYG_PF)
  }
}
###  Addlines Function
addlines_x <- function(datasource) { 
  #makes a function which plots timeseries at certain depth for a dataset
  
  #Make a function which plots lines at a certain depth
  addlines = function(depth,plot.col="random"){
  
  #Get colour palette from specification, or randomize it
  if (plot.col=="random") {
    plot.col=colours()[sample(1:length(colours()),1)]
  }
  else{
    plot.col=plot.col
  }
  
  #
  depthname=names(datasource[,-1])[depth]
  data_plot=(datasource[,-1])[,depth]
  names(data_plot)=depthname
  lines(datasource[,1], data_plot, col = plot.col,lwd=2)
    }
}

indexDepth = function(PFdata){
  items<-data.frame(seq(1,length(PFdata)-1),names(PFdata[,-1]))
  names(items)<-c("Index","Depth")
  items
}

###
TrumpetPlot = function(data, plotyear,dateformat="%d.%m.%Y",Ptz="UTC",dmax=-1){
  
  ## Sanity Checks
  if (is.numeric(dmax)==FALSE){
    stop("Max depth (dmax) is not a valid numeric value")
  }
  
  #Convert time to POSIX
  data[,1]=as.POSIXct(data[,1],tz=Ptz,format=dateformat)
  
  #Get list of depths by using a string match 
  names(data)[-1]=as.numeric(regmatches(names(data),regexpr("\\d+\\.\\d+",names(data))))
  
  #Subset data
  data=data[as.POSIXlt(data[,1])$year+1900==plotyear,]
  
  #Build max, mean and min vectors, ignoring any non-finite values in the matrix
  PFmax  = sapply(data[,-1],function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA))
  PFmin  = sapply(data[,-1],function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA))
  PFmean = sapply(data[,-1],function(x) ifelse( !all(is.na(x)), mean(x, na.rm=T), NA))
  PFstat = data.frame(PFmax,PFmin,PFmean)

  
  #Establish Plotting boundaries
  if (dmax<0) {
    dmax<-c(-max(abs(as.numeric(rownames(PFstat))),na.rm=T),0)
  }
  
  else {
    dmax <- c(-dmax,0)
  }
  #Return output
  
  plot(PFmax,-as.numeric(names(PFmax)),type="l",
       ylim=dmax, xlim=range(PFmax,PFmin,na.rm=T),
       main=paste("Trumpet plot for ",plotyear,sep=""), xlab="Temperature (C)",ylab="Depth (m)"
       )
   lines(PFmin,-as.numeric(names(PFmax)))
   lines(PFmean,-as.numeric(names(PFmax)))
# return(PFstat)
}



