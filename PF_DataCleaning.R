### Data Cleaning.R
###  Nick Brown, 2016
###
### Functions to clean up permafrost data

generateClip = function(minval,maxval=NA){
  ##Generates a clipping function which sets to NA values exceeding the limits 
  ## (if min and max are both specified).  Or sets to NA values which are equal to some
  ## other value.
  
  # Begin #
  if (!is.na(maxval)) {
  doClip = function(input){
    input[input>maxval|input<minval]<-NA
    input
  }
  }
  else if (is.na(maxval)){
    doClip=function(input){
      input[input==minval]<-NA
    }
  }
}



pointOutlierList = function(input,tolerance) {
    ## Generates a logical list corresponding to index data
    ## Entries are FALSE if they differ from each of their neighbours by more than the tolerance value
  
    #Begin#
  
    # Make shifted vectors with identical first or last values
    plus1 = c(input[1],input[-length(input)]) #doubled first
    min1 = c(input[-1],input[length(input)]) # doubled last
    
    #Compare input to neighbours
    P1check=(abs(input-min1)>tolerance)
    P2check=(abs(input-plus1)>tolerance)
    outlierTest = P1check&P2check
    
    #Return Logical vector 
    return(outlierTest)
}



generateRemoveOutliers = function(tolerance="sd") {
  ##Builds a function which removes data points which differ from their neighbours
  ##by more than a specified tolerance
  
  # Begin #
  doRemoveOutliers = function(input){
  if (tolerance == "sd") {
    tolerance<-sd(input,na.rm=T)
  }
  else if (tolerance !="sd"& is.numeric(tolerance)==F){
    tolerance<-sd(input,na.rm=T)
  }
  else {tolerance=tolerance}
  
  Out_list=pointOutlierList(input,tolerance)
  input[Out_list]<-NA
  return(input)
  }
}


