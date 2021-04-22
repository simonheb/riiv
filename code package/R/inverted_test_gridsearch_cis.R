#this function conducts a grid search using ri test to find CIs

inverted_test_gridsearch_cis <-function(
  y,x,z,               #data
  grid,                #prespecified search grid optional
  allowInf=TRUE,          #return Inf for upper and lower CI bound if the grid-ranges cannot be rejected
  level=0.9,           #ci-level
  plot=FALSE,          #plot results?
  refinementsteps=5,   #how many recursive refinement steps?
  refinementsresolution=7,            #steps to divide this into
  seed=1234,           #random seed
  gridres=NULL,        #results
  r=1000,              #RI iterations
  verbose=FALSE,        #debugging
  testfunction=ri_iv_test_p, #this function is defined in src/ri_iv_test_cpp.cpp
  ...
) {
  #loop over values in grid and compute RI-p-values
  for (b in grid) {
    #copute p-value
    pval<-testfunction(b=b,y=y,x=x,z=z,r=r,seed=seed,...)
    #append to grid results
    gridres<-rbind(gridres,c(b,pval))
  }
  #sort by value
  gridres<-gridres[order(gridres[,1]),]
  #round to multiples of 1/r (other values should not be possible anyways but otherwise numeric imprecision can cause problems)
  gridres[,2]<-round(gridres[,2]*r)/r

  #we now do a number of sanity checks on the results

  #which values were we not able to reject?
  unable_to_reject<-which(gridres[,2]>1-level)
  #initialize indicator to break if an erros if ound
  stopnow<-FALSE

  if (length(unable_to_reject)==0) {
    stopmsg<-c("Rejected all values. Redefine grid!\n(probably it needs to be finer or elsewhere)")
    stopnow<-TRUE
  }
  else if (length(unable_to_reject)==nrow(gridres)) {
    stopmsg<-c("Rejected no values. Redefine grid, or accept that the CI is infinite.\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      warning(stopmsg)
      gridres[1,1]<- -Inf
      gridres[nrow(gridres),1]<- Inf
      gridres[1,2]<- NA
      gridres[nrow(gridres),2]<- NA
    } else {
      stopnow<-TRUE
    }
  }
  else if (all(c(1,nrow(gridres)) %in% unable_to_reject)) {
    stopmsg<-c("Can neither reject upper nor lower bound of grid. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      warning(stopmsg)
      gridres[1,1]<- -Inf
      gridres[nrow(gridres),1]<- Inf
      gridres[1,2]<- NA
      gridres[nrow(gridres),2]<- NA
    } else {
      stopnow<-TRUE
    }
  }
  else if (min(unable_to_reject)==1) {
    stopmsg<-c("Cannot reject lower bound. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      gridres[1,1]<- -Inf
      gridres[1,2]<- NA
      warning(stopmsg)
    } else {
      stopnow<-TRUE
    }
  }
  else if (max(unable_to_reject)==nrow(gridres)){
    stopmsg<-c("Cannot reject upper bound. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      gridres[nrow(gridres),1]<- Inf
      gridres[nrow(gridres),2]<- NA
      warning(stopmsg)
    } else {
      stopnow<-TRUE
    }
  }

  #If required produce some output and then end execution on error.
  if (stopnow) {
    if (verbose){
      plot(gridres)
      print(gridres)
    }
    stop(stopmsg)
  }

  #plot if desired
  if (plot & verbose) {
    plotdata<-as.data.frame(gridres)
    names(plotdata)<-c("Treatment Effect","p-value")
    plot(plotdata)
  }
  #if refinementsteps>0, create finer grid around the current estimates for the CI and recursively call this function to test these values
  if (refinementsteps) {
    #initialize empty new grid values
    refined_grid_upper_bound<-refined_grid_lower_bound<-NULL
    #if lower bound is finite:
    if (!is.infinite(gridres[1,1])) {
      refined_grid_lower_bound<-seq(gridres[min(unable_to_reject)-1,1],gridres[min(unable_to_reject),1],length.out=ceiling(refinementsresolution/2)+2)
      refined_grid_lower_bound<-refined_grid_lower_bound[2:(length(refined_grid_lower_bound)-1)]
    }
    #if upper bound is finite:
    if (!is.infinite(gridres[nrow(gridres),1])) {
      refined_grid_upper_bound<-seq(gridres[max(unable_to_reject),1],gridres[max(unable_to_reject)+1,1],length.out=ceiling(refinementsresolution/2)+2)
      refined_grid_upper_bound<-refined_grid_upper_bound[2:(length(refined_grid_upper_bound)-1)]
    }

    #combine to grid
    refined_grid<-c(refined_grid_upper_bound,refined_grid_lower_bound  )

    if (length(refined_grid)>0) {
      #call function
      return(inverted_test_gridsearch_cis(grid=refined_grid,y=y,x=x,z=z,level=level,plot=plot,refinementsteps=refinementsteps-1,allowInf=allowInf,r=r,seed=seed,refinementsresolution=refinementsresolution,gridres=gridres,verbose=verbose,testfunction,...))
    }
  }

  #use gridres to extract CI values, take extra care of Inf and -Inf
  if (is.infinite(gridres[1,1]) & !is.infinite(gridres[nrow(gridres),1]))
    ci<-c(-Inf,gridres[max(unable_to_reject)+1,1])
  else if (!is.infinite(gridres[1,1]) & is.infinite(gridres[nrow(gridres),1]))
    ci<-c(gridres[min(unable_to_reject)-1,1],Inf)
  else if (is.infinite(gridres[1,1]) & is.infinite(gridres[nrow(gridres),1]))
    ci<-c(-Inf,Inf)
  else
    ci<-c(gridres[min(unable_to_reject)-1,1],gridres[max(unable_to_reject)+1,1])

  #plot if desired
  if (plot) {
    plotdata<-as.data.frame(gridres)
    names(plotdata)<-c("Treatment Effect","p-value")
    plot(plotdata)
    lines(ci,rep(1-level,2),col=2)
  }
  return(list(ci=ci,gridres=gridres))
}
