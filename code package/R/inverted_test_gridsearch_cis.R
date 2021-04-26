#this function conducts a grid search using ri test to find CIs.
#it recursively refines the grid between the likely upper and lower estimate of each end of the CI.
# * refinementsteps defines of many recursive steps are done
# * refinementresolution defines how coarse the subgrid in each step is going to be.
inverted_test_gridsearch_cis <-function(
  y,x,z,               #data
  grid,                #prespecified search grid optional
  allowInf=TRUE,          #return Inf for upper and lower CI bound if the grid-ranges cannot be rejected
  level=0.9,           #ci-level
  plot=FALSE,          #plot results?
  refinementsteps=5,   #how many recursive refinement steps?
  refinementsresolution=7,            #steps to divide this into
  seed=1234,           #random seed
  grid_values_and_corresponding_p_values=NULL,        #results from possible previous steps
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
    grid_values_and_corresponding_p_values<-rbind(grid_values_and_corresponding_p_values,c(b,pval))
  }

  #sort by value
  grid_values_and_corresponding_p_values<-grid_values_and_corresponding_p_values[order(grid_values_and_corresponding_p_values[,1]),]

  #round to multiples of 1/r (other values should not be possible but numeric imprecision can cause problems)
  grid_values_and_corresponding_p_values[,2]<-round(grid_values_and_corresponding_p_values[,2]*r)/r

  #which values were we not able to reject?
  unable_to_reject<-which(grid_values_and_corresponding_p_values[,2]>1-level)

  #we now do a number of sanity checks on the results. This basically tries to identify if the interval is going to be open/infinite/closed/...
  #this function also replaces the lowest [highest] grid-value with -Inf [Inf], if it suspects that the interval is not closed.
  grid_values_and_corresponding_p_values<-grid_sanity_checks(grid_values_and_corresponding_p_values,allowInf=allowInf,verbose=verbose,level=level)

  #if refinementsteps>0, create finer grid around the current estimates for the CI-bounds and recursively call this function to test these values
  if (refinementsteps) {
    #initialize empty new grid values
    refined_grid_upper_bound<-refined_grid_lower_bound<-NULL
    #if lower bound is finite:
    if (!is.infinite(grid_values_and_corresponding_p_values[1,1])) {
      refined_grid_lower_bound<-seq(grid_values_and_corresponding_p_values[min(unable_to_reject)-1,1],grid_values_and_corresponding_p_values[min(unable_to_reject),1],length.out=ceiling(refinementsresolution/2)+2)
      refined_grid_lower_bound<-refined_grid_lower_bound[2:(length(refined_grid_lower_bound)-1)]
    }
    #if upper bound is finite:
    if (!is.infinite(grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1])) {
      refined_grid_upper_bound<-seq(grid_values_and_corresponding_p_values[max(unable_to_reject),1],grid_values_and_corresponding_p_values[max(unable_to_reject)+1,1],length.out=ceiling(refinementsresolution/2)+2)
      refined_grid_upper_bound<-refined_grid_upper_bound[2:(length(refined_grid_upper_bound)-1)]
    }

    #combine to new grid-values
    refined_grid<-c(refined_grid_upper_bound,refined_grid_lower_bound  )

    if (length(refined_grid)>0) {
      #call self recursively
      return(inverted_test_gridsearch_cis(grid=refined_grid,y=y,x=x,z=z,level=level,plot=plot,refinementsteps=refinementsteps-1,allowInf=allowInf,r=r,seed=seed,refinementsresolution=refinementsresolution,grid_values_and_corresponding_p_values=grid_values_and_corresponding_p_values,verbose=verbose,testfunction,...))
    }
  }

  #based on grid_values_and_corresponding_p_values to extract CI values, take extra care of Inf and -Inf
  if (is.infinite(grid_values_and_corresponding_p_values[1,1]) & !is.infinite(grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]))
    ci<-c(-Inf,grid_values_and_corresponding_p_values[max(unable_to_reject)+1,1])
  else if (!is.infinite(grid_values_and_corresponding_p_values[1,1]) & is.infinite(grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]))
    ci<-c(grid_values_and_corresponding_p_values[min(unable_to_reject)-1,1],Inf)
  else if (is.infinite(grid_values_and_corresponding_p_values[1,1]) & is.infinite(grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]))
    ci<-c(-Inf,Inf)
  else
    ci<-c(grid_values_and_corresponding_p_values[min(unable_to_reject)-1,1],grid_values_and_corresponding_p_values[max(unable_to_reject)+1,1])

  #plot if desired
  if (plot) {
    #plot p-values
    plot(setNames(as.data.frame(grid_values_and_corresponding_p_values),c("Treatment Effect","p-value")))
    #plot line to show the estimated CI
    pseudoci=ci
    if (pseudoci[1]==-Inf)
      pseudoci[1]<-min(grid_values_and_corresponding_p_values[which(grid_values_and_corresponding_p_values!=-Inf)])
    if (pseudoci[2]==Inf)
      pseudoci[2]<-max(grid_values_and_corresponding_p_values[which(grid_values_and_corresponding_p_values!=Inf)])
    lines(pseudoci,rep(1-level,2),col=4)
    text(pseudoci[2],1-level,paste0(100*(1-level),"% CI",sep=""),pos=4,col=4)
  }
  #return
  return(list(ci=ci,grid_values_and_corresponding_p_values=grid_values_and_corresponding_p_values))
}


grid_sanity_checks<-function(grid_values_and_corresponding_p_values,allowInf,verbose,level) {

  #which values were we not able to reject?
  unable_to_reject<-which(grid_values_and_corresponding_p_values[,2]>1-level)
  #initialize indicator to break if an erros if ound
  stopnow<-FALSE

  if (length(unable_to_reject)==0) {
    stopmsg<-c("Rejected all values. Redefine grid!\n(probably it needs to be finer or elsewhere)")
    stopnow<-TRUE
  }
  else if (length(unable_to_reject)==nrow(grid_values_and_corresponding_p_values)) {
    stopmsg<-c("Rejected no values. Redefine grid, or accept that the CI is infinite.\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      warning(stopmsg)
      grid_values_and_corresponding_p_values[1,1]<- -Inf
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]<- Inf
      grid_values_and_corresponding_p_values[1,2]<- NA
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),2]<- NA
    } else {
      stopnow<-TRUE
    }
  }
  else if (all(c(1,nrow(grid_values_and_corresponding_p_values)) %in% unable_to_reject)) {
    stopmsg<-c("Can neither reject upper nor lower bound of grid. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      warning(stopmsg)
      grid_values_and_corresponding_p_values[1,1]<- -Inf
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]<- Inf
      grid_values_and_corresponding_p_values[1,2]<- NA
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),2]<- NA
    } else {
      stopnow<-TRUE
    }
  }
  else if (min(unable_to_reject)==1) {
    stopmsg<-c("Cannot reject lower bound. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      grid_values_and_corresponding_p_values[1,1]<- -Inf
      grid_values_and_corresponding_p_values[1,2]<- NA
      warning(stopmsg)
    } else {
      stopnow<-TRUE
    }
  }
  else if (max(unable_to_reject)==nrow(grid_values_and_corresponding_p_values)){
    stopmsg<-c("Cannot reject upper bound. Redefine grid, or accept that the CI is infinite!\n(Possibly the grid needs to be wider)")
    if (allowInf) {
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),1]<- Inf
      grid_values_and_corresponding_p_values[nrow(grid_values_and_corresponding_p_values),2]<- NA
      warning(stopmsg)
    } else {
      stopnow<-TRUE
    }
  }

  #If required produce some output and then end execution on error.
  if (stopnow) {
    if (verbose){
      plot(grid_values_and_corresponding_p_values)
      print(grid_values_and_corresponding_p_values)
    }
    stop(stopmsg)
  }
  return(grid_values_and_corresponding_p_values)
}
