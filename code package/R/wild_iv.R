#this function computes a (wild) bootstrap p-value, imposing an null
wild_iv<-function(b_hypothezised,   #the hypothezised treatment effect under the null
                  x,y,z,            #the data x: endogeneous variable, y: outcome, z: instrument
                  r=1000,           #number of repetitions
                  seed=1234, mode=0) { #random seed and mode (mode is meaningless currently)

  #not sure if this is the approach that is also taken by Davidson & MacKinnon, but here I first impose the null and demean all vars
  #@Matt: does this make sense to you?
  set.seed(seed)

  y<-y-b_hypothezised*x #imposing the 0
  y<-y-mean(y) #"y1" in Davidson MacKinnon JBES 2010
  x<-x-mean(x) #"y2" in Davidson MacKinnon JBES 2010
  z<-z-mean(z) #"W"  in Davidson MacKinnon JBES 2010

  #compute the value of the test-statistic in the sample
  # I forgot why I did this after imposing the null, but now it seems to only work that way. odd? or somehow obvious?
  sample_test_statistic<-wild_iv_teststat(y,x,z)

  #structural equation: y = u1 (the null on b*x is already imposed)
  u1<-y

  #reduced form equation x = z*pi + u2
  Mzy <- lm(y ~ z)$residuals #see M_{z}Y in in eq (15) Davidson MacKinnon JBES 2010
  eq15_fit <- lm(x ~ z + Mzy) #equation (15) in Davidson MacKinnon JBES 2010
  u2<-x-eq15_fit$coefficients["z"]*z #based on x = z*. + u2

  #prepare the list where I will save all bootstrap t-statistics
  bs_test_statistic <- matrix(nrow=r)

  #bootstrap iterations
  for (b in 1:r) { #loops are inefficient. since I am going to re-do this in Rcpp anyways, I still use a loop for legibility
    #1: draw rademacher bs weights
    rademacher<-rbinom(n=length(y), size=1, prob=0.5)*2-1
    #2: create boostrap data
    y_bs<-rademacher*u1
    x_bs<-z*eq15_fit$coefficients["z"] + rademacher*u2
    #3: estimate bootstrap t-stat
    bs_test_statistic[b]<-wild_iv_teststat(y_bs,x_bs,z)
  }
  #return p-value
  return(mean(abs(bs_test_statistic)>abs(sample_test_statistic)))
}


#this function defines and computes what test-statistic is used by the WBS.
#currently we used the huber-white-robust t-stat
wild_iv_teststat<-function(y,x,z) {
  return(riiv::huber_white_tsls_tstats(as.matrix(y),as.matrix(x),as.matrix(z))[2])
}
