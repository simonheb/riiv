# load packages
#install.packages("AER")
library(AER)
#install.packages("mvtnorm")
library(mvtnorm)
#install.packages("dplyr")
library(dplyr)
#install.packages("MASS")
library(MASS)

library(riiv)

##############GENERAL SETTINGS
#the CI-level we care about today
ci_level<-.9
#sample size
n <- 200
treatment_effect <-2
mode<-0

##############Illustative Example
#simulate data for an analysis with an omitted variable bias
seed<-as.numeric(Sys.time())
print(seed)
set.seed(seed)
# start off with two correlated covariate
XandC <- mvrnorm(n, c(0,0), matrix(c(1, 0.95, 0.95, 1), 2, 2))
x_star <- XandC[, 1]
c <- XandC[, 2]

#there's also a binary instrument
z <- rbinom(n,1,0.5)
#the observed X combines the x_star and the weak instrument. (it's weak because of the 0.2)
x <- x_star + 0.3*z

# the outcome combines x, c and some noise
y <-  treatment_effect*x + 2*c + rnorm(n)

#as a benchmark, estimate the full information model
r1<-lm(y ~ x + c)
summary(r1)

#compute the naive (biased) estimator
r2<-lm(y ~ x)
summary(r2)

#compute the (weak) IV estimate
r3<-ivreg(y ~ x  | z)
summary(r3)






##############Repeat this 1000 times
cis<-data.frame(est="",lower_ci=0,upper_ci=0, stringsAsFactors=FALSE)
for (i in 1:5000) {
  if (i%%100==0)
    cat(i,"\n")
  cat(".")
  set.seed(i+seed+1)

  XandC <- mvrnorm(n, c(0,0), matrix(c(1, 0.95, 0.95, 1), 2, 2))
  x_star <- XandC[, 1]
  c <- XandC[, 2]
  z <- rbinom(n,1,0.5)
  x <- x_star + 0.2*z
  y <-  treatment_effect*x + 2*c + rnorm(n)

  r1<-lm(y ~ x + c)
  r2<-lm(y ~ x)
  r3<-ivreg(y ~ x  | z)
  r4<-inverted_test_gridsearch_cis(grid=seq(-50+r3$coefficients["x"],50+r3$coefficients["x"],length.out=7),
                                   y=y,x=x,z=z,level=ci_level,mode=mode,testfunction=ri_iv_test_p,r=99)

  reject_the_truth<-ri_iv_test_p(b=treatment_effect,y=y,x=x,z=z,r=999,mode=mode)

  r5<-r4#inverted_test_gridsearch_cis(grid=seq(-20+r3$coefficients["x"],20+r3$coefficients["x"],length.out=7),
         #                          y=y,x=x,z=z,level=ci_level,mode=mode,testfunction=wild_iv,r=999)

  reject_the_truth_wbs<-1#wild_iv(b_hypothezised=treatment_effect,y=y,x=x,z=z,r=999)

  cis<-rbind(cis,
             c("full information CI coverage",confint(r1,level=ci_level)["x",]),
             c("ignoring endogeneity CI (coverage)",confint(r2,level=ci_level)["x",]),
             c("ivreg CI (coverage)",confint(r3,level=ci_level)["x",]),
             c("RIIV-approx. CI (coverage)",r4$ci),
             c("RIIV-approx. CI infinite?",treatment_effect+c(-1e-2,1e-2)+ifelse(any(is.infinite(r4$ci)),1,0)),
             c("RIIV-reject the truth",treatment_effect+c(-1e-2,1e-2)+ifelse(reject_the_truth>1-ci_level,1,0)), #this is a pseudo CI, that contains the truth if we reject it and doesnt  if we can

             c("WBS-CI (coverage)",r5$ci),
             c("WBS-reject the truth", treatment_effect+c(-1e-2,1e-2)+ifelse(reject_the_truth_wbs>1-ci_level,1,0)) #this is a pseudo CI, that contains the truth if we reject it and doesnt  if we can
  )
  cis$upper_ci<-as.numeric(cis$upper_ci)
  cis$lower_ci<-as.numeric(cis$lower_ci)



  #summarize results
  if (i%%100==0)
    cis %>%
    #create indicator for whether the true TE is inside the CI
    mutate(covered=ifelse(lower_ci<treatment_effect  &  treatment_effect<upper_ci,1,0)) %>%
    #summarize and compute average by estimator
    group_by(est) %>%
    summarize(coverage=mean(covered))%>%
    print

}
#summarize results
cis %>%
  #create indicator for whether the true TE is inside the CI
  mutate(covered=ifelse(lower_ci<treatment_effect  &  treatment_effect<upper_ci,1,0)) %>%
  #summarize and compute average by estimator
  group_by(est) %>%
  summarize(coverage=mean(covered))%>%
  print

