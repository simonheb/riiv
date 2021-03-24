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



##############Illustative Example
#simulate data for an analysis with an omitted variable bias
set.seed(0)
# start off with two correlated covariate
XandC <- mvrnorm(n, c(0,0), matrix(c(1, 0.95, 0.95, 1), 2, 2))
x_star <- XandC[, 1]
c <- XandC[, 2]

#there's also a binary instrument
z <- rbinom(n,1,0.5)
#the observed X combines the x_star and the weak instrument. (it's weak because of the 0.2)
x <- x_star + 0.2*z

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

#compute the [confidence band] using RI and plot
r4<-ri_iv_gridsearch(grid=seq(-20+r3$coefficients["x"],20+r3$coefficients["x"],length.out=11),
               y=y,x=x,z=z,level=ci_level,refinementsteps=6,wilcox=F,allowInf=T,plot=T)
print(r4$ci)







##############Repeat this 1000 times
cis<-data.frame(est="",lower_ci=0,upper_ci=0, stringsAsFactors=FALSE)
for (i in 1:1000) {
  set.seed(i)

  XandC <- mvrnorm(n, c(0,0), matrix(c(1, 0.95, 0.95, 1), 2, 2))
  x_star <- XandC[, 1]
  c <- XandC[, 2]
  z <- rbinom(n,1,0.5)
  x <- x_star + 0.2*z
  y <-  treatment_effect*x + 2*c + rnorm(n)

  r1<-lm(y ~ x + c)
  r2<-lm(y ~ x)
  r3<-ivreg(y ~ x  | z)
  r4<-ri_iv_gridsearch(grid=seq(-20+r3$coefficients["x"],20+r3$coefficients["x"],length.out=11),
                       y=y,x=x,z=z,level=ci_level,refinementsteps=6,wilcox=F,allowInf=T,plot=F)

  cis<-rbind(cis,
            c("full information",confint(r1,level=ci_level)["x",]),
            c("naive",confint(r2,level=ci_level)["x",]),
            c("weak instrument | ivreg",confint(r3,level=ci_level)["x",]),
            c("weak instrument | RIIV",r4$ci))

}
#summarize results
    cis %>%
      #create indicator for whether the true TE is inside the CI
      mutate(covered=ifelse(lower_ci<treatment_effect  &  treatment_effect<upper_ci,1,0)) %>%
      #summarize and compute average by estimator
      group_by(est) %>%
      summarize(coverage=mean(covered))%>%
      print



