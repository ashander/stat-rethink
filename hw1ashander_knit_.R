
## @knitr NA
require(bbmle)
require(rethinking)
require(ggplot2)


## @knitr problem1_setup
data(homeworkch2) #


## @knitr prob1a
 p.b <- seq(from=0, to=1, length.out=1000)
 prior <- rep(1/1000, 1000)
 likelihood <- dbinom(sum(birth1)+sum(birth2), size=length(c(birth1,birth2)), prob=p.b) # likelihood of data given 1000 models (binomial success parameter)
 posterior <- prior * likelihood
 posterior <- posterior/sum(posterior)
Max.post <- function(parameter, posterior){
	estimate=parameter[-log(posterior) == min(-log(posterior))]
	estimate
}
 pb.max <- Max.post(p.b, posterior) 


## @knitr prob1a
print(pb.max)


## @knitr prob1b-fxns

# construct table like precis
Precis.bayes<-function (parameter, posterior, param.samples=NA, level=0.95)
{
  # param.samples = an set of samples TODO: change to do resampling ?
  result <- data.frame(est=Max.post(parameter, posterior))
  colnames(result) <- c("Estimate")
  if (length(param.samples)>10){
        ci <- HPDI(param.samples, level)
        result <- cbind(result, as.data.frame(t(ci)))
    }
   result
}


## @knitr prob1b
p.b.sample <- sample(p.b, size=1e4, replace=TRUE, prob=posterior)

CI.data <- t(sapply(CI.types,
                  function(x){unlist(Precis.bayes(p.b, posterior=posterior, param.samples=p.b.sample, level=x))}
         ))
fit.bayes <- as.data.frame(CI.data)
fit.bayes$Interval.width <- as.character(CI.types)


## @knitr prob1btab
fit.bayes


## @knitr prob1c
births <- c(birth1, birth2)
pb.ml <-  mle2(y ~ dbinom(size=length(births), prob=pb), data=list(y=sum(births)), start=list(pb=0.5))
CI.dat.bayes <-  t(sapply(CI.types,
                        function(x){unlist(precis(pb.ml, level=x))}
           ))
fit.ml <- as.data.frame(CI.dat.bayes)
names(fit.ml)[3:4] <- c("lower", "upper")
fit.ml$Interval.width <- as.character(CI.types)
fit.ml$`Std. Error` <- NULL


## @knitr prob1ctab
fit.ml


## @knitr prob1_fig2
dens.compare <- data.frame(probability=p.b, Bayes=-log(posterior), ML = -log(dnorm(p.b, mean=pb.max, sd=sqrt(vcov(pb.ml)))/sum(dnorm(p.b, mean=pb.max, sd=sqrt(vcov(pb.ml))))))
d.c = melt(dens.compare, id.vars='probability')

g <- ggplot(d.c, aes(probability, value, color=variable))
g+geom_line()+ylab("- log density")


## @knitr prob1_fig1
fit.bayes$type <- "Bayes (HDPI)"
fit.ml$type <- "Max. Lik. (CI)"
estimate.ml <- fit.ml[,1] # recenter by the ML estimates
fit.ml[,1:3] <- fit.ml[,1:3] - estimate.ml
fit.bayes[,1:3] <- fit.bayes[,1:3] - estimate.ml
fit.both <- rbind(fit.bayes, fit.ml)
fit.both$Interval.width <- as.factor(fit.both$Interval.width)

g <- ggplot(fit.both)
g + geom_pointrange(aes(Interval.width, Estimate, ymin=lower, ymax=upper, color=type), position=position_dodge(width=0.1))


## @knitr prob2a
# define a function to use in simulations 
Num.boys <- function(p.b.this, n.total){
  births.this <- rbinom(n=n.total, size=1, prob=p.b.this)
  return(sum(births.this))
}

p.b.sims  <- sample(p.b, size=1e4, replace=TRUE, prob=posterior)
boys.sims <- sapply(p.b.sims, function(x){ Num.boys(x, 200)})
dens(boys.sims, adj=1)
abline(v=sum(c(birth1,birth2)),lwd=3, col='red')


## @knitr prob2b
p.b.sims  <- sample(p.b, size=1e4, replace=TRUE, prob=posterior)
fb.boys.sims <- sapply(p.b.sims, function(x){ Num.boys(x, 100)})
dens(fb.boys.sims, adj=1)
abline(v=sum(c(birth1)),lwd=3, col='red')


## @knitr prob2c
birth.pg <- birth2[birth1==0]
cat("first born girls: ", length(birth.pg), ", boys born after girls: ", sum(birth.pg), "\n")


## @knitr prob2c
p.b.sims  <- sample(p.b, size=1e4, replace=TRUE, prob=posterior)
pg.sims <- sapply(p.b.sims, function(x){ Num.boys(x, 49)})
dens(pg.sims, adj=1)
abline(v=sum(c(birth.pg)),lwd=3, col='red')


## @knitr unnamed-chunk-1
require(knitr) ## the package
opts_knit$set(base.url="https://github.com/ashander/stat-rethink/raw/master/")


## @knitr unnamed-chunk-2
opts_knit$set(out.format='gfm')
knit('/Users/jaime/PHD/stat-rethink/hw1ashander_knit_.md', tangle=TRUE) ## to run


