# Homework 1, Statistical rethinking
 Jaime Ashander






# Problem 1  

First load the data, `birth1` and `birth2`, which indicate the gender (male=1, female=0) of reported first (second) born children in 100 families.
All families in this country have a maximum of two children.

Our hypothesis is that births are determined by a random process, with some probability of obtaining a boy `p.b`.
Assuming births are independent samples of this process, we model the data as random binomial trials.
 
```r
data(homeworkch2)  #
```



## 1 a



```r
p.b <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1/1000, 1000)
likelihood <- dbinom(sum(birth1) + sum(birth2), 
    size = length(c(birth1, birth2)), prob = p.b)  # likelihood of data given 1000 models (binomial success parameter)
naive.posterior <- prior * likelihood
naive.posterior <- naive.posterior/sum(naive.posterior)
Max.post <- function(parameter, posterior) {
    estimate = parameter[-log(posterior) == min(-log(posterior))]
    estimate
}
pb.max <- Max.post(p.b, naive.posterior)
```



  The maximimum posterior probability occurs with model: `p.b = `

```
## [1] 0.5546
```




## 1 b




To construct credible intervals, we sample the posterior, order the samples, and choose appropriate elements of from the vector of ordered samples.
Choice of appropriate sample and formatting of output is encapsulated in function `Precis.bayes` (a modified version of `rethinking::precis`). 
  
```r
p.b.sample <- sample(p.b, size = 10000, replace = TRUE, 
    prob = naive.posterior)
p.b.sample <- p.b.sample[order(p.b.sample)]

CI.types = c(0.5, 0.9, 0.95)
CI.data = t(sapply(CI.types, function(x) {
    unlist(Precis.bayes(p.b, posterior = naive.posterior, 
        param.ci = p.b.sample, level = x))
}))
fit.bayes = as.data.frame(CI.data)
names(fit.bayes)[2:3] = c("lower", "upper")
fit.bayes$Interval.width = as.character(CI.types)
```



Credible intervals, based on sampling the posterior density, for probabilty of drawing a boy:
  
```
##   Estimate  lower  upper Interval.width
## 1   0.5546 0.5305 0.5776            0.5
## 2   0.5546 0.4975 0.6116            0.9
## 3   0.5546 0.4855 0.6236           0.95
```



  Maximum posterior probability estimate of `p.b` is also reported.

## 1c

For maximum likelihood estimation we use `mle2`, and `rethinking::precis` to construct confidence intervals. 

```r
births = c(birth1, birth2)
pb.ml = mle2(y ~ dbinom(size = length(births), 
    prob = pb), data = list(y = sum(births)), start = list(pb = 0.5))
```
```
## Warning message: NaNs produced
```
```
## Warning message: NaNs produced
```
```
## Warning message: NaNs produced
```
```r
CI.dat.bayes = t(sapply(CI.types, function(x) {
    unlist(precis(pb.ml, level = x))
}))
fit.ml = as.data.frame(CI.dat.bayes)
names(fit.ml)[3:4] = c("lower", "upper")
fit.ml$Interval.width = as.character(CI.types)
fit.ml$`Std. Error` = NULL
```



  Confidence intervals, based on quadratic approximation of the likelihood surface, and maximum likelihood estimate of `p.b`, probabilty of drawing a boy:
  
```
##   Estimate  lower  upper Interval.width
## 1    0.555 0.5313 0.5787            0.5
## 2    0.555 0.4972 0.6128            0.9
## 3    0.555 0.4861 0.6239           0.95
```


  
```r
fit.bayes$type = "Bayes"
fit.ml$type = "Max. Lik."
fit.both = rbind(fit.bayes, fit.ml)
fit.both$Interval.width = as.factor(fit.both$Interval.width)
#melt.m = melt(fit.ml, id.vars=c('Interval.width',
#   'lower', 'upper'))
#tmp =melt(tot, id.vars=c('type', 'Interval.width',
#   'variable'))

## old plots
#dat = data.frame(prior=prior, likelihood=likelihood,
#   posterior=naive.posterior, probability.boy=p.b)
#g =
#   ggplot(dat)+geom_vline(xintercept=pb.max,color='red')
#   + geom_vline(xintercept=post.50,color='blue')

g = ggplot(fit.both)
g + geom_pointrange(aes(Interval.width, Estimate, 
    ymin = lower, ymax = upper, color = type), position = position_dodge(width = 0.1))
```
![plot of chunk prob1a_fig](https://github.com/ashander/stat-rethink/raw/master/prob1a_fig.png)

  
# Problem 2

## 2a

 Simulating the number of boys born in 200 births using `rbinom`, 10000 replicates

```r
# define a function to use in simulations
Num.boys <- function(p.b.this, n.total) {
    births.this <- rbinom(n = n.total, size = 1, prob = p.b.this)
    return(sum(births.this))
}

p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = naive.posterior)
boys.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 200)
})
dens(boys.sims, adj = 1)
abline(v = sum(c(birth1, birth2)), lwd = 3, col = "red")
```
![plot of chunk prob2a](https://github.com/ashander/stat-rethink/raw/master/prob2a.png)


  The actual total births is right in the center of the distribution of predicted number of boys.
From this use of posterior predictive simulation, the model seems to fit the data well.


## 2b

 Simulating the number of first born boys in 100 births using `rbinom`, 10000 replicates

```r
p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = naive.posterior)
fb.boys.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 100)
})
dens(fb.boys.sims, adj = 1)
abline(v = sum(c(birth1)), lwd = 3, col = "red")
```
![plot of chunk prob2b](https://github.com/ashander/stat-rethink/raw/master/prob2b.png)


  The actual count of first born boys is less than the center of the predicted distribution.
Here, it looks like the model overpredicts the number of first-born boys.


## 2c

 Simulating the number of second born boys with sisters births using `rbinom`, 10000 replicates.

How many first-born girls were there, and how were the siblings born after girls distributed?
Note that the data are ordered by family.

```
## first born girls:  49 , boys born after girls:  39 
```


  
So, 49 girls were born first, and following those births, 39 boys and 10 girls were born. 
Now simulate 49 births and compare to acutal data (focusing as before on boys):

```r
p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = naive.posterior)
pg.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 49)
})
dens(pg.sims, adj = 1)
abline(v = sum(c(birth.pg)), lwd = 3, col = "red")
```
![plot of chunk prob2c](https://github.com/ashander/stat-rethink/raw/master/prob2c.png)


  The actual number of boys born girls (39 boys following 49 girls) exceeds by  _far_ the center of the predicted distribution.
This could be explained if a boy is more likely to be born following a girl.
Yet, the model performs well on births in aggregate.
This could be explained if the assumption of independence between sex of first and second born is invalid in general, and the increased likelihood of bearing girls following boys cancels out the above-demonstrated underprediction of second born boys following girls.



  
# Colophon 

Written using the excellent [knitr](http://yihui.github.com/knitr/).

Need to use not only the options below suggested by knitr docs

```r
require(knitr)  ## the package
opts_knit$set(base.url = "https://github.com/ashander/stat-rethink/raw/master/")
```



but also 

```r
opts_knit$set(out.format = "gfm")
knit("/Users/jaime/PHD/stat-rethink/hw1ashander_knit_.md")  ## to run
```



