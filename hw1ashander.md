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
posterior <- prior * likelihood
posterior <- posterior/sum(posterior)
Max.post <- function(parameter, posterior) {
    estimate = parameter[-log(posterior) == min(-log(posterior))]
    estimate
}
pb.max <- Max.post(p.b, posterior)
```



  The maximimum posterior probability occurs with model: `p.b = `

```
## [1] 0.5546
```




## 1 b




To construct highest posterior density intervals, we use the `HPDI`.
Choice of appropriate sample and formatting of output is encapsulated in function `Precis.bayes` (a modified version of `rethinking::precis`). 
  
```r
p.b.sample <- sample(p.b, size = 10000, replace = TRUE, 
    prob = posterior)

CI.data <- t(sapply(CI.types, function(x) {
    unlist(Precis.bayes(p.b, posterior = posterior, param.samples = p.b.sample, 
        level = x))
}))
fit.bayes <- as.data.frame(CI.data)
fit.bayes$Interval.width <- as.character(CI.types)
```



Credible intervals, based on sampling the posterior density, for probabilty of drawing a boy:
  
```
##   Estimate  lower  upper Interval.width
## 1   0.5546 0.5305 0.5776            0.5
## 2   0.5546 0.4965 0.6126            0.9
## 3   0.5546 0.4825 0.6216           0.95
```



  Maximum posterior probability estimate of `p.b` is also reported.

## 1c

For maximum likelihood estimation we use `mle2`, and `rethinking::precis` to construct confidence intervals. 

```r
births <- c(birth1, birth2)
pb.ml <- mle2(y ~ dbinom(size = length(births), 
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
CI.dat.bayes <- t(sapply(CI.types, function(x) {
    unlist(precis(pb.ml, level = x))
}))
fit.ml <- as.data.frame(CI.dat.bayes)
names(fit.ml)[3:4] <- c("lower", "upper")
fit.ml$Interval.width <- as.character(CI.types)
fit.ml$`Std. Error` <- NULL
```



  Confidence intervals, based on quadratic approximation of the likelihood surface, and maximum likelihood estimate of `p.b`, probabilty of drawing a boy:
  
```
##   Estimate  lower  upper Interval.width
## 1    0.555 0.5313 0.5787            0.5
## 2    0.555 0.4972 0.6128            0.9
## 3    0.555 0.4861 0.6239           0.95
```


  
  
As can be seen from the figure below, the quadratic approximation of the density agrees with the Bayesian posterior in the center of the distribution, but deviates wildly in the tails.
Additionally, the ML posterior is symmetric.

```r
dens.compare <- data.frame(probability = p.b, 
    Bayes = -log(posterior), ML = -log(dnorm(p.b, mean = pb.max, 
        sd = sqrt(vcov(pb.ml)))/sum(dnorm(p.b, mean = pb.max, 
        sd = sqrt(vcov(pb.ml))))))
d.c = melt(dens.compare, id.vars = "probability")

g <- ggplot(d.c, aes(probability, value, color = variable))
g + geom_line() + ylab("- log density")
```
![plot of chunk prob1_fig2](https://github.com/ashander/stat-rethink/raw/master/prob1_fig2.png)


The ML confidence intervals and Bayesian HPDIs are shown below recentered by the ML estimate.
The ML intervals are symmetric, while the Bayesian intervals are slightly asymmetric.
Bayesian estimates are slightly below ML estimates, but the intervals tend to agree approximately in width.
Overall, Bayesian and ML estimate perform similarly for this data set. 

  
```r
fit.bayes$type <- "Bayes (HDPI)"
fit.ml$type <- "Max. Lik. (CI)"
estimate.ml <- fit.ml[, 1]  # recenter by the ML estimates
fit.ml[, 1:3] <- fit.ml[, 1:3] - estimate.ml
fit.bayes[, 1:3] <- fit.bayes[, 1:3] - estimate.ml
fit.both <- rbind(fit.bayes, fit.ml)
fit.both$Interval.width <- as.factor(fit.both$Interval.width)

g <- ggplot(fit.both)
g + geom_pointrange(aes(Interval.width, Estimate, 
    ymin = lower, ymax = upper, color = type), position = position_dodge(width = 0.1))
```
![plot of chunk prob1_fig1](https://github.com/ashander/stat-rethink/raw/master/prob1_fig1.png)



  
# Problem 2

## 2a

  Now, we simulate the number of boys born in 200 births using `rbinom`, 10000 replicates:

```r
# define a function to use in simulations
Num.boys <- function(p.b.this, n.total) {
    births.this <- rbinom(n = n.total, size = 1, prob = p.b.this)
    return(sum(births.this))
}

p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = posterior)
boys.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 200)
})
dens(boys.sims, adj = 1)
abline(v = sum(c(birth1, birth2)), lwd = 3, col = "red")
```
![plot of chunk prob2a](https://github.com/ashander/stat-rethink/raw/master/prob2a.png)


Above is the density of the predicted number of boys, with the actual births (red line) right in the center of the distribution of predicted number of boys.
From this use of posterior predictive simulation, the model seems to fit the data well.


## 2b

We now simulate the number of first born boys in 100 births using `rbinom`, 10000 replicates:

```r
p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = posterior)
fb.boys.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 100)
})
dens(fb.boys.sims, adj = 1)
abline(v = sum(c(birth1)), lwd = 3, col = "red")
```
![plot of chunk prob2b](https://github.com/ashander/stat-rethink/raw/master/prob2b.png)


As above, this plot shows the density of predicted number of boys, but this time for first born boys. The actual count of first born boys (red line) is less than the center of the predicted distribution.
Here, it looks like the model overpredicts the number of first-born boys.


## 2c

To simulate the number of second born boys with sisters as first-born, we must ask:
How many first-born girls were there, and how were the siblings born after girls distributed?
Note that the data are ordered by family.

```
## first born girls:  49 , boys born after girls:  39 
```


  
So, 49 girls were born first, and following those births, 39 boys and 10 girls were born. 
Now simulate 49 births and compare to acutal data (focusing as before on boys,  using `rbinom`, 10000 replicates):

```r
p.b.sims <- sample(p.b, size = 10000, replace = TRUE, 
    prob = posterior)
pg.sims <- sapply(p.b.sims, function(x) {
    Num.boys(x, 49)
})
dens(pg.sims, adj = 1)
abline(v = sum(c(birth.pg)), lwd = 3, col = "red")
```
![plot of chunk prob2c](https://github.com/ashander/stat-rethink/raw/master/prob2c.png)


  The actual number of boys born following girls (red line, 39 boys following 49 girls) exceeds by  _far_ the center of the predicted distribution shown in the density plot above.
This could be explained if a boy is more likely to be born following a girl (non-independence).

Yet, the model performs well on births in aggregate (as shown in 2a).
This could be explained if the assumption of independence between sex of first and second born is invalid in general.
Then, if the increased likelihood of bearing girls following boys cancels out the above-demonstrated underprediction of second born boys following girls, the aggregate data could still fit our model, which assumed independence.

# Full code

[Link to the full code.](https://github.com/ashander/stat-rethink/raw/master/hw1ashander_knit_.R)

  
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



