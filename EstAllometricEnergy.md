Supplement for: Estimating allometric energy allocation between somatic
and gonadal growth
================
Hideyasu Shimadzu and Hui-Yu Wang

This supplement file provides some lines of essential R code used for
the paper.

### Required packages

``` r
library(gam)
library(bbmle)
```

### Data preparation

``` r
data0 <- read.csv("cutlassfish.csv")
data0$soma <- data0$total_w-data0$gonad_w

dt <- 1/(12*10) # the size of time increment: approx 3 days
yr <- seq(0, 6, by=dt) # age cover: 0-6 years
```

### Estimating the probability of maturity, \(\hat{p}_t\)

``` r
#######################
# for site K
locdata.p <- subset(data0, site=="K")
#######################
# for site T
locdata.p <- subset(data0, site=="T")

######################
### estimating maturity probability
#######################
mat <- gam(I(maturity=="M")~lo(newage), data=locdata.p, family="binomial")
p.p <- predict(mat, predict, newdata=data.frame(newage=yr), type="response")
```

This fitting process produces warning messages because of the
concentration of data points at only certain ages. However, this does
not mean the fit has a problem; see Cleveland *et al.* (1991).

Cleveland, W. S., Grosse, E. and Shyu, W. M. (1991). Local Regression
Models in *Statistical Models in S*, Chapman and Hall/CRC.

### Parameter estimation

These are the objective function to be minimised (Section 5.2). The
integration is carried out *via* the trapezium method.

``` r
ssr.v <- function(a, alpha){
  res.v <- pdv - a*pw^alpha
  w <- rep(2, length(res.v))
  w[1] <- w[length(res.v)] <- 1
w%*%(res.v^2)*dt/2
}

ssr.g <- function(b, beta){
  res.w <- pdg - p*b*pw^exp(beta)
  res.w <- res.w[!is.na(res.w)]
  w <- rep(2, length(res.w))
  w[1] <- w[length(res.w)] <- 1
w%*%(res.w^2)*dt/2
}
```

These are the lines of the essential code for parameter estimation.

``` r
#######################
# for site K
hv <- hg <- 0.51 # site"K"
locdata <- subset(data0, site=="K")
#######################
# for site T
hv <- 0.34 # site"T"
hg <- 0.53 # site"T"
locdata <- subset(data0, site=="T")

#######################
### fitting a loess curve
#######################
lo.v <- loess(total_w ~ newage, data=locdata, degree=1, span=hv) # total weight, v
lo.g <- loess(gonad_w ~ newage, data=locdata, degree=1, span=hg) # gonad weight, g

#######################
### predicting a loess curve
#######################
p.v <- predict(lo.v, newdata=data.frame(newage=yr[yr<=4])) # yr<=3 or yr<=4 for site T
p.g <- predict(lo.g, newdata=data.frame(newage=yr[yr<=4])) # yr<=3 or yr<=4 for site T
  x <- p.g
  x[x<0 & !is.na(x)] <- 0
  p.g <- x
p.w <- p.v - p.g # soma weight: w=v-g

######################
### calculating derivatives
#######################
p.dv <- diff(p.v)/dt
p.dw <- diff(p.w)/dt
p.dg <- diff(p.g)/dt

#################
## Estimating a and alpha
#################
pdv <- p.dv
pw <- p.w[-1]
pw <- pw[!is.na(pdv)]
pdv <- pdv[!is.na(pdv)]
pdv <- pdv[pw>0]
pw <- pw[pw>0]

fit.v <- mle2(ssr.v, start=list(a=1, alpha=0.9), lower=c(a=0, alpha=0), method="L-BFGS-B")

#################
## Estimating b and beta
#################
pdg <- p.dg
pw <- p.w[-1]
p <- p.p[-1]
p <- p[!is.na(pdg)]
pw <- pw[!is.na(pdg)]
pdg <- pdg[!is.na(pdg)]
pdg <- pdg[pw>0]
p <- p[pw>0]
pw <- pw[pw>0]

fit.g <- mle2(ssr.g, start=list(b=0.5, beta=0.5), lower=c(b=0, beta=-Inf), method="L-BFGS-B")


a <- coef(fit.v)[1]
alpha <- coef(fit.v)[2]
b <- coef(fit.g)[1]
beta <- exp(coef(fit.g)[2])
```

## Simulation data

### Required package

``` r
library(deSolve)
```

``` r
#######################
### Initial settings
#######################
times <- seq(from=0, to=6, by=1/12)
a <- 4.2
alpha <- 0.6
b <- 0.02
beta <- 1.14
parameters <- c(a, alpha, b, beta)
v <- w <- 7
state <- c(v=v, w=w, g=0)


#######################
### Defining the growth model
#######################
Gmodel <- function(t, state, parms){
  with(as.list(c(state, parms)),{
    p <- 1/(1+exp(6.2-3.2*t))
    dv <- a*w^(alpha)
    dw <- a*w^(alpha) - p*b*w^(beta)
    dg <- p*b*w^(beta)

  list(c(dv, dw, dg), pr=p)
  })
}

Curves <- data.frame(ode(y = state, times = times, func = Gmodel, parms = parameters))
t.pos <- seq(from=0, to=6, by=1/2)
Curves.sub <- Curves[Curves$time%in%t.pos,]


#######################
### Creating the simulation data
#######################
samp <- Curves.sub[sample(1:13, size=250, replace=T),]
samp$ve <- samp$v + with(samp, sapply(time, function(x)rnorm(1, 0, sd=10*sqrt(x))))
samp$ge <- samp$g + with(samp, sapply(time, function(x)rnorm(1, 0, sd=1*sqrt(x))))

p <- 1/(1+exp(6.2-3.2*samp$time))
samp$mat <- sapply(p, function(x)rbinom(1, 1, prob=x))

samp[samp$ve<0, "ve"] <- 0
samp[samp$ge<0|samp$mat==0, "ge"] <- 0

sim.data <- samp
names(sim.data) <- c("newage", "v", "w", "g", "pr", "total_w", "gonad_w", "mat")
```
