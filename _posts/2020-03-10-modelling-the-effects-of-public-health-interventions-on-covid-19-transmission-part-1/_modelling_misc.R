
----

# Other stuff

```
prefecture_plot_lambda <- function(res) {
  plot(res, "lambdas", scale = length(res$local_case_dates) + 1,
       bty="n")
  title(sub=paste("\nEstimated", expression(lambda), "for", 
                            res$prefecture, 
                            "(assuming serial interval mean =",
                            res$si_mean, 
                            ", sd =", 
                            res$si_sd, ")"))
  abline(v = res$local_case_dates, lwd = 3, col = "grey")
  abline(v = res$last_date, col = "blue", lty = 2, lwd = 2)
  points(res$local_case_dates, seq_along(res$local_case_dates), pch = 20, cex = 3)
```

# Initial model 

```{r, eval=FALSE}
# Initialize the network
nw <- network.initialize(1000, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates  + concurrent

# Input the appropriate target statistics for each term
target.stats <- c(37500, 50, 950) # , 950

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 40)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)
```

## Model diagnostics

```{r, eval=FALSE}
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps, 
            dynamic=FALSE, keep.tedgelist = FALSE, 
            nwstats.formula = ~edges + isolates + degree(0:100) + concurrent)
```

```{r, eval=FALSE}
print(dx)
```

```{r, eval=FALSE}
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)
```

## Epidemic model simulation 

```{r}
nsims <- 5
ncores <- 6
nsteps <- 400

act.rates <- c(10, 5, 2)
inf.probs <- c(0.05, 0.02)

# Model parameters
param <- param.net(inf.prob = 0.05, act.rate = 10,
                   ei.rate = 0.125, ir.rate = 0.03)

# Initial conditions
init <- init.net(i.num = 10)

source("_posts/2020-03-03-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/module-fx.R")

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL,
                       verbose.int = 1,
                       verbose=TRUE)
```

## Run the network model simulation with netsim

```{r, eval=FALSE}
sim <- netsim(est, param, init, control)
print(sim)
```

## Plot outcomes

```{r, eval=FALSE}
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("i.num", "e.num", "r.num"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)
```

```{r, eval=FALSE}
plot(sim, y = c("se.flow", "ei.flow", "ir.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 40), legend = TRUE)
```

```{r, eval=FALSE}
# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
kable(df[c(2, 30, 60, 90, 120), ])
```

# Reduced act rate (5 per day instead of 20)

```{r, eval=FALSE}
# Initialize the network
nw <- network.initialize(1000, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates + concurrent

# Input the appropriate target statistics for each term
target.stats <- c(750, 200, 500)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)
```

## Model diagnostics

```{r, eval=FALSE}
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps, 
            dynamic=FALSE, keep.tedgelist = FALSE, sequential = TRUE,
            nwstats.formula = ~edges + isolates + degree(0:10) + concurrent)
```

```{r, eval=FALSE}
print(dx)
```

```{r, eval=FALSE}
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)
```

# Epidemic model simulation 

```{r, eval=FALSE}
# Model parameters
param <- param.net(inf.prob = 0.02, act.rate = 5,
                   ei.rate = 0.125, ir.rate = 0.03)

# Initial conditions
init <- init.net(i.num = 10)

source("_posts/2020-03-03-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/module-fx.R")

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL)
```

## Run the network model simulation with netsim

```{r, eval=FALSE}
sim <- netsim(est, param, init, control)
print(sim)
```

## Plot outcomes

```{r, eval=FALSE}
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("i.num", "e.num", "r.num"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)
```

```{r, eval=FALSE}
plot(sim, y = c("se.flow", "ei.flow", "ir.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 40), legend = TRUE)
```

```{r, eval=FALSE}
# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
kable(df[c(2, 30, 60, 90, 120), ])
```

# Reduced concurrency (50/1000 instead of 500/1000)

```{r, eval=FALSE}
# Initialize the network
nw <- network.initialize(1000, directed = FALSE)

# Define the formation model: edges + isolates (number with degree of 0)
formation = ~edges + isolates + concurrent

# Input the appropriate target statistics for each term
target.stats <- c(750, 200, 50)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)
```

## Model diagnostics

```{r, eval=FALSE}
dx <- netdx(est, nsims = nsims, ncores = ncores, nsteps = nsteps, 
            dynamic=FALSE, keep.tedgelist = FALSE, sequential = TRUE,
            nwstats.formula = ~edges + isolates + degree(0:10) + concurrent)
```

```{r, eval=FALSE}
print(dx)
```

```{r, eval=FALSE}
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)
```

# Epidemic model simulation 

```{r, eval=FALSE}
# Model parameters
param <- param.net(inf.prob = 0.02, act.rate = 20,
                   ei.rate = 0.125, ir.rate = 0.03)

# Initial conditions
init <- init.net(i.num = 10)

source("_posts/2020-03-03-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/module-fx.R")

# Control settings
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL)
```

## Run the network model simulation with netsim

```{r, eval=FALSE}
sim <- netsim(est, param, init, control)
print(sim)
```

## Plot outcomes

```{r, eval=FALSE}
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = c("i.num", "e.num", "r.num"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)
```

```{r, eval=FALSE}
plot(sim, y = c("se.flow", "ei.flow", "ir.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = TRUE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 40), legend = TRUE)
```

```{r, eval=FALSE}
# Average across simulations at beginning, middle, end
df <- as.data.frame(sim)
kable(df[c(2, 30, 60, 90, 120), ])
```
------
    
infected <- 1000
n_rec <- numeric(100)
rec_times <- numeric(0)

for (d in 1:100) {
    n_rec[d] <- rbinom(size=infected, n=1, prob = 1/20)
    infected <- infected - n_rec[d]
    rec_times <- c(rec_times, rep(d, n_rec[d]))
    print(c(d,n_rec[d]))
}

n_rec

weighted.mean(1:100, n_rec)

rec_times

median(rec_times)

plot(n_rec)

cumsum(n_rec)

set.seed(27)

X <- Binomial(1000, 1/14)
X

random(X, 10)

pdf(X, 2L)
log_pdf(X, 2L)

cdf(X, 8L)
quantile(X, 0.5)

cdf(X, quantile(X, 0.7))
quantile(X, cdf(X, 7))
