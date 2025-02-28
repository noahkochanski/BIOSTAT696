library(tidyverse)
library(spBayes)

set.seed(2025)

birth_imp <- readRDS('birth_imp.RDS')
cluster_locations <- read_csv('cluster_locations.csv') %>% rename(clusterid = DHSCLUST) %>% select(clusterid, LATNUM, LONGNUM)

reg_data <- left_join(birth_imp, cluster_locations)

linear_reg_notspatial <- lm(birth_weight ~ birth_weight_type + mother_bmi + sex_of_child + mother_current_age, data = reg_data)
summary(linear_reg_notspatial)

res <- residuals((linear_reg_notspatial))

residual_data <- reg_data %>% 
  select(LATNUM, LONGNUM) %>%
  mutate(resid = res, 
         across(c(LATNUM, LONGNUM), round, digits = 3))
         
residual_data %>%
  ggplot(aes(x = LONGNUM, y = LATNUM, color = resid)) + 
  geom_point(position = "jitter")+ 
  scale_color_viridis_c()



######

X <- reg_data %>% select(birth_weight_type, mother_bmi, sex_of_child, mother_current_age) %>% as.matrix()
y <- reg_data %>% select(birth_weight) %>% as.matrix()
p <- ncol(X)
n <- nrow(X)
coords <- cluster_locations %>% select(LATNUM, LONGNUM) 
coords <- coords + rnorm(n = p*n)


n.samples <- 2000 # length of MCMC chain


starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(0.1, 3/0.1), 
                 "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.1 <- spLM(y~X, coords=as.matrix(coords), starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

###
load("gc2_results.RData")

library(coda)
theta.samples <- as.mcmc(m.1$p.theta.samples)

traceplot(theta.samples)  # Check if chains are mixing well

densplot(theta.samples)   # Inspect posterior distributions

autocorr.plot(theta.samples)

effectiveSize(theta.samples)

geweke.diag(theta.samples)
