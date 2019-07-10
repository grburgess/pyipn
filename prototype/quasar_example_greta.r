### R script to demonstrate Random Fourier Features version of latent Gaussian Process model for quasar time delay estimation
#install_tensorflow(extra_packages = "tensorflow-probability") 

library(tensorflow)
#Sys.setenv(TENSORFLOW_PYTHON="/Users/jburgess/.environs/tflow/bin/python3")
Sys.setenv("CUDA_VISIBLE_DEVICES" = -1) # temporary fix to avoid bug in Greta concerning GPU computations
#use_python("/Users/jburgess/.environs/tflow/bin/python3")
use_virtualenv("/Users/jburgess/.environs/tflow")
library(greta)
#Sys.setenv(TENSORFLOW_PYTHON="/Users/jburgess/.environs/tflow/bin/python3")
### Define imagined data generating process and draw a mock dataset
###
### In this toy model there are two lensed images of the quasar observed at different times over a 50 day window
### Each image has time-series of brightness, one being offset by a fixed delay from the other
### The brightness time-series is modelled as an exponentiated Gaussian process modulating the average brightness of each quasar image

Nobs_quasar_image_A <- 150
observation_times_quasar_image_A <- sort(runif(Nobs_quasar_image_A,0,50))
Nobs_quasar_image_B <- 75
observation_times_quasar_image_B <- sort(runif(Nobs_quasar_image_B,0,50))

true_time_delay <- 1
true_GP_bandwidth <- 1
true_GP_scale <- 0.5
true_expected_counts_quasar_image_A <- 40
true_expected_counts_quasar_image_B <- 30
true_expected_background_counts <- 5

true_cov_matrix <- true_GP_scale^2*exp(-true_GP_bandwidth*as.matrix(dist(sort(c(observation_times_quasar_image_A,observation_times_quasar_image_B+true_time_delay)),upper = TRUE,diag = TRUE))^2)
true_cov_matrix <- true_cov_matrix+diag(Nobs_quasar_image_A+Nobs_quasar_image_B)*0.0001

library(MASS)
true_latent_curve <- exp(mvrnorm(1,rep(0,Nobs_quasar_image_A+Nobs_quasar_image_B),true_cov_matrix))

true_latent_curve_quasar_image_A <- true_latent_curve[which(sort.list(c(observation_times_quasar_image_A,observation_times_quasar_image_B+true_time_delay)) %in% 1:(Nobs_quasar_image_A))]
true_latent_curve_quasar_image_B <- true_latent_curve[which(!(sort.list(c(observation_times_quasar_image_A,observation_times_quasar_image_B+true_time_delay)) %in% 1:(Nobs_quasar_image_A)))]
observed_counts_quasar_image_A <- rpois(Nobs_quasar_image_A,true_expected_background_counts+true_expected_counts_quasar_image_A*true_latent_curve_quasar_image_A)
observed_counts_quasar_image_B <- rpois(Nobs_quasar_image_B,true_expected_background_counts+true_expected_counts_quasar_image_B*true_latent_curve_quasar_image_B)

layout(1:2)
plot(observation_times_quasar_image_A,observed_counts_quasar_image_A,ylim=c(0,150))
title("Mock Dataset: Image A")
plot(observation_times_quasar_image_B,observed_counts_quasar_image_B,ylim=c(0,150))
title("Mock Dataset: Image B")

### Construct Bayesian model with RFF approximation of latent Gaussian process and sample with Greta (HMC)
###
### The Bayesian model mimics the mock data generating process with the parameters assigned sensible priors, close to their true values
### The implementation of the Gaussian process is by way of the RFF approximation: https://people.eecs.berkeley.edu/~brecht/papers/07.rah.rec.nips.pdf
### Not implemented here: RFFs are easily extensible to represent non-stationary Gaussian processes: https://www.sciencedirect.com/science/article/pii/S2211675317302890
### The RFF code here uses a low discrepancy sequence (Quasi-Monte Carlo, rather than completely random), and works with sin/cos pairs (rather than random phases), since these have been shown to perform well
###
### The beauty of the greta package is that you can write your model mostly using base R commands and it will built a graph in tensorflow automatically for you to provide fast HMC

library(randtoolbox)

k_RFF <- 500 # Number of RFF bases (will use sin&cos orthogonal pairs)
omega_qmc_RFF <- t(qnorm(halton(k_RFF,1)))

time_delay <- normal(0,1)
log_GP_bandwidth <- normal(0,1)
log_GP_scale <- normal(0,1)
log_expected_counts_quasar_image_A <- normal(log(50),2)
log_expected_counts_quasar_image_B <- normal(log(50),2)
log_expected_background_counts <- normal(log(5),2)
RFF_slopes <- normal(rep(0,k_RFF*2),1)

GP_bandwidth <- exp(log_GP_bandwidth)
GP_scale <- exp(log_GP_bandwidth)
expected_background_counts <- exp(log_expected_background_counts)

omega_QMC_RFF_scaled <- omega_qmc_RFF*GP_bandwidth
RFFprojections_A <- observation_times_quasar_image_A%*%omega_QMC_RFF_scaled
RFFprojections_B <-  (observation_times_quasar_image_B+time_delay)%*%omega_QMC_RFF_scaled
RFFs_A <- cbind(sin(RFFprojections_A),cos(RFFprojections_A))
RFFs_B <- cbind(sin(RFFprojections_B),cos(RFFprojections_B))
latent_GP_A <- (RFFs_A%*%RFF_slopes)*GP_scale/sqrt(k_RFF)
latent_GP_B <- (RFFs_B%*%RFF_slopes)*GP_scale/sqrt(k_RFF)

latent_signal_A <- exp(latent_GP_A+log_expected_counts_quasar_image_A)
latent_signal_B <- exp(latent_GP_B+log_expected_counts_quasar_image_B)
obs_signal_A <- latent_signal_A+expected_background_counts
obs_signal_B <- latent_signal_B+expected_background_counts

distribution(observed_counts_quasar_image_A) <- poisson(obs_signal_A)
distribution(observed_counts_quasar_image_B) <- poisson(obs_signal_B)

m <- model(time_delay,log_GP_bandwidth,log_GP_scale,log_expected_counts_quasar_image_A,log_expected_counts_quasar_image_B,log_expected_background_counts,RFF_slopes)

optim <- opt(m,max_iterations = 100) ## Take a few simple optimising steps to get towards the bulk of posterior mass before starting HMC

draws <- mcmc(m,n_samples = as.integer(1000),chains=as.integer(1)) ## Single chain here: there are obvious benefits to running multiple chains for a proper analysis

posterior.latent.GP_A <- matrix(0,nrow=25,ncol=(Nobs_quasar_image_A+Nobs_quasar_image_B))
posterior.latent.GP_B <- matrix(0,nrow=25,ncol=(Nobs_quasar_image_A+Nobs_quasar_image_B))
posterior.relative_observation_times_to_A <- matrix(0,nrow=25,ncol=(Nobs_quasar_image_A+Nobs_quasar_image_B))
posterior.relative_observation_times_to_B <- matrix(0,nrow=25,ncol=(Nobs_quasar_image_A+Nobs_quasar_image_B))
for (i in 1:25) {
  drawn.par <- list(
    'time_delay'=draws[[1]][i*10,1],
    'log_GP_bandwidth'=draws[[1]][i*10,2],
    'log_GP_scale'=draws[[1]][i*10,3],
    'log_expected_counts_quasar_image_A'=draws[[1]][i*10,4],
    'log_expected_counts_quasar_image_B'=draws[[1]][i*10,5],
    'log_expected_background_counts'=draws[[1]][i*10,6],
    'RFF_slopes'=draws[[1]][i*10,7:(k_RFF*2+7-1)]
  )
  posterior.latent.GP_A[i,] <- as.numeric(calculate(exp(c(latent_GP_A,latent_GP_B)+log_expected_counts_quasar_image_A)+expected_background_counts,values=drawn.par))
  posterior.latent.GP_B[i,] <- as.numeric(calculate(exp(c(latent_GP_A,latent_GP_B)+log_expected_counts_quasar_image_B)+expected_background_counts,values=drawn.par))
  posterior.relative_observation_times_to_A[i,] <-  c(observation_times_quasar_image_A,observation_times_quasar_image_B+drawn.par$time_delay)
  posterior.relative_observation_times_to_B[i,] <-  c(observation_times_quasar_image_A-drawn.par$time_delay,observation_times_quasar_image_B)
  cat(i,"\n")
}

layout(1:2)
plot(observation_times_quasar_image_A,observed_counts_quasar_image_A,xlim=c(0,50),ylim=c(0,125),xlab='Observation Times A',ylab="Observed Counts A")
for (i in 1:25) {lines(sort(relative_obs_times_to_A),posterior.latent.GP_A[i,sort.list(relative_obs_times_to_A)],col=hsv(i/25*0.6,alpha=0.1))}
title("Posterior Draws of Latent Mean Signal: Image A")

plot(observation_times_quasar_image_B,observed_counts_quasar_image_B,xlim=c(0,50),ylim=c(0,125),xlab='Observation Times B',ylab="Observed Counts B")
for (i in 1:25) {lines(sort(relative_obs_times_to_B),posterior.latent.GP_B[i,sort.list(relative_obs_times_to_B)],col=hsv(i/25*0.6,alpha=0.1))}
title("Posterior Draws of Latent Mean Signal: Image B")

layout(1)
hist(draws[[1]][1,],xlab="",ylab="Freq. in Posterior Samples",)
