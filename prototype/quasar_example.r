### R script to demonstrate RFF version of time delay GP model

## Define imagined data generating process and draw a mock dataset

Nobs_quasar_image_A <- 150
observation_times_quasar_image_A <- sort(runif(Nobs_quasar_image_A,0,10))
Nobs_quasar_image_B <- 75
observation_times_quasar_image_B <- sort(runif(Nobs_quasar_image_B,0,10))

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
plot(true_latent_curve)

true_latent_curve_quasar_image_A <- true_latent_curve[which(sort.list(c(observation_times_quasar_image_A,observation_times_quasar_image_B+true_time_delay)) %in% 1:(Nobs_quasar_image_A))]
true_latent_curve_quasar_image_B <- true_latent_curve[which(!(sort.list(c(observation_times_quasar_image_A,observation_times_quasar_image_B+true_time_delay)) %in% 1:(Nobs_quasar_image_A)))]
observed_counts_quasar_image_A <- rpois(Nobs_quasar_image_A,true_expected_background_counts+true_expected_counts_quasar_image_A*true_latent_curve_quasar_image_A)
observed_counts_quasar_image_B <- rpois(Nobs_quasar_image_B,true_expected_background_counts+true_expected_counts_quasar_image_B*true_latent_curve_quasar_image_B)

plot(observation_times_quasar_image_A,observed_counts_quasar_image_A,ylim=c(0,50))
points(observation_times_quasar_image_B,observed_counts_quasar_image_B,col="red")

## Construct Bayesian model with RFF approximation of latent Gaussian process and prepare derivatives with PyTorch

library(randtoolbox)
library(reticulate)

k_RFF <- 500 # Number of RFF bases (will use sin&cos orthogonal pairs)
omega_qmc_RFF <- t(qnorm(halton(k_RFF,1)))

torch <- import("torch")
np <- import("numpy")
Variable <-  torch$autograd$Variable

omega_qmc_RFF <- torch$from_numpy(omega_qmc_RFF)$float()
observation_times_quasar_image_A <- torch$from_numpy(np$array(observation_times_quasar_image_A))$float()
observation_times_quasar_image_B <- torch$from_numpy(np$array(observation_times_quasar_image_B))$float()

observed_counts_quasar_image_A <- torch$from_numpy(np$array(observed_counts_quasar_image_A))$float()
observed_counts_quasar_image_B <- torch$from_numpy(np$array(observed_counts_quasar_image_B))$float()

time_delay <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
log_GP_bandwidth <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
log_GP_scale <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
log_expected_counts_quasar_image_A <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
log_expected_counts_quasar_image_B <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
log_expected_background_counts <- Variable(torch$zeros(as.integer(1)),requires_grad=TRUE)$float()
RFF_slopes <- Variable(torch$ones(as.integer(2*k_RFF)),requires_grad=TRUE)$float()

optimizer <- torch$optim$Adam(list(time_delay,log_GP_bandwidth,log_GP_scale,log_expected_counts_quasar_image_A,log_expected_counts_quasar_image_B,log_expected_background_counts,RFF_slopes),lr=0.05)

for (j in 1:200) {
GP_bandwidth <- torch$exp(log_GP_bandwidth)
GP_scale <- torch$exp(log_GP_bandwidth)
expected_background_counts <- torch$exp(log_expected_background_counts)

omega_QMC_RFF_scaled <- torch$mul(omega_qmc_RFF,GP_bandwidth$unsqueeze(as.integer(1)))
RFFprojections_A <- torch$matmul(observation_times_quasar_image_A$unsqueeze(as.integer(1)),omega_QMC_RFF_scaled)
RFFprojections_B <- torch$matmul(observation_times_quasar_image_B$add(time_delay)$unsqueeze(as.integer(1)),omega_QMC_RFF_scaled)
RFFs_A <- torch$cat(c(torch$sin(RFFprojections_A),torch$cos(RFFprojections_A)),as.integer(1))
RFFs_B <- torch$cat(c(torch$sin(RFFprojections_B),torch$cos(RFFprojections_B)),as.integer(1))
latent_GP_A <- torch$matmul(RFFs_A,RFF_slopes)$mul(GP_scale$mul(1.0/k_RFF))
latent_GP_B <- torch$matmul(RFFs_B,RFF_slopes)$mul(GP_scale$mul(1.0/k_RFF))

latent_signal_A <- torch$exp(latent_GP_A$add(log_expected_counts_quasar_image_A))
latent_signal_B <- torch$exp(latent_GP_B$add(log_expected_counts_quasar_image_B))
obs_signal_A <- latent_signal_A$add(expected_background_counts)
obs_signal_B <- latent_signal_B$add(expected_background_counts)

likelihood_fn_A <- torch$distributions$Poisson(obs_signal_A)
likelihood_A <- likelihood_fn_A$log_prob(observed_counts_quasar_image_A)
likelihood_fn_B <- torch$distributions$Poisson(obs_signal_B)
likelihood_B <- likelihood_fn_B$log_prob(observed_counts_quasar_image_B)
prior_slopes_fn <- torch$distributions$Normal(0,1)
prior_slopes <- prior_slopes_fn$log_prob(RFF_slopes)
prior_time_delay_fn <- torch$distributions$Normal(0,1)
prior_time_delay <- prior_time_delay_fn$log_prob(time_delay)
prior_log_GP_bandwidth_fn <- torch$distributions$Normal(0,1)
prior_log_GP_bandwidth <- prior_log_GP_bandwidth_fn$log_prob(log_GP_bandwidth)
prior_log_GP_scale_fn <- torch$distributions$Normal(0,1)
prior_log_GP_scale <- prior_log_GP_scale_fn$log_prob(log_GP_scale)
prior_log_expected_background_counts_fn <- torch$distributions$Normal(0,1)
prior_log_expected_background_counts <- prior_log_expected_background_counts_fn$log_prob(log_expected_background_counts)

posterior_log_prob_part_i <- likelihood_A$sum()
posterior_log_prob_part_ii <- torch$add(posterior_log_prob_part_i,likelihood_B$sum())
posterior_log_prob_part_iii <- torch$add(posterior_log_prob_part_ii,prior_slopes$sum())
posterior_log_prob_part_iv <- torch$add(posterior_log_prob_part_iii,prior_time_delay)
posterior_log_prob_part_v <- torch$add(posterior_log_prob_part_iv,prior_log_GP_bandwidth)
posterior_log_prob_part_vi <- torch$add(posterior_log_prob_part_v,prior_log_expected_background_counts)
posterior_log_prob <- torch$add(posterior_log_prob_part_vi,prior_log_GP_scale)$mul(-1)

optimizer$zero_grad()
if (j==1) {posterior_log_prob$backward(retain_graph=TRUE)}
else {posterior_log_prob$backward()}
optimizer$step()
print(posterior_log_prob,"\n")

}

latent_fn_all_Aobs <- exp(c(as.numeric(latent_GP_A$detach()$numpy()),as.numeric(latent_GP_B$detach()$numpy()))+as.numeric(log_expected_counts_quasar_image_A$detach()$numpy()))+as.numeric(expected_background_counts$detach()$numpy())
relative_obs_times_to_A <- c(as.numeric(observation_times_quasar_image_A$detach()$numpy()),as.numeric(observation_times_quasar_image_B$detach()$numpy())+as.numeric(time_delay$detach()$numpy()))
latent_fn_all_Bobs <- exp(c(as.numeric(latent_GP_A$detach()$numpy()),as.numeric(latent_GP_B$detach()$numpy()))+as.numeric(log_expected_counts_quasar_image_B$detach()$numpy()))+as.numeric(expected_background_counts$detach()$numpy())

layout(1:2)
plot(observation_times_quasar_image_A$detach()$numpy(),observed_counts_quasar_image_A$detach()$numpy(),xlim=c(0,10),ylim=c(0,125),xlab='Observation Times A',ylab="Observed Counts A")
lines(sort(relative_obs_times_to_A),latent_fn_all_Aobs[sort.list(relative_obs_times_to_A)],col="red")

plot(observation_times_quasar_image_B$detach()$numpy(),observed_counts_quasar_image_B$detach()$numpy(),xlim=c(0,10),ylim=c(0,125),xlab='Observation Times B',ylab="Observed Counts B")
lines(sort(relative_obs_times_to_A)-as.numeric(time_delay$detach()$numpy()),latent_fn_all_Bobs[sort.list(relative_obs_times_to_A)],col="red")

cat("True time delay = ",true_time_delay,"\n",sep="")
cat("Estimated time delay = ",as.numeric(time_delay$detach()$numpy()),"\n",sep="")
