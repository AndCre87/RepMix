######################################
# Try to implement model with joint P0 (e.g. from Wishart)
# Prior on number of components (e.g., Poisson)
# Follow the algorithm by Argiento and De Iorio
###

rm(list = ls())
library("MASS")
library("salso")
library("mvnfast")
library("MCMCpack")

#Simulated data
set.seed(123)

#Airquality data
y <- as.matrix(airquality[complete.cases(airquality),c(1,2)])
y <- scale(y) #standardise the data
N <- dim(y)[1]
d <- dim(y)[2]
plot(y, pch = 19)


#Set parameters of the NRwp mixture model
M <- 5 #Initial Number of components (allocated and not)
Lambda_M <- 1
gamma_S <- 1

#Unnormalised weights
S <- rep(1, M)
Tsum <- sum(S)
u <- rgamma(1, N, rate = Tsum)

zeta <- 1 #zeta > 0 

#prior on zeta is a gamma(a,b)
update_zeta <- TRUE
a_zeta <- 1
b_zeta <- 1


#Hyperparameters for variances which are component-specific but not repulsive
Psi <- 3 * diag(d)
nu <- d + 4

#Function to compute Selberg integral (normalising constant of joint distribution)
Z_log <- function(zeta, k){
  ind_k <- c(0:(k-1))
  sum_out <- - (k/2 + zeta * k * (k-1)/4) * log(zeta) + k/2 * log(2*pi) +
    sum( lgamma(1 + zeta/2*(ind_k + 1)) ) +
    - k * lgamma(1 + zeta/2)
  return(sum_out)
}
#Try (k = 1 is the normal distribution)
Z_log(zeta, 1)
1/2*log(2*pi/zeta)

#The parameters theta come from a joint distribution. Initalise to marginals
theta <- matrix(rnorm(M * d, 0, sd = sqrt(1/zeta)), M, d)
Sigma <- array(NA, dim = c(M, d, d))
for(m in 1:M){
  Sigma[m,,] <- rwish(nu, Psi)
}

z <- rep(1, N)
K_N <- length(unique(z)) #Number of clusters
M_na <- M - K_N
n_j <- c(as.numeric(table(z)), rep(0, M_na))

#For BD-MH algorithm
p_b <- 0.5
p_d <- 0.5
birth_accept <- 0
birth_count <- 0
death_accept <- 0
death_count <- 0

#Algorithm details
n_burn1 <- 100
n_burn2 <- 2500
thin <- 2
n_save <- 2500
n_tot <- n_burn1 + n_burn2 + thin * n_save

#Adaptive
S_theta <- 0.1 * diag(d)
s_zeta <- 0.01
zeta_accept <- 0

#Output
K_N_out_Gauss <- rep(NA, n_save)
M_out_Gauss <- rep(NA, n_save)
theta_out_Gauss <- vector("list", length = n_save)
Sigma_out_Gauss <- vector("list", length = n_save)
S_out_Gauss <- vector("list", length = n_save)
z_out_Gauss <- matrix(NA, n_save, N)
u_out_Gauss <- rep(NA, n_save)
zeta_out <- rep(NA, n_save)


#Main Gibbs sampler
set.seed(321)
#Progress bar
pb <- txtProgressBar(max = n_tot)
for(it in 1:n_tot){
  
  if(M > 1){
    #Sample allocations z
    for(i in 1:N){
      z_i <- z[i]
      
      #Probability of allocating the variable Z_i
      prob_i <- log(S)
      for(m in 1:M){
        prob_i[m] <- prob_i[m] + dmvn(y[i,], mu = theta[m,], sigma = Sigma[m,,], log = TRUE)
      }
      prob_i <- exp(prob_i - max(prob_i))
      prob_i <- prob_i / sum(prob_i)
      
      z[i] <- sample.int(M, 1, prob = prob_i, replace = TRUE)
      #Update numerosity
      n_j[z_i] <- n_j[z_i] - 1
      n_j[z[i]] <- n_j[z[i]] + 1
    }
    
    #Re-ordering part based on numerosity (easy to do in R)
    #Important to keep the current value of the non-allocated components even if they are all shifted to the end
    n_j_sorted <- sort(n_j, decreasing = TRUE, index.return = TRUE)$ix
    z <- match(z, n_j_sorted)
    n_j <- n_j[n_j_sorted]
    S <- S[n_j_sorted]
    theta <- theta[n_j_sorted,]
    Sigma <- Sigma[n_j_sorted,,]
    K_N <- sum(n_j > 0)
    M_na <- M - K_N
  }
  
  
  #Sample u
  u <- rgamma(1, N, rate = Tsum)
  
  #sample weights S
  S <- apply(matrix(n_j + gamma_S, ncol = 1), 1, function(x){rgamma(1, x, rate = u+1)})
  Tsum <- sum(S)
  
  
  ## sample latent parameters ##
  
  #Allocated components
  for(j in 1:K_N){
    #Indices in this cluster
    index_j <- c(1:N)[which(z == j)]
    
    #Current value of latent
    theta_j <- theta[j,]
    Sigma_j <- Sigma[j,,]
    
    #Propose a new value of theta
    theta_new <- rmvn(1, mu = theta_j, sigma = S_theta)
    
    # Evaluate log-ratio #
    
    #Prior: independent part and #Proposal: log-normal does not cancel (+1)
    log_ratio_theta_j <- - zeta / 2 * sum(theta_new^2 - theta_j^2)
    #Prior: Repulsive part
    for(j_bis in c(1:M)[-j]){
      log_ratio_theta_j <- log_ratio_theta_j + zeta * sum(log(abs(theta_new - theta[j_bis,])))
      log_ratio_theta_j <- log_ratio_theta_j - zeta * sum(log(abs(theta_j - theta[j_bis,])))
    }
    
    #Likelihood
    if(n_j[j] > 1){
      log_ratio_theta_j <- log_ratio_theta_j + sum(apply(y[index_j,], 1, dmvn, mu = theta_new, sigma = Sigma_j, log = TRUE) - apply(y[index_j,], 1, dmvn, mu = theta_j, sigma = Sigma_j, log = TRUE))
    }else{
      log_ratio_theta_j <- log_ratio_theta_j - 0.5 * ((y[index_j,] - theta_new) %*% solve(Sigma_j) %*% t(y[index_j,] - theta_new) - t(y[index_j,] - theta_j) %*% solve(Sigma_j) %*% (y[index_j,] - theta_j))
    }
    
    accept_theta_j <- 1
    if( is.finite(log_ratio_theta_j) ){
      if(log_ratio_theta_j < 0){
        accept_theta_j <- exp(log_ratio_theta_j)
      }
    }else{
      accept_theta_j <- 0
    }
    
    if( runif(1) < accept_theta_j ){
      if(M > 1){
        theta[j,] <- theta_new
      }else{
        theta <- theta_new
      }
    }
    
    #sample latent Sigma's
    if(n_j[j] > 0){
      if(M > 1){
        index_j <- which(z == j)
        nu_post <- nu + n_j[j]
        Psi_post <- Psi + matrix(rowSums(apply(y[index_j,] - matrix(theta[j,], n_j[j], d, byrow = TRUE), 1, function(x){x %*% t(x)})), d, d)
        Sigma[j,,] <- riwish(nu_post, Psi_post)
      }else{
        j <- 1
        index_j <- which(z == j)
        nu_post <- nu + n_j[j]
        Psi_post <- Psi + matrix(rowSums(apply(y[index_j,] - matrix(theta, n_j[j], d, byrow = TRUE), 1, function(x){x %*% t(x)})), d, d)
        Sigma <- riwish(nu_post, Psi_post)
      }
    }else{
      Sigma[j,,] <- riwish(nu, Psi)
    }
  }
  
  #Non allocated components: use BD-MH to also update the number of empty components
  p_b_now <- p_b
  if(M_na == 0){
    p_b_now <- 1
  }
  aux <- runif(1)
  if(aux < p_b_now){#Birth move
    #New M
    M_na_new <- M_na + 1
    M_new <- K_N + M_na_new
    
    #Propose a new value of theta from a normal distribution
    theta_new <- matrix(rnorm(d, 0, sd = sqrt(1/zeta)), 1, d)
    
    # Evaluate log-ratio #
    log_ratio_birth <- 0
    
    #Prior and Mixing measure
    log_ratio_birth <- log_ratio_birth - gamma_S * log(1 + u) + log(Lambda_M) - log(M)
    log_ratio_birth <- log_ratio_birth + d * (Z_log(zeta, M) - Z_log(zeta, M_new))
    #Prior: Repulsive part
    for(j_bis in c(1:M)){
      log_ratio_birth <- log_ratio_birth + zeta * sum(log(abs(theta_new - theta[j_bis,])))
    }
    
    #Proposal
    log_ratio_birth <- log_ratio_birth + d * Z_log(zeta, 1) #proposal is normal
    log_ratio_birth <- log_ratio_birth + log(p_d) - log(p_b_now) - log(M_na_new)
    
    accept_birth <- 1
    if( is.finite(log_ratio_birth) ){
      if(log_ratio_birth < 0){
        accept_birth <- exp(log_ratio_birth)
      }
    }else{
      accept_birth <- 0
    }
    
    birth_accept <- birth_accept + accept_birth
    birth_count <- birth_count + 1
    
    if( runif(1) < accept_birth ){
      M <- M_new
      M_na <- M_na_new
      theta <- rbind(theta, theta_new)
      
      #Add independent latent variables
      S <- c(S, rgamma(1, gamma_S, rate = 1 + u))
      Tsum <- sum(S)
      n_j <- c(n_j, 0)
      Sigma_aux <- array(NA, dim = c(M, d, d))
      for(m in 1:(M-1)){
        Sigma_aux[m,,] <- Sigma[m,,]
      }
      Sigma_aux[M,,] <- rwish(nu, Psi)
      Sigma <- Sigma_aux
    }
  }else{ #Death move
    #New M
    M_na_new <- M_na - 1
    M_new <- K_N + M_na_new
    
    #Select value of theta to be removed
    j_death <- K_N + sample.int(n = M_na, size = 1)
    theta_death <- theta[j_death,]
    
    # Evaluate log-ratio #
    log_ratio_death <- 0
    
    #Prior and Mixing measure
    log_ratio_death <- log_ratio_death + gamma_S * log(1 + u) - log(Lambda_M) + log(M_new)
    log_ratio_death <- log_ratio_death + d * (Z_log(zeta, M) - Z_log(zeta, M_new))
    #Prior: Repulsive part
    for(j_bis in c(1:M)[-j_death]){
      log_ratio_death <- log_ratio_death - zeta * sum(log(abs(theta_death - theta[j_bis,])))
    }
    
    #Proposal
    log_ratio_death <- log_ratio_death - d * Z_log(zeta, 1) #proposal is normal
    log_ratio_death <- log_ratio_death - log(p_d) + log(p_b_now) + log(M_na)
    
    accept_death <- 1
    if( is.finite(log_ratio_death) ){
      if(log_ratio_death < 0){
        accept_death <- exp(log_ratio_death)
      }
    }else{
      accept_death <- 0
    }
    
    death_accept <- death_accept + accept_death
    death_count <- death_count + 1
    
    if( runif(1) < accept_death ){
      M <- M_new
      M_na <- M_na_new
      theta <- theta[-j_death,]
      
      #Remove latent variables
      S <- S[-j_death]
      Tsum <- sum(S)
      n_j <- n_j[-j_death]
      Sigma <- Sigma[-j_death,,]
    }
  }
  
  ##########
  # Update hyperaprameters of repulsive measure jointly because the prior is given together
  
  if(update_zeta){
    #Propose a new value of zeta from a log-normal
    zeta_new <- zeta * exp(sqrt(s_zeta) * rnorm(1))    
    
    # Evaluate log-ratio #
    log_ratio_zeta <- a_zeta * (log(zeta_new) - log(zeta)) - b_zeta *(zeta_new - zeta)
    
    log_ratio_zeta <- log_ratio_zeta + d * (Z_log(zeta, M) - Z_log(zeta_new, M))
    log_ratio_zeta <- log_ratio_zeta - (zeta_new - zeta) / 2 * sum(theta^2)
    
    if(M > 1){
      #Prior: Repulsive part
      for(j1 in c(1:(M-1))){
        for(j2 in c((j1 + 1):M)){
          log_ratio_zeta <- log_ratio_zeta + (zeta_new - zeta) * sum(log(abs(theta[j1,] - theta[j2,])))
        }
      }
    }
    
    accept_zeta <- 1
    if( is.finite(log_ratio_zeta) ){
      if(log_ratio_zeta < 0){
        accept_zeta <- exp(log_ratio_zeta)
      }
    }else{
      accept_zeta <- 0
    }
    
    zeta_accept <- zeta_accept + accept_zeta
    
    if( runif(1) < accept_zeta ){
      zeta <- zeta_new
    }
    
    if(it > n_burn1){
      s_zeta <- s_zeta + it ^(-0.7) * (accept_zeta - 0.234)
      if(s_zeta > exp(50)){
        s_zeta = exp(50)
      }else{
        if(s_zeta < exp(-50)){
          s_zeta = exp(-50)
        }
      }
    }
  }
  
  
  
  
  if(it%%10 == 0){
    print(paste("it = ", it, sep = ""))
    print(paste("M = ", M, sep = ""))
    print(paste("K_N = ", K_N, sep = ""))
    print(paste("M_na = ", M_na, sep = ""))
    print(paste("n_j = ", n_j, sep = ""))
    
    # print(
    # plot(y, pch = 19, col = z)
    # )
  }
  
  #save output
  if((it > (n_burn1 + n_burn2)) & ((it - (n_burn1 + n_burn2)) %% thin == 0)){
    
    iter <- (it - (n_burn1 + n_burn2)) / thin
    
    K_N_out_Gauss[iter] <- K_N
    M_out_Gauss[iter] <- M
    u_out_Gauss[iter] <- u
    theta_out_Gauss[[iter]] <- theta
    Sigma_out_Gauss[[iter]] <- Sigma
    S_out_Gauss[[iter]] <- S
    z_out_Gauss[iter,] <- z
    zeta_out[iter] <- zeta
  }
  
  setTxtProgressBar(pb, it)
}

#Save output
save.image(file = "OUTPUT_MCMC_Gaussian.RData")

print(birth_accept / birth_count)
print(death_accept / death_count)


#Histogram of zeta
pdf("zeta_hist_Gaussian.pdf")
par(mar = c(5,5,5,5))
hist(zeta_out, col = "lightblue", xlab = bquote(zeta), main = "", ylab = bquote(P(zeta~"|"~y)), freq = FALSE, cex.axis = 2, cex.lab = 2)
dev.off()



#Number of components and clusters
pdf("M_K_N_hist_Gaussian.pdf")
plot(table(M_out_Gauss) / n_save, col = "blue", lwd = 4, xlab = "", ylab = "", main = "Number of components and clusters", 
     ylim = c(0,1), xlim = c(0,max(M_out_Gauss)+0.5), axes = FALSE)
axis(2)
axis(1, labels = c(1:max(M_out_Gauss)), at = c(1:max(M_out_Gauss)))
segments(x0 = as.numeric(names(table(K_N_out_Gauss))) + 0.1, y0 = 0, y1 = table(K_N_out_Gauss) / n_save, col = "red", lwd = 4)
legend("topright", legend = c("Number of components", "Number of clusters"), col = c("blue", "red"), lwd = 4, lty = c(1, 1), bty = "n")
dev.off()


#Binder partition
Binder_partition_Gaussian <- c(dlso(z_out_Gauss, loss = "binder"))

K_N_Binder <- length(unique(Binder_partition_Gaussian))
nj_Binder <- table(Binder_partition_Gaussian)
nj_Binder


#Predictive distribution
nx <- 100
ny <- 100
xx <- seq(min(y[,1])-1, max(y[,1])+1, length = nx)
yy <- seq(min(y[,2])-1, max(y[,2])+1, length = ny)

f_d <- matrix(0, nx, ny)
for(it in 1:n_save){
  S_norm <- S_out_Gauss[[it]]
  S_norm <- S_norm / sum(S_norm)
  for(j in 1:M_out_Gauss[it]){
    if(M_out_Gauss[it] > 1){
      mean_it <- theta_out_Gauss[[it]][j,]
      cov_it <- Sigma_out_Gauss[[it]][j,,]
    }else{
      mean_it <- theta_out_Gauss[[it]]
      cov_it <- Sigma_out_Gauss[[it]]
    }
    
    for(ix in 1:nx){
      f_d[ix,] <- f_d[ix,] + S_norm[j] * dmvn(cbind(rep(xx[ix], ny), yy), mu = mean_it, sigma = cov_it)
    }
  }
}
f_d <- f_d / n_save


pdf("AirQuality_contour_Gaussian.pdf")
par(mar = c(5,5,5,5))
contour(xx, yy, f_d, main = "", xlab = colnames(y)[1], ylab = colnames(y)[2], cex.main = 2, cex = 1.25, cex.lab = 2)
points(y, pch = 19, col = Binder_partition_Gaussian, cex = 1.25)
dev.off()


pdf("zeta_post_Gaussian_kde.pdf")
zeta_out_Gauss <- zeta_out
par(mar = c(5,5,5,5))
plot(density(zeta_out_Gauss), col = "#FF6666", cex.main = 2, xlab = bquote(zeta), main = "", ylab = bquote(P(zeta~"|"~y)), cex.axis = 2, cex.lab = 2, lwd = 3)
dev.off()


pdf("zeta_post_Gaussian.pdf")
zeta_out_Gauss <- zeta_out
par(mar = c(5,5,5,5))
hist(zeta_out_Gauss, col = "lightblue", cex.main = 2, xlab = bquote(zeta), main = "", ylab = bquote(P(zeta~"|"~y)), cex.axis = 2, cex.lab = 2, lwd = 3)
dev.off()

quantile(zeta_out_Gauss, probs = c(0.025, 0.975))
mean(zeta_out_Gauss)
median(zeta_out_Gauss)

