#detach("package:vegan", unload=TRUE)
#rm(list = ls())
library(MASS)
library(mvtnorm)
library(gdata)
library(ellipse)

######Example 1, repeat the results for 100 times
alpha = c(1/2,1/2)

mu_1 <- c(0,0)
mu_2 <- c(0,20)

d = length(mu_1)

covar_1 <- matrix(c(1,0,0,1), nrow = 2, byrow = TRUE)
covar_2 <- matrix(c(9,0,0,9), nrow = 2, byrow = TRUE)

n <- 400
n_1  = n*alpha[1]
n_2  = n*alpha[2]


EM_mu <- list()
EM_sigma <- list()

EM_mu_w <- list()
EM_sigma_w <- list()
EM_iteration <- list()

cluster <- list()
mu <- list()
sigma <- list()

mu_w <- list()
sigma_w <- list()
iteration <- list()

loop = 1

for (z in 1:loop){

distr_1 <- mvrnorm(n_1,mu=mu_1,Sigma=covar_1)
distr_2 <- mvrnorm(n_2,mu=mu_2,Sigma=covar_2)

X <- rbind(distr_1,distr_2)

data_list <- list()
for (i in 1:n) data_list[[i]] <- c(X[i,1],X[i,2])

plot(X[,1], X[,2])

####EM clustering algorithm for normal mixture
c = 2
e = 0.00001
s = 1

##inital values
#Step 1
Z_old = matrix(runif(c*n),ncol = c)
Z_old <- t(apply(Z_old, 1, function(i) i/sum(i)))

Z <- Z_old

##Estimation
#Step 2
while(s>0){
  alpha_hat <- (t(c(rep(1,n)))%*%Z)/n
  
  z_row_inv <- diag(c)
  for (i in 1:c) z_row_inv[i,i] <- 1/sum(Z[,i])
  
  mu_hat <- z_row_inv%*%(t(Z)%*%X)
  mu_list <- as.list(as.data.frame(t(mu_hat)))
  
  #Step 3
  covar_hat <- list()
  for (i in 1:c) {
    num_sum = matrix(c(rep(0,d*d)),ncol = d)
    list_1 <- list()
    for (j in 1:n){
      list_1[[j]] <- X[j,] - mu_list[[i]]
      add_term = Z[j,i]*list_1[[j]]%*%t(list_1[[j]])
      num_sum = num_sum + add_term
    }
    covar_hat[[i]] = num_sum/sum(Z[,i])
  }
  
  #Define a function to calculate the posterior distribution with weight alpha
  Posterior <- function(A,v){
    alpha_hat[v]*dmvnorm(A,mean = mu_list[[v]], sigma = covar_hat[[v]])
  }
  
  for (i in 1:c){
    Z[,i] <- apply(X,MARGIN = 1,FUN = Posterior,v=i)
  }
  
  Z <- t(apply(Z, 1, function(i) i/sum(i)))
  
  #Step 4
  if (norm(Z-Z_old, type = "I") < e) {
    t=s
    s = 0}
  else {
    Z_old = Z
    s = s + 1
  }
}

EM_mu[[z]] <- mu_list
EM_sigma[[z]] <- covar_hat

EM_mu_w[[z]] <- c(rep(0,d))
EM_sigma_w[[z]] <- matrix(c(rep(0,d*d)), nrow = d, byrow = TRUE)
EM_iteration[[z]] <- t-1

for (i in 1:c){
  EM_mu_w[[z]] <- EM_mu_w[[z]] + alpha_hat[[i]]*mu_list[[i]]
  EM_sigma_w[[z]] <- EM_sigma_w[[z]] + alpha_hat[[i]]*covar_hat[[i]]
}


####A robust EM clustering algorithm for normal mixture
e = 0.00001
t = 1

##initial values
#Step 1
b = 1
c = n

c_list <- list()
c_list[[1]] <- c

alpha_hat <- c(rep(1/n,n))
mu_hat <- X
mu_list <- as.list(as.data.frame(t(X)))

#Step 2
D <- list()

M <- 0*diag(n)
for (i in 1:d){
  M_diff <- outer(X[,d],X[,d],"-")^2
  M = M + M_diff
}

for (i in 1:n){
  D[[i]] <- sort(M[i,])
}

gamma = 0.0001
Q <- diag(d)

covar_hat <- list()
for (i in 1:c){
  covar_hat[[i]] <- D[[i]][ceiling(sqrt(c))]*diag(d)
}

#Step 3
Z = matrix(runif(c*n),ncol = c)

Posterior <- function(A,v){
  alpha_hat[v]*dmvnorm(A,mean = mu_list[[v]], sigma = covar_hat[[v]])
}

for (i in 1:c){
  Z[,i] <- apply(X,MARGIN = 1,FUN = Posterior,v=i)
}

Z <- t(apply(Z, 1, function(i) i/sum(i)))

##Estimation
#Step 4
#while(t>0){
z_row_inv <- diag(c)
for (i in 1:c) z_row_inv[i,i] <- 1/sum(Z[,i])

mu_hat <- z_row_inv%*%(t(Z)%*%X)
mu_list <- as.list(as.data.frame(t(mu_hat)))

#Step 5
while(t>0){
  alpha_hat_old <- alpha_hat
  
  alpha_EM <- (t(c(rep(1,n)))%*%Z)/n
  
  alpha_old_M <- diag(c)
  for (i in 1:c){
    alpha_old_M[i,i] <- alpha_hat[i]
  }
  
  E <- (t(alpha_hat_old)%*%log(alpha_hat_old))[1]*t(c(rep(1,c)))
  
  alpha_hat <- alpha_EM + b*(log(alpha_hat_old)-E)%*%alpha_old_M
  
  #Step 6
  eta  = min(1,0.5^(floor(d/2-1)))
  
  b1 <- sum(exp(-1*eta*n*abs(alpha_hat-alpha_hat_old)))/c
  b2 <- (1-max(alpha_EM))/(-1*max(alpha_hat_old)*(t(alpha_hat_old)%*%log(alpha_hat_old))[1])
  
  b <- min(b1,b2)
  
  #Step 7
  c <- length(alpha_hat[alpha_hat > 1/n])
  c_list[[t+1]] <- c
  
  Z <- Z[,which(alpha_hat > 1/n)]
  Z <- t(apply(Z, 1, function(i) i/sum(i)))
  
  mu_hat_old <- mu_hat[which(alpha_hat > 1/n),]
  mu_hat <- mu_hat[which(alpha_hat > 1/n),]
  mu_list <- as.list(as.data.frame(t(mu_hat)))
  
  alpha_hat <- alpha_hat[which(alpha_hat > 1/n)]
  alpha_hat <- alpha_hat/sum(alpha_hat)
  
  v <- max(t-60,1)
  
  if ((t >= 61) && (c_list[[v]]-c_list[[t]] == 0)) {b = 0}
  
  #Step 8
  covar_hat <- list()
  for (i in 1:c) {
    num_sum = matrix(c(rep(0,d*d)),ncol = d)
    list_1 <- list()
    for (j in 1:n){
      list_1[[j]] <- X[j,]-mu_list[[i]]
      add_term = Z[j,i]*list_1[[j]]%*%t(list_1[[j]])
      num_sum = num_sum + add_term
    }
    covar_hat[[i]] = (1-gamma)*num_sum/sum(Z[,i])+gamma*Q
  }
  
  #Step 9
  Z_old <- Z
  
  Posterior <- function(A,v){
    alpha_hat[v]*dmvnorm(A,mean = mu_list[[v]], sigma = covar_hat[[v]])
  }
  
  for (i in 1:c){
    Z[,i] <- apply(X,MARGIN = 1,FUN = Posterior,v=i)
  }
  
  Z <- t(apply(Z, 1, function(i) i/sum(i)))
  
  #Step 10
  z_row_inv <- diag(c)
  for (i in 1:c) z_row_inv[i,i] <- 1/sum(Z[,i])
  
  mu_hat <- z_row_inv%*%(t(Z)%*%X)
  mu_list <- as.list(as.data.frame(t(mu_hat)))
  
  #Step 11
  #if (norm(Z-Z_old, type = "I") < e) t = 0
  #else {
  #  Z_old = Z
  #  t = t + 1
  #}
  if (max((mu_hat-mu_hat_old)^2%*%matrix(rep(1,d),ncol=1))<e) {t=0
  } else {t=t+1}
}

cluster[[z]] <- c
mu[[z]] <- mu_list
sigma[[z]] <- covar_hat

mu_w[[z]] <- c(rep(0,d))
sigma_w[[z]] <- matrix(c(rep(0,d*d)), nrow = d, byrow = TRUE)

for (i in 1:c){
  mu_w[[z]] <- mu_w[[z]] + alpha_hat[[i]]*mu_list[[i]]
  sigma_w[[z]] <- sigma_w[[z]] + alpha_hat[i]*covar_hat[[i]]
  iteration[[z]] <- length(c_list)-1
}

}

######Result
EM_mean <- c(rep(0,d))
EM_covariance <- matrix(c(rep(0,d*d)), nrow = d, byrow = TRUE)

for (i in 1:loop){
  EM_mean = EM_mean + EM_mu_w[[i]]/loop
  EM_covariance = EM_covariance + EM_sigma_w[[i]]/loop
}

mean(unlist(EM_iteration))

mu_c <- mu_w[which(unlist(cluster)==2)]
covariance_c <- sigma_w[which(unlist(cluster)==2)]

mean <- c(rep(0,d))
covariance <- matrix(c(rep(0,d*d)), nrow = d, byrow = TRUE)

for (i in 1:length(unlist(cluster)[unlist(cluster)==2])){
  mean = mean + mu_c[[i]]/length(unlist(cluster)[unlist(cluster)==2])
  covariance = covariance + covariance_c[[i]]/length(unlist(cluster)[unlist(cluster)==2])
}

length(unlist(cluster)[unlist(cluster)==2])
mean(unlist(iteration))

####Scatter plot and confidence ellipse (sqrt(5.991*lambda))
#Standard EM
plot(X[,1], X[,2],main = "Standard EM",xlab = "",ylab = "",col="red")
lines(ellipse(x=EM_sigma[[1]][[1]], centre=EM_mu[[1]][[1]],level=0.95)[,1],ellipse(x=EM_sigma[[1]][[1]], centre=EM_mu[[1]][[1]],level=0.95)[,2],type = "l",col="blue")
lines(ellipse(x=EM_sigma[[1]][[2]], centre=EM_mu[[1]][[2]],level=0.95)[,1],ellipse(x=EM_sigma[[1]][[2]], centre=EM_mu[[1]][[2]],level=0.95)[,2],type = "l",col="blue")
points(EM_mu[[1]][[1]][1],EM_mu[[1]][[1]][2],type="p",pch=19)
points(EM_mu[[1]][[2]][1],EM_mu[[1]][[2]][2],type="p",pch=19)

#Robust EM
plot(X[,1], X[,2],main = "Robust EM",xlab = "",ylab = "",col="red")
lines(ellipse(x=sigma[[1]][[1]], centre=mu[[1]][[1]],level=0.95)[,1],ellipse(x=sigma[[1]][[1]], centre=mu[[1]][[1]],level=0.95)[,2],type = "l",col="blue")
lines(ellipse(x=sigma[[1]][[2]], centre=mu[[1]][[2]],level=0.95)[,1],ellipse(x=sigma[[1]][[2]], centre=mu[[1]][[2]],level=0.95)[,2],type = "l",col="blue")
points(mu[[1]][[1]][1],mu[[1]][[1]][2],type="p",pch=19)
points(mu[[1]][[2]][1],mu[[1]][[2]][2],type="p",pch=19)