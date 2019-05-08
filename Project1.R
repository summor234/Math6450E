load("C:/Users/hlchengac/Desktop/Math6450E/Project1/GWAS.RData") 

library("numDeriv")
library("scalreg")
library("qut")
library("mvtnorm")

#Due to memory issue, only consider first 200 biomarkers
Y <- data.matrix(L$y[1:500,])
X <- data.matrix(L$X[1:500,1:31000])

rm(L)

names(Y) <- c("Y1", "Y2", "Y3", "Y4")

n <- nrow(X)
q <- ncol(X)

########Linear mixed model

###Standardized genotype matrix
W <- X
for (i in 1:ncol(X)){
  p = sum(W[,i])/(2*nrow(X))
  W[,i] <- (W[,i] - p)/sqrt((2*p*(1-p)*nrow(X)))
}

# W <- data.matrix(W)
# X <- data.matrix(X)
# Y <- data.matrix(Y)

#Create vectors to save the estiamtes
LMM <- list()
LMM_covar <- list() 
her_se <- c(rep(0,ncol(Y)))

Lasso_e <- c(rep(0,ncol(Y)))

RCV <- matrix(c(rep(0,4*50)),ncol=50)

####The fix effect design matrix
##PCA
pca <- prcomp(W, center = TRUE,scale. = TRUE)
pca_x <- data.frame(pca$x)
pca_x <- pca_x[,c(1:10)]
pca_x$X0 <- 1
pca_x <- pca_x[c("X0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
pca_x <- as.matrix(pca_x)

####Implement the MLE algorithm with Brent method for optimization
K = W%*%t(W)
spectral <- eigen(K)
U <- spectral$vectors
S <-diag(nrow(K))

for (j in 1:length(spectral$values)){
 S[j,j] <- spectral$values[j]
}

I <- diag(nrow(U))

for(a in 1:ncol(Y)){
U_X = t(U)%*%pca_x
U_Y = t(U)%*%Y
n = nrow(X)

fr <- function(d,j) {
  -log(det(S+d*I))-n*log(t(U_Y[,j]-U_X%*%solve(t(U_X)%*%solve(S+d*I)%*%U_X)%*%t(U_X)%*%solve(S+d*I)%*%U_Y[,j])%*%solve(S+d*I)
                         %*%(U_Y[,j]-U_X%*%solve(t(U_X)%*%solve(S+d*I)%*%U_X)%*%t(U_X)%*%solve(S+d*I)%*%U_Y[,j])/n)
}

delta <- optim(0, fr,j=a, method = "Brent", lower = exp(-10), upper = exp(10))
delta <- delta$par

b <- solve(t(U_X)%*%solve(S+delta*I)%*%U_X)%*%t(U_X)%*%solve(S+delta*I)%*%U_Y[,a]
var_u = t(U_Y[,a] - U_X%*%b)%*%solve(S+delta*I)%*%(U_Y[,a] - U_X%*%b)
var_e = delta*var_u

theta_0 <- c(b,var_u,var_e)
LMM[[a]] <- theta_0

####Calculate the standard error
#Loglikelihood of Y

# density <- function(theta,y){
#     log(dmvnorm(y, mean = pca_x%*%theta[1:11], sigma = K*theta[12]+theta[13]*I))
#   }

density <- function(theta,y){
  U_Y_T = t(U)%*%y
  -1/2*(log(det(S+theta[13]/theta[12]*I))+n*log(1/n*t(U_Y_T-(U_X)%*%theta[1:11])%*%solve(S+theta[13]/theta[12]*I)%*%(U_Y_T-U_X%*%theta[1:11])))
}

#Calculate the observed fisher information matrix
fisher_obs <- solve(-1*hessian(density,y=Y[,a],x=theta_0))
LMM_covar[[a]] <- fisher_obs

##Standard error of heritability by delta method
gradient <- c(rep(0,11),var_e/(var_e+var_u)^2,-1*var_u/(var_u+var_e)^2)

her_se[a] <- sqrt(as.numeric(t(gradient)%*%fisher_obs%*%gradient/nrow(Y)))

}

LMM
LMM_covar
her_se

########Scaled Lasso
for (a in 1:ncol(Y)){
Lasso_e[a]<-scalreg(pca_x,Y[,a],lam0 = "univ")$hsigma

}

Lasso_e

########Refitted Cross-validation (repeat 50 times)
for (a in 1:ncol(Y)){
for (i in 1:50){
RCV[a,i] <-sigmarcv(Y[,a],pca_x,cv=FALSE)$sigmahat
}
}

sum(RCV[1,]/50)
sum(RCV[2,]/50)
sum(RCV[3,]/50)
sum(RCV[4,]/50)

sort(RCV[1,])
sort(RCV[2,])
sort(RCV[3,])
sort(RCV[4,])
