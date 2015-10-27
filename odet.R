#### odet.R use Cross-validation to select lambda in solving HLASSO

######### HLasso ##########
##### min 1/2 * || y - d*theta - gamma ||_{2,w}^2 + lambda ||gamma||_{1,w} 
##### is equal to
##### (1) min 1/2 * || Uc^T*sqrt(w)*y  - Uc^T*sqrt(w)*gamma ||_2^2 + lambda ||gamma||_{1,w} 
##### (2) min || Ud^T*sqrt(w)*y  - Ud^T*sqrt(w)*\hat{gamma} - Ud * sqrt(w)*d ||_2^2
##### Y=Uc^T*sqrt(w)*y; X=Uc^T*sqrt(w)

##### (1) is used to detect outliers. CV is used here to determine lambda

#### Input
####      d: grad matrix. for kth edge (i,j),
####		d[k,i]=1, d[k,j]=-1
####      y: Comparison result
####      w: weight of each comparison

#### Output
####      success:     1 for success and 0 for fail
####      intercept:   whether intercept is involved
####      theta:       Estimator of score in HLASSO
####      outInd:      Outliers detected in CV, described as index of row in d or y
####      X:           X=U_2^T * sqrt(W)
####      Y:           X=U_2^T * sqrt(W) *y
####      cvout:       Detail of croos-validation result. see "quadrupen" for detail.

odet <- function(d,y=rep(1,dim(d)[1]),w=rep(1,dim(d)[1]),intercept=FALSE,plot=TRUE){
library('quadrupen')

m <- dim(d)[1]
n <- dim(d)[2]
index <- 1:m
eps <- 1e-10
if (intercept){         ### Add intercept to matrix d
	d0=cbind(d,matrix(rep(1,m),ncol=1))
	n=n+1
}

s <- svd(diag(sqrt(w))%*%d,m,n)     ### svd for sqrt(W)*d 
r <- sum(s$d>eps)
if(r!=n-1){
	print("Not connected or problematic d")
	return(list(success=0))
}

### split HLASSO to two optimization problem
U <- s$u[,1:(n-1)]
Sigma <- diag(s$d[1:(n-1)])
V <- s$v[,1:(n-1)]
Uc <- s$u[,n:m]
Ud <- s$u[,1:(n-1)]
Y <- t(Uc)%*%diag(sqrt(w))%*%y
X <-t(Uc)%*%diag(sqrt(w))

### Cross-validation for Outlier Detection, get \gamma
cvout <- crossval(X,Y,lambda2=0,normalize=FALSE,intercept=FALSE,penscale=w)
if (plot){
	plot(cvout)
}

### Plug in \gamma to get the estimator of score theta
A=Sigma%*%t(V)
b=t(Ud)%*%(sqrt(w)*(y-cvout@beta.min))
theta <- lm(b~A-1)
theta=theta$coefficients
theta[is.na(theta)]=0
if (intercept){
	inter=theta[n]
	theta <- theta[1:(n-1)]-mean(theta[1:(n-1)])
}
else{
	inter=0
	theta <- theta-mean(theta)
}

### O is outlier index
O=index[cvout@beta.min!=0]

return(list(success=1,theta=theta,intercept=inter,outInd=O,X=X,Y=Y,cvout=cvout))
}