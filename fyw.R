library('quadrupen')

####### generate simulation data
n=16   ##number of nodes
compare=3000   ## number of comparisons
p=0.1   ## pecentages of reversed comparisons
data=matrix(rep(0,2*compare),ncol=2)   ##  each row (i,j) means a comparison i>j
	##  ground truth is n>n-1>...>1
for (i in 1:compare){
	data[i,1]=floor(runif(1)*16)+1
	data[i,2]=floor(runif(1)*15)+1
	data[i,2]=data[i,2]+(data[i,2]>=data[i,1])
}
data=cbind(apply(data,1,max),apply(data,1,min))

Realout=sample(1:compare)
Realout=Realout[1:round(p*compare)]  ##Index of Reversed comparisons
data[Realout,]=data[Realout,2:1]
m_out=length(Realout)   ## number of reversed comparisons

w=rep(1,compare)   ## weight of each comparison, here 1 for all.

###### merge comparisons on (i,j) with same y_{ij} together to reduce computation 
###### Compute Comparison matrix Z 
###### Z[i,j]=w means there are together w weight of comparisons believe i>j
Z=matrix(rep(0,n*n),ncol=n)
for (k in 1:compare){
    a=as.numeric(data[k,])
    Z[a[1],a[2]]=Z[a[1],a[2]]+w[k]
}
Z

####### Compute d & y & w in merged data
####### y is always 1. w is given by matrix Z
k=0
m=sum(sum(Z!=0))
d=matrix(rep(0,m*n),ncol=n)
w=rep(0,m)
for (i in 1:n){
    for (j in 1:n){
        if (i!=j & Z[i,j]!=0){
            k=k+1
            w[k]=Z[i,j]
            d[k,i]=1
            d[k,j]=-1
        }
    }
}
y=rep(1,m)

######### Least Square Score theta_LS##########
ls_model =lm(y~d-1,weights=w)
theta_LS =ls_model$coefficients
theta_LS[is.na(theta_LS)]=0
theta_LS =as.numeric(theta_LS-mean(theta_LS))

######### HLasso ##########
##### min 1/2 * || y - d*theta - gamma ||_{2,w}^2 + lambda ||gamma||_{1,w} 
##### is equal to
##### (1) min 1/2 * || Uc^T*sqrt(w)*y  - Uc^T*sqrt(w)*gamma ||_2^2 + lambda ||gamma||_{1,w} 
##### (2) min || Ud^T*sqrt(w)*y  - Ud^T*sqrt(w)*\hat{gamma} - Ud * sqrt(w)*d ||_2^2
##### Y=Uc^T*sqrt(w)*y; X=Uc^T*sqrt(w)

m <- dim(d)[1]
n <- dim(d)[2]
index <- 1:m
eps <- 1e-10

s <- svd(diag(sqrt(w))%*%d,m,n)
r <- sum(s$d>eps)
if(r!=n-1){
	print("Not connected or problematic d")
	return(list(success=0))
}

Sigma <- diag(s$d[1:(n-1)])
V <- s$v[,1:(n-1)]
Uc <- s$u[,n:m]
Ud <- s$u[,1:(n-1)]
Y <- t(Uc)%*%diag(sqrt(w))%*%y
X <-t(Uc)%*%diag(sqrt(w))

##### Solve (1); lasso path
lasso=elastic.net(X,Y,lambda2=0,normalize=FALSE,intercept=FALSE,penscale=w)

##### Determine lambda
##Cross-valisation
cvout <- crossval(X,Y,lambda2=0,normalize=FALSE,intercept=FALSE,penscale=w)
plot(cvout)

##### Solve (2)
A=Sigma%*%t(V)
b=t(Ud)%*%(sqrt(w)*(y-cvout@beta.min))
theta_HLasso <- lm(b~A-1)
theta_HLasso=theta_HLasso$coefficients
theta_HLasso[is.na(theta_HLasso)]=0
theta_HLasso=as.numeric(theta_HLasso-mean(theta_HLasso))


############################################################
##### Or simply delete p% of first outliers along lasso path


######## Order of gamma added into lasso path  ##########
######## matrix "a" store the order added. 
######## format of a: lam i j weight  
######## means pair (i,j) with weight w is first added into model when lambda decrease to lam
coefs = lasso@coefficients
k_th = rep(dim(coefs)[1],m)   
for (i in dim(coefs)[1]:1){
	k_th[coefs[i,]!=0]=i
}
temp= sort(k_th,index=TRUE)
k_th=temp$x
order=temp$ix

num=1:n
a=matrix(rep(0,4*m),ncol=4)
for (i in 1:m){
	a[i,2]=num[d[order[i],]==1]
	a[i,3]=num[d[order[i],]==-1]
}
a[,1]=lasso@lambda1[k_th]
a[,4]=w[order]
a

##########  The first p_outlier% outliers.  ###########
p_outlier=0.05

Top_weight=0
Out_weight=p_outlier*compare

for (i in 1:m){
	Top_weight=Top_weight+a[i,4]
	if (Top_weight>Out_weight){
		break
	}
}
m_out=i  #### top m_out rows of "a" are p_outlier% outliers
out_p = a[1:m_out,]
out_p_Ind=order[1:m_out]

########## delete top p% edges and run Least Square, get theta_LassoLS  ##########

ls_model =lm(y[-out_p_Ind]~d[-out_p_Ind,]-1,weights=w[-out_p_Ind])
theta_LassoLS =ls_model$coefficients
theta_LassoLS[is.na(theta_LassoLS)] = 0
theta_LassoLS = as.numeric(theta_LassoLS-mean(theta_LassoLS))

###### score of HLasso when p% edges are regnized as outliers#############

A=Sigma%*%t(V)
b=t(Ud)%*%(sqrt(w)*(y-coefs[k_th[m_out],]))
theta_pLasso <- lm(b~A-1)
theta_pLasso=theta_pLasso$coefficients
theta_pLasso[is.na(theta_pLasso)]=0
theta_pLasso=as.numeric(theta_pLasso-mean(theta_pLasso))

######### p% outlier in lasso path ################
labels <- rep("The rest", m)
labels[out_p_Ind] <- c("First p%")
plot(lasso,labels=labels,xvar="lambda")

############  End  ################
