#' @title A function using to generate the rescaled Epanechnikov kernel sampler
#' @description A function using to generate the rescaled Epanechnikov kernel sampler
#' @param n the number of samples
#' @return a random sample
#' @examples {
#' y=numeric(10000)
#' for(i in 1:10000){
#' y[i]=select(3)
#' }
#' print(y[1:10])
#' }
#' @importFrom stats runif
#' @export
select=function(n){
  U=runif(n,-1,1)
  if(abs(U[3])>=abs(U[2])&&abs(U[3])>=abs(U[1])) y=U[2]
  else y=U[3]
  return(y)
}

#' @title A function using a Monte Carlo simulation by the antithetic variate approach and by the simple Monte Carlo method
#' @description  A function using a Monte Carlo simulation by the antithetic variate approach and by the simple Monte Carlo method
#' @param R the number of samples
#' @param antithetic the choice of the method
#' @return the estimation
#' @examples 
#' \dontrun{
#' set.seed(12345)
#' m=1000
#' MC1=MC2=numeric(m)
#' for (i in 1:m) {
#'  MC1[i] <- MC.theta(R = 1000, anti = FALSE)
#'  MC2[i] <- MC.theta(R = 1000)
#' }
#' }
#' @export
MC.theta <- function(R = 10000, antithetic = TRUE) {
  u <- runif(R/2,0,1)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)
  theta=mean(exp(u))
  return(theta)
}

#' @title A function to compute the skewness
#' @description A function to compute the skewness
#' @param x the samples
#' @return the skewness of the samples
#' @examples {
#' x=c(1,3,6,10)
#' sk(x)
#' }
#' @export
sk=function(x){
  xbar=mean(x)
  m3=mean((x-xbar)^3)
  m2=mean((x-xbar)^2)
  return(m3/m2^1.5)
}

#' @title A function to do Count Five Test
#' @description A function to do Count Five Test
#' @param x samples
#' @param y paired samples
#' @return return 1 (reject) or 0 (do not reject H0)
#' @examples {
#' x=rnorm(100,0,1)
#' y=rnorm(100,0,1)
#' count5test(x,y)
#' }
#' @export
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

#' @title A function to do F Test
#' @description A function to do F Test
#' @param x samples
#' @param y paired samples
#' @param n size of the sample
#' @return return 1 (reject) or 0 (do not reject H0)
#' @examples {
#' x=rnorm(100,0,1)
#' y=rnorm(100,0,1)
#' Ftest(x,y,100)
#' }
#' @importFrom stats var
#' @importFrom stats qf
#' @export
Ftest=function(x,y,n){
  F=var(y)/var(x)
  return(as.integer(F<qf(0.055/2,n-1,n-1)||F>qf(1-0.055/2,n-1,n-1)))
}

#' @title A function to do Paired-t Test
#' @description A function to do Paired-t Test
#' @param x samples
#' @param n size of the sample
#' @return return 1 (reject) or 0 (do not reject H0)
#' @examples {
#' x=rnorm(100,0,1)
#' y=rnorm(100,0,1)
#' Ttest(x-y,100)
#' }
#' @export
#' @importFrom stats sd
#' @importFrom stats qt
Ttest=function(x,n){
  T=abs(sqrt(n)*mean(x)/sd(x))
  return(as.integer(T>qt(1-0.055/2,n-1)))
}

#' @title A function to compute the skewness of a  multivariant sample
#' @description A function to compute the skewness of a  multivariant sample
#' @param x the samples
#' @param sigma the variance of the distribution
#' @param n size of the sample
#' @return the skewness of the samples
#' @export
skm=function(x,sigma,n){
  xbar=matrix(rep(colMeans(x),n),nrow=n,byrow=TRUE)
  sum=sum(((x-xbar)%*%solve(sigma)%*%t(x-xbar))^3)
  return(sum/(n^2))
}

#' @title A function to generate a sample with Multivariate normal distribution
#' @name mvrnorm
#' @description A function to generate a sample with Multivariate normal distribution
#' @examples{
#' sigma1=matrix(c(1,0,0,1),2,2)
#' MASS::mvrnorm(1,rep(0,2),sigma1)
#' }
#' @import MASS
#' @importFrom MASS mvrnorm
#' @useDynLib StatComp20038
NULL

#' @title A function to compute the Covariance
#' @description A function to compute the Covariance
#' @param x the samples
#' @param i the row of the samples
#' @return the covariance of the samples
#' @export
#' @importFrom stats cor
b.cor=function(x,i)cor(x[i,1],x[i,2])

#' @title A Metropolis sampler using R
#' @description A Metropolis sampler using R
#' @param sigma Variances  used for the proposal distribution
#' @param x0 Initial value
#' @param N the number of samples
#' @return a random sample of size \code{n} and the rejection rate
#' @examples {
#' rw=rw.Metropolis(2,25,1000)
#' plot(rw$x,type="l",lwd=2,xlab="sigma=2",ylab="X")
#' }
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @export
rw.Metropolis=function(sigma,x0,N){
  f=function(x) 0.5*exp(-abs(x))
  x=numeric(N)
  x[1]=x0
  u=runif(N)
  k=0
  for(i in 2:N){
    y=rnorm(1,x[i-1],sigma)
    if(u[i]<=f(y)/f(x[i-1])) x[i]=y
    else{
      x[i]=x[i-1]
      k=k+1
    }
  }
  return(list(x=x,k=k))
}

#' @title A function using the Gelman-Rubin Method
#' @description A function using the Gelman-Rubin Method
#' @param psi the statistic psi
#' @return the G-R statistc
#' @importFrom stats var
#' @export
Gelman.Rubin=function(psi){
  psi=as.matrix(psi)
  n=ncol(psi)
  k=nrow(psi)
  psi.means=rowMeans(psi)
  B=n*var(psi.means)
  psi.w=apply(psi,1,"var")
  W=mean(psi.w)
  v.hat=W*(n-1)/n+(B/n)
  r.hat=v.hat/W
  return(r.hat)
}

#' @title A function using to generate the M-H sampler 
#' @description A function using to generate the M-H sampler 
#' @param sigma Variances  used for the proposal distribution
#' @param N the number of samples
#' @param X1 Initial value
#' @return a random sample of size \code{n}
#' @importFrom stats runif rnorm dnorm
#' @export
laplace.chain=function(sigma,N,X1){
  f=function(x) 0.5*exp(-abs(x))
  x=rep(0,N)
  x[1]=X1
  u=runif(N)
  for(i in 2:N){
    xt=x[i-1]
    y=rnorm(1,xt,sigma)
    r1=f(y)*dnorm(xt,y,sigma)
    r2=f(xt)*dnorm(y,xt,sigma)
    r=r1/r2
    if(u[i]<=r)x[i]=y
    else x[i]=xt
  }
  return(x)
}

#' @title A function  counts the maximum number of extreme points of each sample with respect to the range of the other sample
#' @description A function  counts the maximum number of extreme points of each sample with respect to the range of the other sample
#' @param x one sample
#' @param y the other sample
#' @return the maximum number of extreme points of each sample with respect to the range of the other sample
#' @examples {
#' x=1:10
#' y=2:11
#' maxout(x,y)
#' }
#' @export
maxout=function(x,y){
  X=x-mean(x)
  Y=y-mean(y)
  outx=sum(X>max(Y))+sum(X<min(Y))
  outy=sum(Y>max(X))+sum(Y<min(X))
  return(max(c(outx,outy)))
}

#' @title A function to compute boot statistic
#' @description A function to compute boot statistic
#' @param z data
#' @param ix the ixth data
#' @param sizes the sizes of the samples
#' @param k The maximum number of nearest neighbours to compute
#' @return boot statistic
#' @import RANN
#' @importFrom RANN nn2
#' @export
Tn=function(z,ix,sizes,k){
  n1=sizes[1]
  n2=sizes[2]
  n=n1+n2
  if(is.vector(z)) z=data.frame(z,0)
  z=z[ix,]
  NN=nn2(data=z,k=k+1)
  block1=NN$nn.idx[1:n1,-1]
  block2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(block1<n1+0.5)
  i2=sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}

#' @title A function to compute the results of the NN methods
#' @description  A function to compute the results of the NN methods
#' @param z data
#' @param sizes the sizes of the samples
#' @param k The maximum number of nearest neighbours to compute
#' @return the results of the NN methods
#' @importFrom boot boot
#' @export
eqdist.nn=function(z,sizes,k){
  R=999
  boot.obj=boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value=mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

#' @title A function to generate BA models
#' @description A function to generate BA models
#' @param m0 Number of nodes in the initial network
#' @param m The number of connected edges of newly added nodes at each moment
#' @param N sizes of the network
#' @return a BA model network
#' @examples 
#' \dontrun{
#' BA(10,2,10000)
#' }
#' @export
BA=function(m0,m,N){
  d=numeric(N)
  p=numeric(N)
  d[1:m0]=m0-1  #generate a complete graph
  for(i in 1:(N-m0)){
    p=d/sum(d)
    n=sample(1:N,size=m,replace=FALSE,prob=p) #Preferred connection
    d[n]=d[n]+1
    d[m0+i]=m
  }
  return(d)
}

#' @title A function to generate hidden variable network models
#' @description A function to generate hidden variable network models
#' @param N sizes of the network
#' @param beta a param in the connection function 
#' @return a hidden variable network model
#' @examples 
#' \dontrun{
#' HV(10000,0.7)
#' }
#' @importFrom stats rexp runif
#' @export
HV=function(N,beta){
  delta=log(N)
  x=rexp(N,1)
  f=function(m,n){
    p=1/(1+exp(-beta*(m+n-delta))) #Connection function
    return(p)
  }
  d=numeric(N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(runif(1)<f(x[i],x[j])){
        d[i]=d[i]+1
        d[j]=d[j]+1}
      else{d[i]=d[i]
      d[j]=d[j]}
    }
  }
  return(d)
}