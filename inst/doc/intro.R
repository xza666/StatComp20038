## -----------------------------------------------------------------------------
x=rnorm(50)
y=rnorm(50)
opar <- par()
par(bg="lightyellow", col.axis="blue", mar=c(4, 4, 2.5, 0.25))
plot(x, y, xlab="Fifty random values", ylab="Ten other values",
xlim=c(-2, 2), ylim=c(-2, 2), pch=22, col="red", bg="yellow",
bty="l", tcl=-.25, las=1, cex=1.5)
title("How to customize a plot with R (bis)", font.main=3, adj=1)

## -----------------------------------------------------------------------------
knitr::kable(head(airquality),pad=0)

## -----------------------------------------------------------------------------
knitr::kable(head(airquality),format="html",pad=0)

## -----------------------------------------------------------------------------
set.seed(12345)
n=1000
a=b=2
u=runif(n)
x=b/((1-u)^(1/a))
hist(x,prob=TRUE,breaks=30,main=expression(Pareto(2,2)))
y=seq(2,35,0.1)
lines(y,8*y^(-3))

## -----------------------------------------------------------------------------
select=function(n){
  U=runif(3,-1,1)
  if(abs(U[3])>=abs(U[2])&&abs(U[3])>=abs(U[1])) y=U[2]
  else y=U[3]
  return(y)
}

## -----------------------------------------------------------------------------
y=numeric(10000)
for (i in 1:10000){
  y[i]=select(3)
  i=i+1
}

## -----------------------------------------------------------------------------
hist(y,prob=TRUE,main=expression(f(x)==0.75*(1-x^2)))
a=seq(-1,1,0.01)
lines(a,0.75*(1-a^2))

## -----------------------------------------------------------------------------
set.seed(12345)
n=1000
u=runif(n)
r=4
b=2
x=b/((1-u)^(1/r))-b
hist(x,prob=TRUE,breaks=30,main="pareto cdf with r=4 and beta=2")
y=seq(0,6,0.02)
lines(y,64*(2+y)^(-5))

## -----------------------------------------------------------------------------
set.seed(12345)
m=1e5
x=runif(m,min=0,max=pi/3)
theta.hat=mean(sin(x))*pi/3
print(c(theta.hat,1-cos(pi/3)))

## -----------------------------------------------------------------------------
MC.theta <- function(R = 10000, antithetic = TRUE) {
  u <- runif(R/2,0,1)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)
  theta=mean(exp(u))
  return(theta)
}

## -----------------------------------------------------------------------------
set.seed(12345)
m=1000
MC1=MC2=numeric(m)
for (i in 1:m) {
  MC1[i] <- MC.theta(R = 1000, anti = FALSE)
  MC2[i] <- MC.theta(R = 1000)
}
c(var(MC1),var(MC2),(var(MC1)-var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
c(mean(MC1),mean(MC2),exp(1)-1)

## -----------------------------------------------------------------------------
set.seed(500)
x=seq(1,10,0.01)
g=x^2*exp(-x^2/2)/(sqrt(2*pi))
plot(x,g,ylim=c(0,0.5))
f1=exp(-(x-1.5)^2/2)/(sqrt(2*pi))
lines(x, f1, lty = 3, lwd = 2,col=3)
f2=x*exp(-x/2)/(gamma(2)*4)
lines(x, f2, lty = 3, lwd = 2,col=4)

## -----------------------------------------------------------------------------
set.seed(12345)
m=10000
theta.hat=se=numeric(3)
g=function(x) {
   exp(-1/(2*x^2))/(sqrt(2*pi)*x^4)*(x>0)*(x<1)
}
h=function(x){
  x^2*exp(-x^2/2)/(sqrt(2*pi))*(x>1)
}
x=runif(m)
theta.hat[1]=mean(g(x))
se[1]=sd(g(x))
x=rnorm(m,1.5,1)
fg=h(x)/dnorm(x,1.5,1)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
x <- rgamma(m,2,0.5)
fg <- h(x) / dgamma(x,2,0.5)
theta.hat[3] <- mean(fg)
se[3] <- sd(fg)
res <- rbind(theta_hat=round(theta.hat,5), se=round(se,5))
colnames(res) <- paste0('f',0:2)
print(res)

## -----------------------------------------------------------------------------
M=10000
N=50
T=numeric(5)
est=numeric(N)
g=function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
}
for(i in 1:N){
  for(j in 1:5){
    u=runif(M)
    x=-log(exp(-(j-1)/5)-u*(exp(-(j-1)/5)-exp(-j/5)))
    fg=g(x)/(exp(-x)/(exp(-(j-1)/5)-exp(-j/5)))
    T[j]=mean(fg)
  }
  est[i]=sum(T)
}
round(c(mean(est),sd(est)),5) 

## -----------------------------------------------------------------------------
set.seed(12345)
n=10000
m=1000
alpha=0.05
UCL=matrix(0,m,2)
level=numeric(m)
for(i in 1:m){
  x=rlnorm(n,0,2)
  x=log(x)
  UCL[i,1]=mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  UCL[i,2]=mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  if(UCL[i,1]<0&&UCL[i,2]>0)level[i]=1
  else level[i]=0
}
print(mean(level))

## -----------------------------------------------------------------------------
set.seed(12345)
n=20
alpha=0.05
m=1000
UCL=matrix(0,m,2)
level=numeric(m)
for(i in 1:m){
  x=rchisq(n,df=2)
  UCL[i,1]=mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  UCL[i,2]=mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  if(UCL[i,1]<2&&UCL[i,2]>2)level[i]=1
  else level[i]=0
}
print(mean(level))

## -----------------------------------------------------------------------------
set.seed(12345)
n=30
m=10000
alpha=seq(0.1,10,0.1)
M=length(alpha)
power=numeric(M)
cv=qnorm(0.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sk=function(x){
  xbar=mean(x)
  m3=mean((x-xbar)^3)
  m2=mean((x-xbar)^2)
  return(m3/m2^1.5)
}
for(i in 1:M){
  sktests=numeric(m)
  for(j in 1:m){
    x=rbeta(n,alpha[i],alpha[i])
    sktests[j]=as.integer(abs(sk(x))>=cv)
  }
  power[i]=mean(sktests)
}
plot(alpha,power,type="b",xlab=bquote(alpha))
lines(alpha,power)

## -----------------------------------------------------------------------------
set.seed(12345)
n=30
m=10000
v=seq(1,100,1)
M=length(v)
power=numeric(M)
cv=qnorm(0.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sk=function(x){
  xbar=mean(x)
  m3=mean((x-xbar)^3)
  m2=mean((x-xbar)^2)
  return(m3/m2^1.5)
}
sktests=numeric(m)
for(i in 1:M){
  sktests=numeric(m)
  for(j in 1:m){
    x=rt(n,v[i])
    sktests[j]=as.integer(abs(sk(x))>=cv)
  }
  power[i]=mean(sktests)
}
plot(v,power)
lines(v,power,lty=3)

## -----------------------------------------------------------------------------
#Count Five Test
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
#F Test
Ftest=function(x,y,n){
  F=var(y)/var(x)
  return(as.integer(F<qf(0.055/2,n-1,n-1)||F>qf(1-0.055/2,n-1,n-1)))
}
#Paired-t Test
Ttest=function(x,n){
  T=abs(sqrt(n)*mean(x)/sd(x))
  return(as.integer(T>qt(1-0.055/2,n-1)))
}

## -----------------------------------------------------------------------------
sigma1=1
sigma2=1.5
m=10000
#Small sizes
power1=power2=numeric(m)
for(i in 1:m){
  x=rnorm(20,0,sigma1)
  y=rnorm(20,0,sigma2)
  power1[i]=count5test(x,y)
  power2[i]=Ftest(x,y,20)
}
print(c(mean(power1),mean(power2)))
power=power1-power2
print(Ttest(power,m))

## -----------------------------------------------------------------------------
#Medium sizes
power1=power2=numeric(m)
for(i in 1:m){
  x=rnorm(200,0,sigma1)
  y=rnorm(200,0,sigma2)
  power1[i]=count5test(x,y)
  power2[i]=Ftest(x,y,200)
}
print(c(mean(power1),mean(power2)))
power=power1-power2
print(Ttest(power,m))

## -----------------------------------------------------------------------------
#Large sizes
power1=power2=numeric(m)
for(i in 1:m){
  x=rnorm(2000,0,sigma1)
  y=rnorm(2000,0,sigma2)
  power1[i]=count5test(x,y)
  power2[i]=Ftest(x,y,2000)
}
print(c(mean(power1),mean(power2)))
power=power1-power2
print(Ttest(power,m))

## -----------------------------------------------------------------------------
skm=function(x,sigma,n){
  xbar=matrix(rep(colMeans(x),n),nrow=n,byrow=TRUE)
  sum=sum(((x-xbar)%*%solve(sigma)%*%t(x-xbar))^3)
  return(sum/(n^2))
}
library(MASS)
m=1000
n=c(10,20,30,50,100,500)
cv=qchisq(0.95,4)
p.reject=numeric(length(n))
for(i in 1:length(n)){
  sktests=numeric(m)
  for(j in 1:m){
    sigma=matrix(c(1,0,0,1),2,2)
    x=mvrnorm(n[i],rep(0,2),sigma)
    S=matrix(c(var(x[,1]),0,0,var(x[,2])),2,2)
    sigma.hat=(n[i]-1)*S/n[i]
    b=skm(x,sigma.hat,n[i])
    sktests[j]=as.integer(n[i]*b/6>cv)
  }
  p.reject[i]=mean(sktests)  
}
print(p.reject)

## -----------------------------------------------------------------------------
n=20
m=1000
sigma1=matrix(c(1,0,0,1),2,2)
sigma2=matrix(c(100,0,0,100),2,2)
epsilon=c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N=length(epsilon)
pwr=numeric(N)
cv=qchisq(0.95,4)
for(j in 1:N){
  e=epsilon[j]
  sktests=numeric(m)
  for(i in 1:m){
    x=matrix(rep(0,40),nrow=20,ncol=2)
    for(k in 1:n){
    if(runif(1)<=1-e) x[k,]=mvrnorm(1,rep(0,2),sigma1)
    else x[k,]=mvrnorm(1,rep(0,2),sigma2)}
    S=matrix(c(var(x[,1]),0,0,var(x[,2])),2,2)
    sigma.hat=(n-1)*S/n
    b=skm(x,sigma.hat,n)
    sktests[i]=as.integer(n*b/6>cv)
  }
  pwr[j]=mean(sktests)
}
plot(epsilon,pwr)
lines(epsilon,pwr,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
b.cor=function(x,i)cor(x[i,1],x[i,2])
n=nrow(law)
theta.hat=b.cor(law,1:n)
theta.jack=numeric(n)
for(i in 1:n){
  theta.jack[i]=b.cor(law,(1:n)[-i])
}
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
se.jack=sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),7)

## -----------------------------------------------------------------------------
library(boot)
theta.boot=function(dat,ind){
  x=dat[ind,1]
  mean(x)
}
x=as.matrix(aircondit$hours)
boot.obj=boot(x,statistic = theta.boot,R=2000)
print(boot.ci(boot.out = boot.obj,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
library(boot)
n=nrow(scor)
lambda.hat=eigen(cov(scor))$values
theta.hat=lambda.hat[1]/sum(lambda.hat)
theta.jack=numeric(n)
for(i in 1:n){
  x=scor[-i,]
  lambda=eigen(cov(x))$values
  theta.jack[i]=lambda[1]/sum(lambda)
}
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
se.jack=sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),7)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n=length(magnetic)
e1=e2=e3=e4=c(NULL)
for(i in 1:(n-1)){
  for(j in (i+1):n){
    y=magnetic[c(-i,-j)]
    x=chemical[c(-i,-j)]
    
    J1=lm(y~x)
    yhat1=J1$coef[1]+J1$coef[2]*chemical[i]
    yhat2=J1$coef[1]+J1$coef[2]*chemical[j]
    e1=c(e1,magnetic[i]-yhat1,magnetic[j]-yhat2)
    
    J2=lm(y~x+I(x^2))
    yhat3= J2$coef[1] + J2$coef[2] * chemical[i] + J2$coef[3] * chemical[i]^2
    yhat4= J2$coef[1] + J2$coef[2] * chemical[j] + J2$coef[3] * chemical[j]^2
    e2=c(e2,magnetic[i]-yhat3,magnetic[j]-yhat4)
    
    J3 <- lm(log(y) ~ x)
    logyhat5 <- J3$coef[1] + J3$coef[2] * chemical[i]
    yhat5 <- exp(logyhat5)
    logyhat6 <- J3$coef[1] + J3$coef[2] * chemical[j]
    yhat6 <- exp(logyhat6)
    e3=c(e3,magnetic[i]-yhat5,magnetic[j]-yhat6)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat7 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
    yhat7 <- exp(logyhat7)
    logyhat8 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
    yhat8 <- exp(logyhat8)
    e4=c(e4,magnetic[i]-yhat7,magnetic[j]-yhat8)
  }
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
maxout=function(x,y){
  X=x-mean(x)
  Y=y-mean(y)
  outx=sum(X>max(Y))+sum(X<min(Y))
  outy=sum(Y>max(X))+sum(Y<min(X))
  return(max(c(outx,outy)))
}
n1=20
n2=30
mu1=mu2=0
sigma1=sigma2=1
m=1000
p=numeric(m)
for(i in 1:m){
  x=rnorm(n1,mu1,sigma1)
  y=rnorm(n2,mu2,sigma2)
  z=c(x,y)
  R=999
  reps=numeric(R)
  t0=maxout(x,y)
  for(j in 1:R){
    k=sample(1:(n1+n2),size=n1,replace=FALSE)
    x1=z[k]
    y1=z[-k]
    reps[j]=maxout(x1,y1)
  }
  p[i]=mean(c(t0,reps)>=t0)
}
print(mean(p))

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
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
library(energy)
library(Ball)

## -----------------------------------------------------------------------------
m=1000
k=3
p=2
set.seed(12345)
n1=n2=50
R=999
n=n1+n2
N=c(n1,n2)
eqdist.nn=function(z,sizes,k){
  boot.obj=boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
  ts=c(boot.obj$t0,boot.obj$t)
  p.value=mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*p,0,1),ncol=p)
  y=matrix(rnorm(n2*p,0,1.5),ncol=p)
  z=rbind(x,y)
  p.values[i,1]=eqdist.nn(z,N,k)$p.value
  p.values[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.va
}
alpha=0.1
pow=colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
p.values=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*p,0,1.5),ncol=p)
  y=matrix(rnorm(n2*p,0.3,1),ncol=p)
  z=rbind(x,y)
  p.values[i,1]=eqdist.nn(z,N,k)$p.value
  p.values[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.va
}
alpha=0.1
pow=colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
p.values=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rt(n1*p,df=1),ncol=p)
  y=cbind(rnorm(n2,0,1),rnorm(n2,0.3,1))
  z=rbind(x,y)
  p.values[i,1]=eqdist.nn(z,N,k)$p.value
  p.values[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.va
}
alpha=0.1
pow=colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
m=1000
k=3
p=2
set.seed(12345)
n1=10
n2=100
R=999
n=n1+n2
N=c(n1,n2)
p.values=matrix(NA,m,3)
for(i in 1:m){
  x=matrix(rnorm(n1*p,0,1.5),ncol=p)
  y=matrix(rnorm(n2*p,0,1),ncol=p)
  z=rbind(x,y)
  p.values[i,1]=eqdist.nn(z,N,k)$p.value
  p.values[i,2]=eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3]=bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.va
}
alpha=0.1
pow=colMeans(p.values<alpha)
print(pow)

## -----------------------------------------------------------------------------
set.seed(12345)
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

## -----------------------------------------------------------------------------
N=2000
sigma=c(0.05,0.5,2,16)
x0=25
rw1=rw.Metropolis(sigma[1],x0,N)
rw2=rw.Metropolis(sigma[2],x0,N)
rw3=rw.Metropolis(sigma[3],x0,N)
rw4=rw.Metropolis(sigma[4],x0,N)
plot(rw1$x,type="l",lwd=2,xlab="sigma=0.05",ylab="X")
plot(rw2$x,type="l",lwd=2,xlab="sigma=0.5",ylab="X")
plot(rw3$x,type="l",lwd=2,xlab="sigma=2",ylab="X")
plot(rw4$x,type="l",lwd=2,xlab="sigma=16",ylab="X")

## -----------------------------------------------------------------------------
print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))

## -----------------------------------------------------------------------------
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
sigma=0.8
k=4
n=15000
b=1000
x0=c(-10,-5,5,10)
X=matrix(0,nrow=k,ncol=n)
for(i in 1:k){
  X[i,]=laplace.chain(sigma,n,x0[i])
}
psi=t(apply(X,1,cumsum))
for(i in 1:nrow(psi)){
  psi[i,]=psi[i,]/(1:ncol(psi))
}

## -----------------------------------------------------------------------------
for(i in 1:k){
  plot(psi[i,(b+1):n],type="l",xlab=i,ylab=bquote(psi))
}

## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
rhat=rep(0,n)
for(j in (b+1):n){
  rhat[j]=Gelman.Rubin(psi[,1:j])
}
plot(rhat[(b+1):n],type="l",xlab="",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
k=c(4:25,100,500,1000)
n=length(k)
res=numeric(n)
for(i in 1:n){
  f=function(x) pt(sqrt(x^2*(k[i]-1)/(k[i]-x^2)),df=k[i]-1)-pt(sqrt(x^2*k[i]/(k[i]+1-x^2)),df=k[i])
  res0=uniroot(f,c(0.01,min(sqrt(k[i])-0.01,4)))
  res[i]=unlist(res0)[1]
}
res=rbind(k,res)
print(res)

## -----------------------------------------------------------------------------
na=444
nb=132
noo=361
nab=63
n=na+nb+noo+nab
p=q=r=L=rep(0,20)
p[1]=q[1]=r[1]=1/3
for(i in 2:20){
  naa=na*(p[i-1])^2/((p[i-1])^2+2*p[i-1]*r[i-1])
  nao=na*(2*p[i-1]*r[i-1])/((p[i-1])^2+2*p[i-1]*r[i-1])
  nab=nab
  nbb=nb*(q[i-1])^2/((q[i-1])^2+2*q[i-1]*r[i-1])
  nbo=nb*(2*q[i-1]*r[i-1])/((q[i-1])^2+2*q[i-1]*r[i-1])
  noo=noo
  p[i]=(2*naa+nao+nab)/(2*n)
  q[i]=(2*nbb+nbo+nab)/(2*n)
  r[i]=(2*noo+nao+nbo)/(2*n)
  L[i]=naa*log((p[i])^2)+nao*log(2*p[i]*r[i])+nbb*log((q[i])^2)+nbo*log(2*q[i]*r[i])+nab*log(p[i]*q[i])+noo*log((r[i])^2)
}
print(cbind(p,q,r,L))

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
lapply(formulas,lm)

## -----------------------------------------------------------------------------
res=list(0,0,0,0)
for(i in 1:length(formulas)){
  res[[i]]=lm(formulas[[i]])
}
print(res)

## -----------------------------------------------------------------------------
set.seed(12345)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(trials,function(x) return(x$p.value))

## -----------------------------------------------------------------------------
sapply(trials,"[[",3)

## -----------------------------------------------------------------------------
Mv=function(f,n,type,...){
  M=Map(f,...)
  return(vapply(M,cbind,type(n)))
}

## -----------------------------------------------------------------------------
lapply(mtcars, function(x) x / mean(x))[1:3]

## -----------------------------------------------------------------------------
M=Map("/",mtcars,colMeans(mtcars))

## -----------------------------------------------------------------------------
n=nrow(mtcars)
vapply(M,cbind,numeric(n))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
sourceCpp(code="#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix rwMC(double sigma,double x0,int N) {
  NumericMatrix mat(N, 2);
  mat(0,0)=x0;
  mat(0,1)=0;
  double y=0,u=0;
  for(int i = 2; i < N+1; i++) {
    y=rnorm(1,mat(i-2,0),sigma)[0];
    u=runif(1,0,1)[0];
    if(u<=exp(-abs(y))/exp(-abs(mat(i-2,0)))){
      mat(i-1,0)=y;
      mat(i-1,1)=mat(i-2,1);
    }
    else{
      mat(i-1,0)=mat(i-2,0);
      mat(i-1,1)=mat(i-2,1)+1;
    }
  }
  return(mat);
}")

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
N=2000
sigma=c(0.05,0.5,2,16)
x0=25
rw1=rw.Metropolis(sigma[1],x0,N)
rw2=rw.Metropolis(sigma[2],x0,N)
rw3=rw.Metropolis(sigma[3],x0,N)
rw4=rw.Metropolis(sigma[4],x0,N)
plot(rw1$x,type="l",lwd=2,xlab="sigma=0.05",ylab="X")
plot(rw2$x,type="l",lwd=2,xlab="sigma=0.5",ylab="X")
plot(rw3$x,type="l",lwd=2,xlab="sigma=2",ylab="X")
plot(rw4$x,type="l",lwd=2,xlab="sigma=16",ylab="X")
print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))

## -----------------------------------------------------------------------------
rwc1=rwMC(sigma[1],x0,N)
rwc2=rwMC(sigma[2],x0,N)
rwc3=rwMC(sigma[3],x0,N)
rwc4=rwMC(sigma[4],x0,N)
plot(rwc1[,1],type="l",lwd=2,xlab="sigma=0.05",ylab="X")
plot(rwc2[,1],type="l",lwd=2,xlab="sigma=0.5",ylab="X")
plot(rwc3[,1],type="l",lwd=2,xlab="sigma=2",ylab="X")
plot(rwc4[,1],type="l",lwd=2,xlab="sigma=16",ylab="X")
print(c(1-rwc1[N,2]/N,1-rwc2[N,2]/N,1-rwc3[N,2]/N,1-rwc4[N,2]/N))

## -----------------------------------------------------------------------------
qqplot(rw1$x,rwc1[,1])
qqplot(rw2$x,rwc2[,1])
qqplot(rw3$x,rwc3[,1])
qqplot(rw4$x,rwc4[,1])

## -----------------------------------------------------------------------------
ts1=microbenchmark(rw=rw.Metropolis(sigma[1],x0,N),rwc=rwMC(sigma[1],x0,N))
summary(ts1)[,c(1,3,5,6)]
ts2=microbenchmark(rw=rw.Metropolis(sigma[2],x0,N),rwc=rwMC(sigma[2],x0,N))
summary(ts2)[,c(1,3,5,6)]
ts3=microbenchmark(rw=rw.Metropolis(sigma[3],x0,N),rwc=rwMC(sigma[3],x0,N))
summary(ts3)[,c(1,3,5,6)]
ts4=microbenchmark(rw=rw.Metropolis(sigma[4],x0,N),rwc=rwMC(sigma[4],x0,N))
summary(ts4)[,c(1,3,5,6)]

