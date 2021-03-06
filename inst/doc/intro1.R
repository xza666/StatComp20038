## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
set.seed(12345)
N=10000
m0=10
m=2
d=BA(m0,m,N)
Fk=numeric(max(d))
for(i in 1:max(d)){
  Fk[i]=sum(as.integer(d>=i))
}
x=1:max(d)
plot(log(x),log(Fk[1:max(d)]))

## -----------------------------------------------------------------------------
x=2:54
plot(log(x),log(Fk[2:54]))
fit=lm(log(Fk[2:54])~log(x))
print(fit)
coef=round(fit$coef,3)
abline(a=coef[1],b=coef[2],col="red",lwd=2)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
set.seed(12345)
N=10000
beta=0.7
d=HV(N,beta)
Fk=numeric(max(d))
for(i in 1:max(d)){
  Fk[i]=sum(as.integer(d>=i))
}
x=1:max(d)
plot(log(x),log(Fk[1:max(d)]))

## -----------------------------------------------------------------------------
x=55:1050
plot(log(x),log(Fk[55:1050]))
fit=lm(log(Fk[55:1050])~log(x))
print(fit)
coef=round(fit$coef,3)
abline(a=coef[1],b=coef[2],col="red",lwd=2)

