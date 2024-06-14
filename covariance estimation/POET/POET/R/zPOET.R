POET<-function(Y,K=-Inf,C=-Inf, thres='soft', matrix='cor') {

p=nrow(Y)
n=ncol(Y)
Y<-Y-t(t(apply(t(Y),2,mean)))%*%matrix(1,1,n) # Y is de-meaned

if (K==-Inf){
K1=0.25*(POETKhat(Y)$K1HL+POETKhat(Y)$K2HL+POETKhat(Y)$K1BN+POETKhat(Y)$K2BN)
K=floor(K1)+1
}

if (K>0){
V<-eigen(t(Y)%*%Y)$vectors
V=as.matrix(V)
Dd<-eigen(t(Y)%*%Y)$values
Dd=as.vector(Dd)
W<-sort(diag(Dd),index.return=TRUE)$x
W=as.matrix(W)
Id<-sort(diag(Dd),index.return=TRUE)$ix
Id=as.matrix(Id)
F<-sqrt(n)*V[,1:K]  #F is n by K
LamPCA=Y%*%F/n
uhat=Y-LamPCA%*%t(F)  # p by n
Lowrank=LamPCA%*%t(LamPCA)
    rate=1/sqrt(p)+sqrt((log(p))/n)
} else {
uhat=Y   # Sigma_y itself is sparse

    rate=sqrt((log(p))/n)
    Lowrank=matrix(0,p,p)
}



SuPCA=uhat%*%t(uhat)/n
SuDiag=diag(diag(SuPCA))
if (matrix =='cor'){
    R=solve(SuDiag^(1/2))%*%SuPCA%*%solve(SuDiag^(1/2))
}
if (matrix  =='vad') {
    R=SuPCA
}

if (C==-Inf){
C1=POETCmin(Y,K,thres,matrix)
C=C1+0.1
}

uu=array(0,dim=c(p,p,n))
roottheta=array(0,dim=c(p,p))
lambda=array(0,dim=c(p,p))
for (i in 1:p)  # adaptive threshold
{for (j in 1:i)
{uu[i,j,]=uhat[i,]*uhat[j,]
roottheta[i,j]=sd(uu[i,j,])
lambda[i,j]=roottheta[i,j]*rate*C
lambda[j,i]=lambda[i,j]
}
}



Rthresh=matrix(0,p,p)

if (thres=='soft') {
    for (i in 1:p) {
        for (j in 1:i) {
            if (abs(R[i,j])<lambda[i,j]  && j<i) {
                Rthresh[i,j]=0
            } 
            else {if (j==i) { 
                    Rthresh[i,j]=R[i,j]
                }
                else {
                    Rthresh[i,j]=sign(R[i,j])*(abs(R[i,j])-lambda[i,j]);
                }
            }
            Rthresh[j,i]=Rthresh[i,j];
        }
    }
}



if (thres=='hard') {
    for (i in 1:p) {
        for (j in 1:i) {
            if (abs(R[i,j])<lambda[i,j]  && j<i) {
                Rthresh[i,j]=0
            }
            else  {
                Rthresh[i,j]=R[i,j]
            }
            Rthresh[j,i]=Rthresh[i,j]
        }
    }
}



if (thres=='scad') {
    for (i in 1:p) {
        for (j in 1:i) {
             if (j==i) {Rthresh[i,j]=R[i,j]
             }
             else {if (abs(R[i,j])<lambda[i,j]){ 
                        Rthresh[i,j]=0 
                   }
                 else {if (abs(R[i,j])<2*lambda[i,j]) {
                         Rthresh[i,j]=sign(R[i,j])*(abs(R[i,j])-lambda[i,j])
                      }
                     else {if (abs(R[i,j])<3.7*lambda[i,j]) { 
                             Rthresh[i,j]=((3.7-1)*R[i,j]-sign(R[i,j])*3.7*lambda[i,j])/(3.7-2)
                              }
                         else {Rthresh[i,j]=R[i,j]
                         }
                     }
                 }
             }
             Rthresh[j,i]=Rthresh[i,j]
        }
    }
}

SigmaU=matrix(0,p,p)
if (matrix =='cor') {
    SigmaU=SuDiag^(1/2)%*%Rthresh*SuDiag^(1/2)
}
if (matrix =='vad')     {
SigmaU=Rthresh
}


SigmaY=SigmaU+Lowrank


result <-list(SigmaU=SigmaU,SigmaY=SigmaY,factors=t(F),loadings=LamPCA) 
return(result)
}



POETCmin<-function(Y,K,thres,matrix) {
p=nrow(Y)
n=ncol(Y)

mineig<-function(Y,K,C,thres,matrix){
SigmaU=POET(Y,K,C,thres,matrix)$SigmaU
f=min(eigen(SigmaU)$values)
result<-c(f)
return(result)
}

f<-function(x) mineig(Y,K,x,thres,matrix)

if (f(50)*f(-50)<0){
r<-uniroot(f,c(-50,50),tol=0.001)
result<-max(0,c(r)$root)
return(result)
}

else {
C=0
result<-C
return(result)
}
}

POETKhat<-function(Y){
p=nrow(Y)
n=ncol(Y)
Y<-Y-t(t(apply(t(Y),2,mean)))%*%matrix(1,1,n) # Y is de-meaned

#Hallin and Liska method

c=seq(0.05,5,length=100)
re=20
rmax=10
IC=array(0,c(2,re,rmax,100))  
gT1HL=array(1,c(re))
gT2HL=array(1,c(re))
pi=array(1,c(re))
ni=array(1,c(re))

for (i in 1:re) {    #generate the subsets, "re" of them
pi[i]=min(i*floor(p/re)+min(p,5),p)
ni[i]=min(i*floor(n/re)+min(n,5),n)
if (i==re){
pi[i]=p
ni[i]=n
}
Yi=Y[1:pi[i],1:ni[i]]
frob=array(0,c(rmax))
penal=array(0,c(rmax))

for (k in 1:min(pi[i],ni[i],rmax)) {
V<-eigen(t(Yi)%*%Yi)$vectors   
V=as.matrix(V)
Dd<-eigen(t(Yi)%*%Yi)$values
Dd=as.vector(Dd)

F<-V[,1:k]


LamPCA=Yi%*%F/ni[i]
uhat=Yi-LamPCA%*%t(F)  # pi by ni
frob[k]=sum(diag(uhat%*%t(uhat)))/(pi[i]*ni[i])
gT1HL[i]=log((pi[i]*ni[i])/(pi[i]+ni[i]))*(pi[i]+ni[i])/(pi[i]*ni[i])
gT2HL[i]=log(min(pi[i],ni[i]))*(pi[i]+ni[i])/(pi[i]*ni[i])



for (l in 1:100){    # only fills in the ICs up to k, which may be <rmax
IC[1,i,k,l]=log(frob[k])+c[l]*k*gT1HL[i]
IC[2,i,k,l]=log(frob[k])+c[l]*k*gT2HL[i]
}

}
}

rhat=array(0,c(2,re,100))
for (i in 1:re){
for (l in 1:100) {
m=min(pi[i],ni[i],rmax);
temp1=which.min(IC[1,i,1:m,l])
rhat[1,i,l]=temp1
temp2=which.min(IC[2,i,1:m,l])
rhat[2,i,l]=temp2
}
}

Sc1=array(0,c(100))
Sc2=array(0,c(100))

for (l in 1:100){
Sc1[l]=sd(rhat[1,,l])
Sc2[l]=sd(rhat[2,,l])

}
 
c1vec=which(Sc1==0)
ctemp1=c1vec[1]
c1=c[ctemp1]         # constant we choose in the penalty function
K1HL=rhat[1,1,ctemp1]   #all Ks in that row are equal

c2vec=which(Sc2==0)
ctemp2=c2vec[1]
c2=c[ctemp2]       
K2HL=rhat[2,1,ctemp2]


# Bai and Ng method - penalty corresponds to c=1

c=1
rmax=10
IC=array(0,c(2,rmax))  
frob=array(0,c(rmax))
penal=array(0,c(rmax))

for (k in 1:rmax) {
V<-eigen(t(Y)%*%Y)$vectors   
V=as.matrix(V)
Dd<-eigen(t(Y)%*%Y)$values
Dd=as.vector(Dd)

F<-V[,1:k]

LamPCA=Y%*%F/n
uhat=Y-LamPCA%*%t(F)  # p by n
frob[k]=sum(diag(uhat%*%t(uhat)))/(p*n)
gT1BN=log((p*n)/(p+n))*(p+n)/(p*n)
gT2BN=log(min(p,n))*(p+n)/(p*n)
IC[1,k]=log(frob[k])+k*gT1BN
IC[2,k]=log(frob[k])+k*gT2BN
}
K1BN=which.min(IC[1,])
K2BN=which.min(IC[2,])


result<-list(K1HL=K1HL,K2HL=K2HL,K1BN=K1BN,K2BN=K2BN,IC=IC)
return(result)

}







