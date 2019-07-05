library(VineCopula)
library(condMVNorm)
library(tmvtnorm)
library(numDeriv)


LogLik.vect.data.uni = function(theta,h2,data.uni)
{
shape = exp(theta[1])
scale = exp(theta[2])

indices1 = (data.uni$ind.id==1)  ### probands and thus observed & left-truncated
indices2 = (data.uni$ind.id>1)&(data.uni$delta==-1) ### non-proband & left-censored
indices3 = (data.uni$ind.id>1)&(data.uni$delta==0) ### non-proband & right-censored
indices4 = (data.uni$ind.id>1)&(data.uni$delta==1) ### non-proband & observed

N = dim(data.uni)[1]

LogLik = rep(NA,N)

S1 = pweibull(data.uni$obs.time.proband,shape=shape,scale=scale,lower.tail=FALSE)
S2 = pweibull(data.uni$obs.time,shape=shape,scale=scale,lower.tail=FALSE)
f2 = dweibull(data.uni$obs.time,shape=shape,scale=scale,log=TRUE)

LogLik[indices1] = dweibull(data.uni$obs.time[indices1],shape=shape,scale=scale,log=TRUE) - pweibull(data.uni$assessmentAge[indices1],shape=shape,scale=scale,log.p=TRUE)
LogLik[indices2] = log(1-BiCopHfunc1(u1=S1[indices2],u2=S2[indices2],family=1,par=h2*data.uni$KinshipWithProband[indices2]))
LogLik[indices3] = log(BiCopHfunc1(u1=S1[indices3],u2=S2[indices3],family=1,par=h2*data.uni$KinshipWithProband[indices3]))
LogLik[indices4] = log(BiCopPDF(u1=S1[indices4],u2=S2[indices4],family=1,par=h2*data.uni$KinshipWithProband[indices4])) + f2[indices4]
LogLik.fam = NULL

for(fam in unique(data.uni$fam.id))
{
LogLik.fam = c(LogLik.fam,sum(LogLik[data.uni$fam.id==fam]))
}

return(LogLik.fam)
}

Neg.LogLik.data.uni = function(theta,h2,data.uni) {-sum(LogLik.vect.data.uni(theta,h2,data.uni))}

############################################################################################################

Compute.stat.terms = function(theta,h2,data.biv)
{
shape = exp(theta[1])
scale = exp(theta[2])

L1 = qnorm(pweibull(data.biv$obs.time1,shape=shape,scale=scale,lower.tail=FALSE))
L2 = qnorm(pweibull(data.biv$obs.time2,shape=shape,scale=scale,lower.tail=FALSE))

L.proband = qnorm(pweibull(data.biv$obs.time.proband,shape=shape,scale=scale,lower.tail=FALSE))
L.lim.proband = qnorm(pweibull(data.biv$assessment.age.proband,shape=shape,scale=scale,lower.tail=FALSE))

N = dim(data.biv)[1]

mult = data.biv$kinship.ff-(data.biv$kinship.mf+data.biv$kinship.fm+data.biv$kinship.mm)

res = rep(NA,N)

indices0 = (mult==0)

indices1 = ((data.biv$ind.id1==1) & (data.biv$delta2==-1) & (mult!=0))
indices2 = ((data.biv$ind.id1==1) & (data.biv$delta2==0) & (mult!=0))
indices3 = ((data.biv$ind.id1==1) & (data.biv$delta2==1) & (mult!=0))

indices4 = ((data.biv$ind.id1>1) & (data.biv$delta1==-1) & (data.biv$delta2==-1) & (mult!=0))
indices5 = ((data.biv$ind.id1>1) & (data.biv$delta1==-1) & (data.biv$delta2==0) & (mult!=0))
indices6 = ((data.biv$ind.id1>1) & (data.biv$delta1==-1) & (data.biv$delta2==1) & (mult!=0))

indices7 = ((data.biv$ind.id1>1) & (data.biv$delta1==0) & (data.biv$delta2==-1) & (mult!=0))
indices8 = ((data.biv$ind.id1>1) & (data.biv$delta1==0) & (data.biv$delta2==0) & (mult!=0))
indices9 = ((data.biv$ind.id1>1) & (data.biv$delta1==0) & (data.biv$delta2==1) & (mult!=0))

indices10 = ((data.biv$ind.id1>1) & (data.biv$delta1==1) & (data.biv$delta2==-1) & (mult!=0))
indices11 = ((data.biv$ind.id1>1) & (data.biv$delta1==1) & (data.biv$delta2==0) & (mult!=0))
indices12 = ((data.biv$ind.id1>1) & (data.biv$delta1==1) & (data.biv$delta2==1) & (mult!=0))

res[indices0] = 0

if (sum(indices1)>0) {
res[indices1] = sapply(1:(sum(indices1)),fct1,L1[indices1],L2[indices1],L.lim.proband[indices1],data.biv$kinship[indices1],h2)}

if (sum(indices2)>0) {
res[indices2] = sapply(1:(sum(indices2)),fct2,L1[indices2],L2[indices2],L.lim.proband[indices2],data.biv$kinship[indices2],h2)}

if (sum(indices3)>0) {
res[indices3] = sapply(1:(sum(indices3)),fct3,L1[indices3],L2[indices3],L.lim.proband[indices3],data.biv$kinship[indices3],h2)}


if (sum(indices4)>0) {
res[indices4] = sapply(1:(sum(indices4)),fct4,L1[indices4],L2[indices4],L.proband[indices4],data.biv$kinship[indices4],data.biv$kinship01[indices4],data.biv$kinship02[indices4],h2)}

if (sum(indices5)>0) {
res[indices5] = sapply(1:(sum(indices5)),fct5,L1[indices5],L2[indices5],L.proband[indices5],data.biv$kinship[indices5],data.biv$kinship01[indices5],data.biv$kinship02[indices5],h2)}

if (sum(indices6)>0) {
res[indices6] = sapply(1:(sum(indices6)),fct6,L1[indices6],L2[indices6],L.proband[indices6],data.biv$kinship[indices6],data.biv$kinship01[indices6],data.biv$kinship02[indices6],h2)}

if (sum(indices7)>0) {
res[indices7] = sapply(1:(sum(indices7)),fct7,L1[indices7],L2[indices7],L.proband[indices7],data.biv$kinship[indices7],data.biv$kinship01[indices7],data.biv$kinship02[indices7],h2)}

if (sum(indices8)>0) {
res[indices8] = sapply(1:(sum(indices8)),fct8,L1[indices8],L2[indices8],L.proband[indices8],data.biv$kinship[indices8],data.biv$kinship01[indices8],data.biv$kinship02[indices8],h2)}

if (sum(indices9)>0) {
res[indices9] = sapply(1:(sum(indices9)),fct9,L1[indices9],L2[indices9],L.proband[indices9],data.biv$kinship[indices9],data.biv$kinship01[indices9],data.biv$kinship02[indices9],h2)}

if (sum(indices10)>0) {
res[indices10] = sapply(1:(sum(indices10)),fct10,L1[indices10],L2[indices10],L.proband[indices10],data.biv$kinship[indices10],data.biv$kinship01[indices10],data.biv$kinship02[indices10],h2)}

if (sum(indices11)>0) {
res[indices11] = sapply(1:(sum(indices11)),fct11,L1[indices11],L2[indices11],L.proband[indices11],data.biv$kinship[indices11],data.biv$kinship01[indices11],data.biv$kinship02[indices11],h2)}

if (sum(indices12)>0) {
res[indices12] = sapply(1:(sum(indices12)),fct12,L1[indices12],L2[indices12],L.proband[indices12],data.biv$kinship[indices12],data.biv$kinship01[indices12],data.biv$kinship02[indices12],h2)}


res = res*mult

res.fam = NULL
for(fam in unique(data.biv$fam.id))
{
res.fam = c(res.fam,sum(res[data.biv$fam.id==fam]))
}
return(res.fam)
}



fct1 = function(i,L1,L2,L.lim.proband,kinship,h2)
{
v = matrix(c(1,h2*kinship[i],h2*kinship[i],1),ncol=2)
m=c(0,0)
R1 = condMVN(mean=m,sigma=v,dependent.ind=2,given.ind=1,X.given=L1[i])
condMean1 = R1$condMean
condVar1   = R1$condVar
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,lower=L2[i])
condMean2 = R2$tmean
observed.stat = L1[i]*condMean2
R3 = mtmvnorm(mean=m,sigma=v,lower=c(L.lim.proband[i],-Inf))
condMean3 = R3$tmean
condVar3 = R3$tvar
expected.stat = prod(R3$tmean)+R3$tvar[1,2]
return(observed.stat-expected.stat)
}

fct2 = function(i,L1,L2,L.lim.proband,kinship,h2)
{
v = matrix(c(1,h2*kinship[i],h2*kinship[i],1),ncol=2)
m=c(0,0)
R1 = condMVN(mean=m,sigma=v,dependent.ind=2,given.ind=1,X.given=L1[i])
condMean1 = R1$condMean
condVar1   = R1$condVar
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,upper=L2[i])
condMean2 = R2$tmean
observed.stat = L1[i]*condMean2
R3 = mtmvnorm(mean=m,sigma=v,lower=c(L.lim.proband[i],-Inf))
condMean3 = R3$tmean
condVar3 = R3$tvar
expected.stat = prod(R3$tmean)+R3$tvar[1,2]
return(observed.stat-expected.stat)
}

fct3 = function(i,L1,L2,L.lim.proband,kinship,h2)
{
observed.stat = L1[i]*L2[i]
v = matrix(c(1,h2*kinship[i],h2*kinship[i],1),ncol=2)
m=c(0,0)
R = mtmvnorm(mean=m,sigma=v,lower=c(L.lim.proband[i],-Inf))
condMean = R$tmean
condVar = R$tvar
expected.stat = prod(R$tmean)+R$tvar[1,2]
return(observed.stat-expected.stat)
}

fct4 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,lower=c(L1[i],L2[i]))
condMean2 = R2$tmean
condVar2 = R2$tvar
observed.stat = prod(condMean2) + condVar2[1,2]
return(observed.stat-expected.stat)
}


fct5 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,lower=c(L1[i],-Inf),upper=c(Inf,L2[i]))
condMean2 = R2$tmean
condVar2 = R2$tvar
observed.stat = prod(condMean2) + condVar2[1,2]
return(observed.stat-expected.stat)
}

fct6 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = condMVN(mean=condMean1,sigma=condVar1,dependent.ind=1,given.ind=2,X.given=L2[i])
condMean2 = R2$condMean
condVar2 = R2$condVar
R3 = mtmvnorm(mean=condMean2,sigma=condVar2,lower=L2[i])
observed.stat = R3$tmean*L2[i]
return(observed.stat-expected.stat)
}


fct7 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,lower=c(-Inf,L2[i]),upper=c(L1[i],Inf))
condMean2 = R2$tmean
condVar2 = R2$tvar
observed.stat = prod(condMean2) + condVar2[1,2]
return(observed.stat-expected.stat)
}

fct8 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = mtmvnorm(mean=condMean1,sigma=condVar1,upper=c(L1[i],L2[i]))
condMean2 = R2$tmean
condVar2 = R2$tvar
observed.stat = prod(condMean2) + condVar2[1,2]
return(observed.stat-expected.stat)
}

fct9 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = condMVN(mean=condMean1,sigma=condVar1,dependent.ind=1,given.ind=2,X.given=L2[i])
condMean2 = R2$condMean
condVar2 = R2$condVar
R3 = mtmvnorm(mean=condMean2,sigma=condVar2,upper=L1[i])
observed.stat = R3$tmean*L2[i]
return(observed.stat-expected.stat)
}

fct10 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = condMVN(mean=condMean1,sigma=condVar1,dependent.ind=2,given.ind=1,X.given=L1[i])
condMean2 = R2$condMean
condVar2 = R2$condVar
R3 = mtmvnorm(mean=condMean2,sigma=condVar2,lower=L2[i])
observed.stat = L1[i]*R3$tmean
return(observed.stat-expected.stat)
}


fct11 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
R2 = condMVN(mean=condMean1,sigma=condVar1,dependent.ind=2,given.ind=1,X.given=L1[i])
condMean2 = R2$condMean
condVar2 = R2$condVar
R3 = mtmvnorm(mean=condMean2,sigma=condVar2,upper=L2[i])
observed.stat = L1[i]*R3$tmean
return(observed.stat-expected.stat)
}

fct12 = function(i,L1,L2,L.proband,kinship,kinship01,kinship02,h2)
{
m = c(0,0,0)
v = matrix(c(1,h2*kinship01[i],h2*kinship02[i],h2*kinship01[i],1,h2*kinship[i],h2*kinship02[i],h2*kinship[i],1),ncol=3)
R1 = condMVN(mean=m,sigma=v,dependent.ind=c(2,3),given.ind=1,X.given=L.proband[i])
condMean1 = R1$condMean
condVar1 = R1$condVar
expected.stat = prod(condMean1) + condVar1[1,2]
observed.stat = L1[i]*L2[i]
return(observed.stat-expected.stat)
}

Compute.stat.terms.means = function(theta,h2,data.biv) {mean(Compute.stat.terms(theta,h2,data.biv),na.rm=TRUE)}

Compute.p.value=function(theta.estim,h2,data.uni,data.biv)
{
RR1 = genD(LogLik.vect.data.uni,theta.estim,h2=h2,data.uni=data.uni)$D
l = t(RR1[,1:2])
B = solve(matrix(apply(RR1[,3:5],2,mean)[c(1,2,2,3)],ncol=2))
A = jacobian(Compute.stat.terms.means,theta.estim,h2=h2,data.biv=data.biv)
stat.terms = Compute.stat.terms(theta.estim,h2,data.biv)
add.term = A%*%B%*%l
Z.obs=sum(stat.terms,na.rm=TRUE)/sqrt(sum((stat.terms-add.term)^2,na.rm=TRUE))
p.value = 2*(1-pnorm(abs(Z.obs)))
return(p.value)
}

