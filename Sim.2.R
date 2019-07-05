library(truncdist)
library(condMVNorm)

Generate.T.proband = function(shape,scale,assessmentAgeProband) {rtrunc(1,spec="weibull",a=0,b=assessmentAgeProband,shape=shape,scale=scale)}

Generate.T.nonproband = function(shape,scale,T.proband,K1,K2,h1,h2)
{
len = dim(K1)[1]
V = 2*h2*(h1*K1+(1-h1)*K2) + (1-h2)*diag(rep(1,len))
U.proband = pweibull(T.proband,shape=shape,scale=scale,lower.tail=FALSE)
Z.proband = qnorm(U.proband)
Z.nonproband = rcmvnorm(1,mean=rep(0,len),sigma=V,dependent.ind=(2:len),given.ind=1,X.given=Z.proband)
U.nonproband = pnorm(Z.nonproband)
T.nonproband = qweibull(U.nonproband,shape=shape,scale=scale,lower.tail=FALSE)
T.nonproband
}

Generate.T.family = function(shape,scale,assessmentAgeProband,K1,K2,h1,h2)
{
T.proband = Generate.T.proband(shape,scale,assessmentAgeProband)
T.nonproband = Generate.T.nonproband(shape,scale,T.proband,K1,K2,h1,h2)
c(T.proband,T.nonproband)
}

Generate.data = function(shape,scale,kin1,kin2,h1,h2,fam.id,ind.id,assessmentAge,unique.fam.id)
{
Generated.T = NULL
compteur=1
for(fam in unique.fam.id)
{
K1 = kin1[[compteur]]
K2 = kin2[[compteur]]
assessmentAgeProband = assessmentAge[(fam.id==fam)&(ind.id==1)]
T.family = Generate.T.family(shape,scale,assessmentAgeProband,K1,K2,h1,h2)
Generated.T = c(Generated.T,T.family)
compteur = compteur + 1 
}

obs.time = round(pmin(Generated.T,assessmentAge))
delta = as.numeric(Generated.T<assessmentAge)

indices.to.choose.from = which((delta==1)&(ind.id>1))
indices.left.censored = sample(indices.to.choose.from,2)
obs.time[indices.left.censored] = assessmentAge[indices.left.censored]
delta[indices.left.censored]= rep(-1,2)

data.uni.sim = data.frame(fam.id=fam.id,ind.id=ind.id,obs.time=obs.time,delta=delta,assessmentAge=assessmentAge)

return(data.uni.sim)
}