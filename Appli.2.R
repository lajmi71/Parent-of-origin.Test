
source("Sim.2.R")
source("Create.kinship.2.R")
source("Prog.2.R")
source("Create.data.sets.2.R")

fam.id = data.uni.ori$fam.id
ind.id = data.uni.ori$ind.id
assessmentAge = data.uni.ori$assessmentAge

unique.fam.id = unique(fam.id)

shape = 3.5
scale = 140
theta=c(log(shape),log(scale))
h1=0.7
h2.real=0.2
h2.test=0.6

data.uni.sim = Generate.data(shape,scale,kin1,kin2,h1,h2.real,fam.id,ind.id,assessmentAge,unique.fam.id)
data=Create.data.sets(data.uni.sim,data.biv.ori)
data.uni=data[[1]]
data.biv=data[[2]]
theta.estim = optim(theta,Neg.LogLik.data.uni,h2=h2.test,data.uni=data.uni)$par
p.value=Compute.p.value(theta.estim,h2.test,data.uni,data.biv)

