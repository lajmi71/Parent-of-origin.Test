
Create.data.sets=function(data.uni,data.biv)
{
len.uni = dim(data.uni)[1]
obs.time.proband = 0
KinshipWithProband = 0

for(i in 1:len.uni)
{
if(data.uni$ind.id[i]==1)
{
obs.time.proband[i]=data.uni$assessmentAge[i]
KinshipWithProband[i]=1
}
else
{
ind.proband.uni = which((data.uni$fam==data.uni$fam[i])&(data.uni$ind.id==1))
obs.time.proband[i]=data.uni$obs.time[ind.proband.uni]
ind.proband.biv = which((data.biv$fam.id==data.uni$fam.id[i])&(data.biv$ind.id1==1)&(data.biv$ind.id2==data.uni$ind.id[i]))
KinshipWithProband[i] = data.biv$kinship[ind.proband.biv]
} 
}
data.uni$obs.time.proband=obs.time.proband
data.uni$KinshipWithProband=KinshipWithProband

len.biv = dim(data.biv)[1]
ind1=NULL
ind2=NULL
for(i in 1:len.biv)
{
ind1 = c(ind1,which((data.uni$fam.id==data.biv$fam.id[i])&(data.uni$ind.id==data.biv$ind.id1[i])))
ind2 = c(ind2,which((data.uni$fam.id==data.biv$fam.id[i])&(data.uni$ind.id==data.biv$ind.id2[i])))
}

data.biv$obs.time1 = data.uni$obs.time[ind1]
data.biv$delta1 = data.uni$delta[ind1]
data.biv$obs.time2 = data.uni$obs.time[ind2]
data.biv$delta2 = data.uni$delta[ind2]


obs.time.proband = 0
assessment.age.proband = 0

for(i in 1:len.biv)
{
ind = which((data.uni$fam.id==data.biv$fam.id[i])&(data.uni$ind.id==1))
obs.time.proband[i] = data.uni$obs.time[ind]
assessment.age.proband[i] = data.uni$assessmentAge[ind]
}

data.biv$obs.time.proband = obs.time.proband
data.biv$assessment.age.proband = assessment.age.proband



kinship01=0
kinship02=0
for(i in 1:len.biv)
{
if(data.biv$ind.id1[i]==1)
{
kinship01[i] = 1
kinship02[i] = data.biv$kinship[i]
}
else
{
ind1 = which((data.biv$fam.id==data.biv$fam.id[i])&(data.biv$ind.id1==1)&(data.biv$ind.id2==data.biv$ind.id1[i]))
kinship01[i] = data.biv$kinship[ind1]
ind2 = which((data.biv$fam.id==data.biv$fam.id[i])&(data.biv$ind.id1==1)&(data.biv$ind.id2==data.biv$ind.id2[i]))
kinship02[i] = data.biv$kinship[ind2]
}
}

data.biv$kinship01 = kinship01
data.biv$kinship02 = kinship02


return(list(data.uni,data.biv))

}


