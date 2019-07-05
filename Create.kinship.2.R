


data.uni.ori = read.table("data.uni.ori.2.txt",header=TRUE)
data.biv.ori = read.table("data.biv.ori.2.txt",header=TRUE)

kin1 = NULL
kin2 = NULL

compteur=1

for (fam in unique(data.uni.ori$fam.id))
{
len = sum(data.uni.ori$fam.id==fam)
ind = data.uni.ori$ind.id[data.uni.ori$fam.id==fam]

V1 = diag(rep(0.5,len))
V2 = diag(rep(0.5,len))

for(i in 1:(len-1))
{
for(j in (i+1):len)
{
cond1 = (data.biv.ori$fam.id==fam)&(data.biv.ori$ind.id1==ind[i])&(data.biv.ori$ind.id2==ind[j])
cond2 = (data.biv.ori$fam.id==fam)&(data.biv.ori$ind.id2==ind[i])&(data.biv.ori$ind.id1==ind[j])
cond = cond1 | cond2
indice = which(cond)
V1[i,j] = (data.biv.ori$kinship.ff)[indice]
V1[j,i] = V1[i,j]
V2[i,j] = (data.biv.ori$kinship.fm)[indice]+(data.biv.ori$kinship.mf)[indice]+(data.biv.ori$kinship.mm)[indice]
V2[j,i] = V2[i,j]
}
}
kin1[[compteur]]=V1
kin2[[compteur]]=V2
compteur=compteur+1
}

