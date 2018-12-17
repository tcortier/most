library("ppcor")
library(lsmeans)
#Pour tester l'effet sur un paramètre suelement
library(multcomp)
library(pROC)
library(plyr)
library(nlme)

#Importation

setwd("~/Documents/master2/most/Projet_2018-2019/")
foodweb<-read.table("135FoodWebs.csv",header=TRUE,sep=";")
zoo<-read.table("135Zoocoenoses.csv",header=TRUE,sep=";")

#Creation d'une nouvelle table avec les informations des deux tables
zf<-merge(foodweb, zoo, by.x = "Web.ID",by.y =  "Web.ID")
zf$N.tot<-as.numeric(zf$N.tot)
zf$Vegetation.cover....<-as.numeric(zf$Vegetation.cover....)
zf$Soil.Ni<-as.numeric(zf$Soil.Ni)
zf$Soil.Ptot<-as.numeric(zf$Soil.Ptot)
zf$Soil.Cd<-as.numeric(zf$Soil.Cd)
ind<-c(95:138,140:173,175:194,196:232)

for (i in 10:34){
  v<-c()
  for (j in ind){
    e=zf[zf$Web.ID==j,][1,i]
    v<-c(v,e)
  }
  zoo[names(zf)[i]]<-v
}

#Etude de la répartition des espèces par site

stat<-as.data.frame(table(zf$Genus.Morphon,zf$Web.ID))

count(zf,var="Web.ID")

par(mfrow = c(1,1))

plot(zoo$Soil.pH~zoo$Description,las=2)
boxplot(zoo$Soil.pH~zoo$Description,las=2)
boxplot(zf$Log.Abundance.~zf$Description,las=2)
boxplot(zf$Log.averageMass.~zf$Description,las=2)
boxplot(zf$Log.Biomass.~zf$Description,las=2)
boxplot(zf$C.tot~zf$Description,las=2)
boxplot(zf$N.tot~zf$Description,las=2)
boxplot(zf$Mean.rainfall~zf$Description,las=2)
boxplot(zf$Average.T~zf$Description,las=2)
boxplot(zf$Vegetation.cover....~zf$Description,las=2)
boxplot(zoo$Average.T~zoo$Description,las=2)
boxplot(zoo$Average.Tmax~zoo$Description,las=2)


#PCA results and soil quality separation

library(FactoMineR)
library(factoextra)

names(zoo)
str(zf)


#On fait une PCA sur tous les types d'habitats
respcapro = FactoMineR::PCA(X = zoo, # the data set used. Rows are individuals and columns are numeric variables.
                         scale.unit = TRUE,  #if TRUE, the data are scaled to unit variance before the analysis. 
                         quali.sup = c(1:29,34,35,48,49,50,51,53),
                         graph = F, #if TRUE, the graphs are displayed.
                         ncp = 5) #indexes of the annotation columns

eig <- get_eig(respcapro)
fviz_screeplot(respcapro, addlabels = TRUE)
ind <- get_pca_ind(respcapro)
ind$coord
ind$contrib
fviz_pca_ind(respcapro, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
fviz_pca_biplot(respcapro, axes = c(1,2), habillage = 'Description',select.var = list(contrib = 6),invisible = 'quali')
#corrplot(respcapro$var$cos2)
barplot(respcapro$var$contrib[,1],las=2)
barplot(respcapro$var$contrib[,2],las=2)
respcapro$var$contrib
fviz_pca_var(respcapro, axes = c(1,2), choix = 'var', select.var = list(contrib = 6))

zoo$PCA.dim.1<-respcapro$ind$coord[,1]

#On va essayer de faire la même chose mais qu'avec les fermes
fermes=zoo[(zoo$Ecosystem.Type.ID==1)|(zoo$Ecosystem.Type.ID==2)|(zoo$Ecosystem.Type.ID==3)|(zoo$Ecosystem.Type.ID==4),]
respcafermes = FactoMineR::PCA(X = fermes, # the data set used. Rows are individuals and columns are numeric variables.
                            scale.unit = TRUE,  #if TRUE, the data are scaled to unit variance before the analysis. 
                            quali.sup = c(1:29,34,35,48,49,50,51,53,54,55,56),
                            graph = F, #if TRUE, the graphs are displayed.
                            ncp = 5) #indexes of the annotation columns

eig <- get_eig(respcafermes)
fviz_screeplot(respcafermes, addlabels = TRUE)
ind <- get_pca_ind(respcafermes)
ind$coord
ind$contrib
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
fviz_pca_biplot(respcafermes, axes = c(1,2), habillage = 'Description',select.var = list(contrib = 6),invisible = 'quali')
#corrplot(respcapro$var$cos2)

barplot(respcafermes$var$contrib[,1],las=2)
barplot(respcafermes$var$contrib[,2],las=2)
respcafermes$var$contrib
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
barplot(respcafermes$var$contrib[,2],las=2)
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
fviz_pca_var(respcafermes, axes = c(1,2), choix = 'var', select.var = list(contrib = 8))



count(zf$Feeding.Preference)
count(zf$Genus.Morphon)

##Test statistique du phosphate selon les fermes

mod.phos<-lm(Soil.phosphate~Description,data=fermes)

par(mfrow = c(2,2))
plot(mod.phos)
summary(mod.phos)
anova(mod.phos)
pairwise.t.test(fermes$Soil.phosphate,fermes$Description,p.adjust.method = "holm")
ggplot(fermes, aes(Description,Soil.phosphate))+ geom_dotplot(binaxis = "y",stackdir = "center",alpha=0.2,col=2)+geom_violin(alpha=0.2,col=2,draw_quantiles = c(0.5))

#Test avec l'azote

ggplot(fermes, aes(Description,Total.N.input))+ geom_dotplot(binaxis = "y",stackdir = "center",alpha=0.2,col=2)+geom_violin(alpha=0.2,col=2,draw_quantiles = c(0.5))
fermes["log.total.N"]<-fermes
mod.N<-lm(Total.N.input~Description,data=fermes)

par(mfrow = c(2,2))
plot(mod.N)
summary(mod.N)
anova(mod.N)

##Influence avec le nombre de taxas

boxplot(zf$Taxa.S~zf$Description,las=2)

mod1<-lm(Taxa.S~Description,data=zoo)
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)
pairwise.t.test(zoo$Taxa.S,zoo$Description,p.adjust.method = "bonferroni")

#Création d'une variable par site représentantant la nombre de feeding preference différent selon le sampling et le nombre d'invertébrés différents
#Dans un premier temps, nous allons conostruire un vecteurs avec tous les cathégories d'invertébrés différentes. 

ind<-c(95:138,140:173,175:194,196:232)
shannon.trophic<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]$Trophic.ID
  w=zf[zf$Web.ID==j,]$Log.Biomass.
  feed<-list(0,0,0,0,0,0,0,0,0)
  somme=0
  sh=0
  for (i in 1:length(v)){
    feed[[v[i]%%10]]<-feed[[v[i]%%10]]+exp(w[i])
    somme=somme+exp(w[i])
  }
  for (k in 1:9){
    if (feed[[k]]!=0){
      r=feed[[k]]/somme
      sh=sh-r*log2(r)
    }
    
  }
  shannon.trophic<-c(shannon.trophic,sh)
}

zoo["shannon.trophic"]<-shannon.trophic

ggplot(zoo, aes(Description,shannon.trophic ))+ geom_dotplot(binaxis = "y",stackdir = "center",alpha=0.2,col=2)+geom_violin(alpha=0.2,col=2)

plot(zoo$shannon.trophic~zoo$PCA.dim.1,col=zoo$Description)
plot(zoo$Taxa.S~zoo$PCA.dim.1,col=zoo$Description)
ggplot(zoo, aes(PCA.dim.1,shannon.trophic ))+geom_smooth(aes(color = Description), model = lm)
ggplot(zoo, aes(PCA.dim.1,shannon.trophic ))+geom_smooth( model = lm)

#Clairement dans les analyses, il faut mettre la température en facteur et peut être la distance

#Biomassse des mématodes en fonction des caractéristiques du sol. 

ind<-c(95:138,140:173,175:194,196:232)
nem<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]
  s=0
  for (i in 1:length(v$Web.ID)){
    if (v$Trophic.ID[i]%%10==1){
      s=s+exp(v$Log.Biomass.[i])
      print(s)
    }
  }
  nem<-c(nem,log(s))
}

zoo["biom.nematodes"]<-nem

ind<-c(95:138,140:173,175:194,196:232)
nem<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]
  s=0
  for (i in 1:length(v$Web.ID)){
    if (v$Trophic.ID[i]%%10==1){
      s=s+exp(v$Log.Biomass.[i])
      print(s)
    }
  }
  nem<-c(nem,log(s))
}

plot(zoo$biom.nematodes~zoo$Description,las=2)

variables<- c(2,30:33,36:47,52,54)
zoo0<-zoo[zoo$Ecosystem.Type.ID!=6,]
zoo2<-zoo0[variables]


pcor(zoo2)

mod1<-lm(biom.nematodes ~ .,data=zoo2)
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)

mod2<-step(mod1,scope=~ . ,direction="both")

moddef<-lm(biom.nematodes ~ Description + Airborne.N + Total.N.input + Soil.Cd + 
             Soil.Hg + Soil.Pb + Soil.Zn + Max.rainfall,data=zoo2 )
par(mfrow = c(2,2))
plot(moddef)
summary(moddef)
anova(moddef)
par(mfrow = c(1,1))

plot(biom.nematodes~Soil.Hg,col=Description,data=zoo0)
legend(x='bottomright',legend=levels(zoo$Description),col=zoo$Description,pch=21,cex=1)
plot(zoo$biom.nematodes~zoo$Soil.phosphate,col=zoo0$Description)

#La biomasse de nématodes en fonction des sites

ggplot(zoo[zoo$Ecosystem.Type.ID!=6,], aes(Description,biom.nematodes ))+ geom_dotplot(binaxis = "y",stackdir = "center",alpha=0.2,col=2)+geom_violin(alpha=0.2,col=2)


mod1<-lm(biom.nematodes ~ Description,data=zoo0)
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)

#Modèle avec le type de sol en effet mixte

mod.mixte<-lme(biom.nematodes ~  Airborne.N + Total.N.input + Soil.Cd + 
                 Soil.Hg + Soil.Pb + Soil.Zn + Max.rainfall,data=zoo2,random =~1|Description,method = "ML")
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)
anova(mod1)

#Modèle avec le type d'utilisation et le résultat de la PCA sur les paramètres du sol

plot(zoo$biom.nematodes~zoo$PCA.dim.1,col=zoo$Description)

ggplot(zoo, aes(PCA.dim.1,biom.nematodes ))+geom_smooth( model = lm)

#Analyses simples

ind<-c(95:138,140:173,175:194,196:232)
bio<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]
  s=0
  for (i in 1:length(v$Web.ID)){
    s=s+exp(v$Log.Biomass.[i])
  }
  bio<-c(bio,log(s))
}

zoo["biom.totale"]<-bio

plot(zoo$biom.totale~zoo$Description,las=2)

ind<-c(95:138,140:173,175:194,196:232)
bio<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]
  s=0
  for (i in 1:length(v$Web.ID)){
    s=s+exp(v$Log.averageMass.[i])
  }
  bio<-c(bio,log(s))
}

ind<-c(95:138,140:173,175:194,196:232)
bio<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]
  s=0
  for (i in 1:length(v$Web.ID)){
    s=s+exp(v$Log.Abundance.[i])
  }
  bio<-c(bio,s)
}

zoo["abundance"]<-bio

ggplot(zoo, aes(PCA.dim.1, biom.totale))+geom_smooth( model = lm)

ggplot(zoo, aes(Total.N.input, biom.totale))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
ggplot(zoo, aes(Total.N.input,Taxa.S))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
ggplot(zoo, aes(Soil.pH,Taxa.S))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
ggplot(zoo, aes(Airborne.N,biom.totale))+geom_smooth(aes(color = Description), model = lm,se=FALSE)

plot(zoo$Total.N.input~zoo$Description,las=2)

variables<- c(2,30:33,36:47,52,58)
zoo2<-zoo[variables]


pcor(zoo2)

mod1<-lm(average.masse ~ .,data=zoo2)
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)
anova(mod1)

zoo0<-zoo[zoo$Ecosystem.Type.ID!=6,]
fermes=zoo[(zoo$Ecosystem.Type.ID==1)|(zoo$Ecosystem.Type.ID==2)|(zoo$Ecosystem.Type.ID==3)|(zoo$Ecosystem.Type.ID==4),]

ggplot(zoo, aes(Soil.Pb,biom.totale))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
ggplot(zoo, aes(Vegetation.cover....,biom.totale))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
ggplot(zoo, aes(N.tot,average.masse))+geom_smooth(aes(color = Description), method = "loess",se=FALSE)
par(mfrow = c(1,1))
boxplot(zoo$average.masse~zoo$N.tot)
boxplot(zoo$biom.totale~zoo$N.tot)
boxplot(zoo$biom.nematodes~zoo$N.tot)
plot(fermes$Absolute.Re~fermes$N.tot)
plot(zoo$biom.totale~zoo$Manure.N)
plot(fermes$shannon.trophic~fermes$Manure.N)
plot(zoo$shannon.trophic~zoo$Total.N.input)
plot(fermes$shannon.trophic~fermes$Manure.N,col=fermes$Description)
ggplot(fermes, aes(Total.N.input,biom.totale))+geom_smooth(aes(color = Description),se=FALSE)
plot(fermes$biom.totale~fermes$N.tot,col=fermes$Description)
plot(fermes$abundance~fermes$Total.N.input,col=fermes$Description)
plot(zoo$Slope.of.AMR~zoo$Description)

plot(fermes$Taxa.S~fermes$Soil.phosphate)

mod.mixte<-lme(Taxa.S ~  Airborne.N + Total.N.input + Soil.Cd + 
                 Soil.Hg + Soil.Pb + Soil.Zn + Max.rainfall,data=zoo,random =~1|Description,method = "ML")
par(mfrow = c(2,2))
plot(mod.mixte)
summary(mod.mixte)
anova(mod.mixte)

ggplot(zoo, aes(Soil.phosphate,Taxa.S))+geom_smooth(aes(color = Description),se=FALSE)
