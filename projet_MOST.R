library("ppcor")
library(lsmeans)
#Pour tester l'effet sur un paramètre suelement
library(multcomp)
library(pROC)
library(plyr)
library(nlme)
library(ggplot2)
library(ncf)
library(lme4)
library(MASS)

#Importation

setwd("~/Documents/master2/most/Projet_2018-2019/")
foodweb<-read.table("food.csv",header=TRUE,sep=";")
zoo<-read.table("zoo.csv",header=TRUE,sep=";")

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

zoo$Latitude=jitter(zoo$Latitude)
zoo$Longitude=jitter(zoo$Longitude)
#Analyse descriptive de la composition des sols


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
ggplot(zoo, aes(Description, Soil.Zn ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilZn.png", width = 600)
ggplot(zoo, aes(Description, Soil.Cu ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilCu.png", width = 600)
ggplot(zoo, aes(Description, Soil.pH ))+geom_boxeplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilpH.png", width = 600)
ggplot(zoo, aes(Description, Soil.phosphate ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilphosphate.png", width = 600)

#Modèle Zn
mod.Zn<-gls(Soil.Zn~Description,data=zoo, method = "ML")
par(mfrow = c(2,2))
plot(mod.Zn)
summary(mod.Zn)
anova(mod.Zn)
Znpair=pairwise.t.test(zoo$Soil.Zn,zoo$Description,p.adjust.method = "bonferroni")
Znpair
lsmeans(mod.Zn,pairwise~Description,adjust="bonferroni")

cres = spline.correlog(x = zoo$Latitude, y = zoo$Longitude, z = resid(mod.Zn), resamp = 0)
plot(cres)

mod.Zn.rand<-gls(Soil.Zn~Description, data = zoo, correlation  = corSpatial(form = ~Latitude + Longitude, type ="exponential", nugget = T), method = "ML")
anova(mod.Zn,mod.Zn.rand)
lsmeans(mod.Zn.rand,pairwise~Description,adjust="bonferroni")


#Modèle pH
mod.pH<-lm(Soil.pH~Description,data=zoo)
par(mfrow = c(2,2))
plot(mod.)
summary(mod.pH)
anova(mod.pH)
pHpair=pairwise.t.test(zoo$Soil.pH,zoo$Description,p.adjust.method = "bonferroni")
pHpair
mod.phosphate<-lm(Soil.phosphate~Description,data=zoo)
par(mfrow = c(2,2))
plot(mod.phosphate)
summary(mod.phosphate)
anova(mod.phosphate)
phosphatepair=pairwise.t.test(zoo$Soil.phosphate,zoo$Description,p.adjust.method = "bonferroni")
phosphatepair

mod.Cu<-lm(Soil.Cu~Description,data=zoo)
par(mfrow = c(2,2))
plot(mod.Cu)
summary(mod.Cu)
anova(mod.Cu)
Cupair=pairwise.t.test(zoo$Soil.Cu,zoo$Description,p.adjust.method = "bonferroni")
Cupair

Cu<-lme(Soil.Cu~Description,data = zoo,random =~1|Average.T,method = "ML")
par(mfrow = c(2,2))
plot(Cu)
summary(Cu)
anova(Cu)
anova(mod.Cu,Cu)

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
dev.print(device = png, file = "figures/PCAall.png", width = 600)
fviz_pca_biplot(respcapro, axes = c(1,2), habillage = 'Description',select.var = list(contrib = 6),invisible = 'quali')

#corrplot(respcapro$var$cos2)
barplot(respcapro$var$contrib[,1],las=2)
barplot(respcapro$var$contrib[,2],las=2)
respcapro$var$contrib
fviz_pca_var(respcapro, axes = c(1,2), choix = 'var', select.var = list(contrib = 6))
dev.print(device = png, file = "figures/PCAallfeatures.png", width = 600)

zoo$PCA.dim.1<-respcapro$ind$coord[,1]

#On va essayer de faire la même chose mais qu'avec les fermes
fermes=zoo[(zoo$Ecosystem.Type.ID==1)|(zoo$Ecosystem.Type.ID==2)|(zoo$Ecosystem.Type.ID==3)|(zoo$Ecosystem.Type.ID==4),]
respcafermes = FactoMineR::PCA(X = fermes, # the data set used. Rows are individuals and columns are numeric variables.
                            scale.unit = TRUE,  #if TRUE, the data are scaled to unit variance before the analysis. 
                            quali.sup = c(1:29,34,35,48,49,50,51,53),
                            graph = F, #if TRUE, the graphs are displayed.
                            ncp = 5) #indexes of the annotation columns

eig <- get_eig(respcafermes)
fviz_screeplot(respcafermes, addlabels = TRUE)
ind <- get_pca_ind(respcafermes)
ind$coord
ind$contrib
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
dev.print(device = png, file = "figures/PCAfermes.png", width = 600)
fviz_pca_biplot(respcafermes, axes = c(1,2), habillage = 'Description',select.var = list(contrib = 6),invisible = 'quali')
#corrplot(respcapro$var$cos2)

#Faisons les tests sur les paramètres intéressants pour les fermes

mod.phosphate.fermes<-lm(Soil.phosphate~Description,data=fermes)
par(mfrow = c(2,2))
plot(mod.phosphate.fermes)
summary(mod.phosphate.fermes)
anova(mod.phosphate.fermes)
phosphatefermespair=pairwise.t.test(fermes$Soil.phosphate,fermes$Description,p.adjust.method = "bonferroni")
phosphatefermespair



mod.Ni.fermes<-lm(Soil.Ni~Description,data=fermes)
par(mfrow = c(2,2))
plot(mod.Ni.fermes)
summary(mod.Ni.fermes)
anova(mod.Ni.fermes)
Nifermespair=pairwise.t.test(fermes$Soil.Ni,fermes$Description,p.adjust.method = "bonferroni")
Nifermespair




barplot(respcafermes$var$contrib[,1],las=2)
barplot(respcafermes$var$contrib[,2],las=2)
respcafermes$var$contrib
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
barplot(respcafermes$var$contrib[,2],las=2)
fviz_pca_ind(respcafermes, axes = c(1,2), habillage = 'Description',invisible = 'quali',label='none')
fviz_pca_var(respcafermes, axes = c(1,2), choix = 'var', select.var = list(contrib = 6))
dev.print(device = png, file = "figures/PCAfermesfeatures.png", width = 600)

#Les boxplot des éléments chimiques du sol dans les fermes

ggplot(fermes, aes(Description, Soil.Ni ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilNifermes.png", width = 600)
ggplot(fermes, aes(Description, Soil.phosphate ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/soilPfermes.png", width = 600)
ggplot(fermes, aes(Description, P.pore.water ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))


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

####Section sur la biodiversité
zoor=zoo[zoo$Description!="forest",]
#Dans cette partie on enlève la forêt
##Influence avec le nombre de taxas
par(mfrow = c(1,1))
boxplot(zf$Taxa.S~zf$Description,las=2)

plot(zf$Soil.pH,zf$Taxa.S,las=2,col=zf$Description)
plot(zf$Soil.phosphate,zf$Taxa.S,las=2,col=zf$Description)

mod1<-lm(Taxa.S~Description,data=zoor)
mod2<-lme(Taxa.S~Description,random =~1|Average.T,data=zoor,method="ML")
mod5<-gls(Taxa.S~Soil.Cd,data=zoor,correlation = corAR1(0.5,~Average.T),method="ML")
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)
summary(mod1bis)
anova(mod1)
summary(mod2)
anova(mod2,mod1)
modpH<-glm(Taxa.S~Soil.pH,family = poisson(link = "log"),data = zoor)
summary(mopH)
anova(modpH)
mod1<-lm(Taxa.S~Description,data=zoor)
par(mfrow = c(2,2))
plot(mod1)
summary(mod1)
pairwise.t.test(zoor$Taxa.S,zoor$Description,p.adjust.method = "bonferroni")



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
    feed[[v[i]%/%10]]<-feed[[v[i]%/%10]]+exp(w[i])
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

ind<-c(95:138,140:173,175:194,196:232)
shannon.trophic.bis<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]$Trophic.ID
  w=zf[zf$Web.ID==j,]$Log.Abundance.
  feed<-list(0,0,0,0,0,0,0,0,0)
  somme=0
  sh=0
  for (i in 1:length(v)){
    feed[[v[i]%/%10]]<-feed[[v[i]%/%10]]+exp(w[i])
    somme=somme+exp(w[i])
  }
  for (k in 1:9){
    if (feed[[k]]!=0){
      r=feed[[k]]/somme
      sh=sh-r*log2(r)
    }
    
  }
  shannon.trophic.bis<-c(shannon.trophic.bis,sh)
}

ind<-c(95:138,140:173,175:194,196:232)
shannon.taxa<-c()
for (j in ind){
  v=zf[zf$Web.ID==j,]$Trophic.ID
  w=zf[zf$Web.ID==j,]$Log.Biomass.
  tax<-list(0,0,0,0,0)
  somme=0
  sh=0
  for (i in 1:length(v)){
    cat(v[i]%%10)
    tax[[v[i]%%10]]<-tax[[v[i]%%10]]+exp(w[i])
    somme=somme+exp(w[i])
  }
  for (k in 1:5){
    if (tax[[k]]!=0){
      r=tax[[k]]/somme
      sh=sh-r*log2(r)
    }
    
  }
  shannon.taxa<-c(shannon.taxa,sh)
}

ind<-c(95:138,140:173,175:194,196:232)
shannon.taxa.bis<-c()
for (j in ind){
  w=zf[zf$Web.ID==j,]$Log.Abundance.
  tax<-c()
  somme=0
  sh=0
  for (i in 1:length(w)){
    tax<-c(tax,exp(w[i]))
    somme=somme+exp(w[i])
  }
  for (k in 1:length(tax)){
    r=tax[k]/somme
    sh=sh-r*log2(r)
    
  }
  shannon.taxa.bis<-c(shannon.taxa.bis,sh)
}

ind<-c(95:138,140:173,175:194,196:232)
biomasse.tot<-c()
for (j in ind){
  w=zf[zf$Web.ID==j,]$Log.Biomass.
  biom=0
  for (i in 1:length(w)){
    biom=biom+exp(w[i])
  }
  biomasse.tot<-c(biomasse.tot,biom)
}


zoo["shannon.trophic.bis"]<-shannon.trophic
zoo["shannon.trophic"]<-shannon.trophic.bis
zoo["shannon.taxa"]<-shannon.taxa.bis
zoo["shannon.taxa.bis"]<-shannon.taxa
zoo["biomasse.tot"]<-biomasse.tot
zoor=zoo[zoo$Description!="forest",]


ggplot(zoo, aes(Description,shannon.trophic ))+ geom_dotplot(binaxis = "y",stackdir = "center",alpha=0.2,col=2)+geom_violin(alpha=0.2,col=2)
ggplot(zoor, aes(Description,shannon.trophic ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggplot(zoor, aes(Description,shannon.trophic.bis ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggplot(zoor, aes(Description,shannon.taxa ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggplot(zoor, aes(Description,shannon.taxa.bis ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggplot(zoor, aes(Description,biomasse.tot))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))

plot(zoo$shannon.trophic~zoo$PCA.dim.1,col=zoo$Description)
plot(zoo$Taxa.S~zoo$PCA.dim.1,col=zoo$Description)
ggplot(zoo, aes(PCA.dim.1,shannon.trophic ))+geom_smooth(aes(color = Description), model = lm)
ggplot(zoo, aes(PCA.dim.1,shannon.trophic ))+geom_smooth( model = lm)

#Les tests sans prendre en compte de facteurs aléatoires. 

modsoiltype<-glmmPQL(Taxa.S~Description,random=~1|grp,family = poisson(link = "log"),data = zoor)
modsoiltype0=glmmPQL(Taxa.S~1,random=~1|grp,family = poisson(link = "log"),data = zoor)
plot(modsoiltype)
summary(modsoiltype)
anova(modsoiltype,modsoiltype0)
lsmeans(modsoiltype,pairwise~Description,adjust="bonferroni")
ggplot(zoor, aes(Description,Taxa.S ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/boxtaxa.png", width = 600)
plot(zoor$Description,zoor$Taxa.S)

mod1chimrand<-glmmPQL(Taxa.S~Description,random=~1|grp,family = poisson(link = "log"),correlation= corSpatial(form = ~Latitude + Longitude, type ="exponential", nugget = T),data = zoor)

#Les tests sans facteurs aléatoires

soiltypeshannontaxa<-gls(shannon.taxa~Description,data = zoor,method = "ML")
par(mfrow=c(2,2))
plot(soiltypeshannontaxa)
summary(soiltypeshannontaxa)
anova(soiltypeshannontaxa)
lsmeans(soiltypeshannontaxa,pairwise~Description,adjust="bonferroni")
ggplot(zoor, aes(Description,shannon.taxa ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/boxtaxashannon.png", width = 600)

modsoiltypeshannontaxarand<-gls(shannon.taxa~Description, data = zoor, correlation  = corSpatial(form = ~Latitude + Longitude, type ="exponential", nugget = T), method = "ML")
anova(soiltypeshannontaxa,modsoiltypeshannontaxarand)

soiltypeshannontrophic<-lm(shannon.trophic~Description,data = zoor)
par(mfrow=c(2,2))
plot(soiltypeshannontrophic)
summary(soiltypeshannontrophic)
anova(soiltypeshannontrophic)
lsmeans(soiltypeshannontrophic,pairwise~Description,adjust="bonferroni")
ggplot(zoor, aes(Description,shannon.trophic ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/boxtrophicshannon.png", width = 600)

soiltypesbiom<-lm(biomasse.tot~Description,data = zoor)
par(mfrow=c(2,2))
plot(soiltypesbiom)
summary(soiltypesbiom)
anova(soiltypesbiom)
lsmeans(soiltypesbiom,pairwise~Description,adjust="bonferroni")
ggplot(zoor, aes(Description,biomasse.tot ))+geom_boxplot(notch=TRUE)+theme(axis.text.x = element_text(angle = 25, hjust = 1))
dev.print(device = png, file = "figures/boxbiom.png", width = 600)


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

#La biomasse en fonction de la composition chimique du sol et du type de sol

mod1chim<-glm(Taxa.S~Description+Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
               Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,family = poisson(link = "log"),data = zoor)

qqnorm(rstandard(mod1chim))
summary(mod1chim)
anova(mod1chim)



step(mod1chim,Taxa.S~Description+Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
       Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,direction="both")

mod1chimred<-glm(formula = Taxa.S ~ Description + Soil.Cu + Soil.Pb, family = poisson(link = "log"), data = zoor)

plot(mod1chimred)
summary(mod1chimred)
anova(mod1chimred)

par(mfrow=c(1,1))
plot(zoor$Soil.Cu,zoor$Taxa.S,col=zoor$Description)
interaction.plot(zoor$Soil.Cu,zoor$Description,zoor$Taxa.S,fixed=TRUE,col=2:3:6)
test<-gls(Taxa.S ~ Description , data=zoo, correlation =corSpher(form = ~ Latitude + Longitude, nugget = TRUE))

plot(zoo$Soil.Cu,zoo$Taxa.S,col=zoo$Description)

#Modèle de la richesse trophique en fonction de la composition du sol et de la nature du sol

mod2chim<-lm(shannon.trophic~Description+Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
                Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,family = poisson(link = "log"),data = zoor)
modnull<-lm(shannon.trophic~1,family = poisson(link = "log"),data = zoor)

qqnorm(rstandard(mod2chim))
summary(mod2chim)
anova(mod2chim)

step(mod2chim,shannon.trophic~Description+Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
       Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,direction="both")

mod2chimred<-glm(formula = Taxa.S ~ Description + Soil.Hg + Soil.Ptot + Soil.Zn, family = poisson(link = "log"), data = zoor)

#Modèle du nombre de taxons en fonction de la composition chimique du sol

mod3chim<-glm(Taxa.S~Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
                Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,family = poisson(link = "log"),data = zoor)
modnull<-glm(Taxa.S~1,family = poisson(link = "log"),data = zoor)

qqnorm(rstandard(mod3chim))
summary(mod3chim)
anova(mod3chim,modnull)



step(mod3chim,Taxa.S~Soil.pH+Soil.Cd+Soil.phosphate+Soil.Cr+Soil.Cu+Soil.Hg+Soil.Ni+Soil.Pb+
       Soil.Ptot+Soil.Zn+Max.rainfall+N.tot+Average.T+P.pore.water,direction="both")

mod3chimred<-glm(formula = Taxa.S ~ Description + Soil.Cu + Soil.Pb, family = poisson(link = "log"), data = zoor)

plot(mod3chimred)
summary(mod3chimred)
anova(mod3chimred)


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

ggplot(zoo, aes(Total.N.input, biomasse.tot))+geom_smooth(aes(color = Description), model = lm,se=FALSE)
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

#Distance

station.dists <- as.matrix(dist(cbind(zoo$Longitude, zoo$Latitude)))
station.dists.inv <- 1/(station.dists*10)
station.dists.inv[is.infinite(station.dists.inv)] <- 0

mod.mixte<-gls(Taxa.S ~ Description+Longitude,data=zoo,correlation = corAR1(0.5,~Longitude),method = "ML")

a=corExp(0.5,~dist(cbind(zoo$Longitude, zoo$Latitude)))
m<-lm(Taxa.S ~ Description, data = zoo)
cres = spline.correlog(x = zoo$Latitude, y = zoo$Longitude, z = resid(m), resamp = 0)
plot(cres)
zoo$grp = rep(1, dim(zoo)[1])
summary(m)         
