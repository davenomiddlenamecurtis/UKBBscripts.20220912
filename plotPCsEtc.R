wd="C:/Users/dave/OneDrive/sharedseq/UKBB"
ppi=600

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)

setwd(wd)
ethnicityCodes=data.frame(read.table("UKBB.ethnicityCodes.tsv",header=TRUE,sep="\t"))
# Use colors as here: https://www.researchgate.net/figure/Principal-components-analysis-A-above-and-spatial-ancestry-analysis-B-opposite-A_fig2_264390976
"
Brown - Africa
Blue - European
Yellow - East Asian
Green - South Asian
"

ethnicityCodes$colors=c(
"blue", # White
"darkblue", # British
"darkgoldenrod", # White and Black Caribbean
"green", # Indian
"brown", # Caribbean
"burlywood", # Mixed
"cyan", # Irish
"darkorange", # White and Black African
"limegreen", # Pakistani
"brown", # African
"palegreen2", # Asian or Asian British
"deepskyblue", # Any other white background
"lightgreen", # White and Asian
"olivedrab3", # Bangladeshi
"brown", # Any other Black background
"orange4", # Black or Black British
"orange", # Any other mixed background
"lightgreen", # Any other Asian background
"yellow", # Chinese
"magenta", # Other ethnic group
"white", # Do not know
"white" # Prefer not to answer
)

PCs=data.frame(read.table("ukb41465.exomes.allchr.eigenvec",header=TRUE,sep="\t"))
ethnicityTable=data.frame(read.table("UKBB.ethnicity.txt",header=TRUE,sep="\t"))
colnames(ethnicityTable)=c("IID","ethnicityCode")
ethnicityCodes$Meaning=as.character(ethnicityCodes$Meaning)
for (r in 1:nrow(ethnicityCodes)) {
  rows=ethnicityTable$ethnicityCode==ethnicityCodes$Coding[r]
  ethnicityTable$Ethnicity[rows]=ethnicityCodes$Meaning[r]
  }

ancestry=merge(PCs,ethnicityTable,by="IID")
scheme=scale_color_manual(
  breaks=ethnicityCodes$Meaning,
  values=ethnicityCodes$colors)

for (p in 1:3) {
  g=ggplot(ancestry,aes_q(x=ancestry[,p*2+1],y=ancestry[,p*2+2],color=factor(ancestry$Ethnicity)))+geom_point(size=1)+ theme_bw()
  g=g+xlab(colnames(ancestry)[p*2+1])+ylab(colnames(ancestry)[p*2+2])
  filename=sprintf("%s.v.%s.by.ethnicity.png",colnames(ancestry)[p*2+1],colnames(ancestry)[p*2+2])
#  png(filename,width=12*ppi, height=6*ppi, res=ppi)
#  print(g+scheme + guides(colour = guide_legend(override.aes = list(size=5)))+labs(colour="Ethnicity")) # this makes the key markers bigger
#  dev.off()
  }

centreTable=data.frame(read.table("UKBB.centre.txt",header=TRUE,sep="\t"))
countryOfBirthTable=data.frame(read.table("UKBB.countryOfBirth.txt",header=TRUE,sep="\t"))
sexTable=data.frame(read.table("UKBB.sex.txt",header=TRUE,sep="\t"))
northTable=data.frame(read.table("UKBB.north.txt",header=TRUE,sep="\t"))
eastTable=data.frame(read.table("UKBB.east.txt",header=TRUE,sep="\t"))

allPCs="PC1"
PCnames="PC1"
for (p in 2:20) { 
  allPCs=sprintf("%s + PC%d",allPCs,p)
  PCnames=append(PCnames,sprintf("PC%d",p))
  }
target="north"
fit=lm(as.formula(sprintf("toPlot$%s ~ %s",target,allPCs)),toPlot) 
summary (fit)
target="east"
fit=lm(as.formula(sprintf("toPlot$%s ~ %s",target,allPCs)),toPlot) 
summary (fit)
# get_tidy(fit,target)

townsendTable=data.frame(read.table("UKBB.townsend.txt",header=TRUE,sep="\t"))
YOBTable=data.frame(read.table("UKBB.YOB.txt",header=TRUE,sep="\t"))

ggerrorplot(toPlot,x="Ethnicity",y=PCnames,merge=TRUE)+rotate_x_text()

BMITable=na.omit(data.frame(read.table("UKBB.BMI.txt",header=TRUE,sep="\t")))
toPlot=merge(na.omit(ancestry),BMITable,by="IID")
g=ggerrorplot(toPlot,x="Ethnicity",y="BMI",merge=TRUE,palette="jco")+rotate_x_text()
png("BMI.v.ethnicity.png",width=12*ppi, height=6*ppi, res=ppi)
print(g)
dev.off()

BMIvEthnicity=data.frame(matrix(ncol=3,nrow=nrow(ethnicityCodes)))
colnames(BMIvEthnicity)=c("Ethnicity","Number","Mean")
ethnicityCodes$Meaning=as.character(ethnicityCodes$Meaning)
for (r in 1:nrow(ethnicityCodes)) {
  these=subset(toPlot,Ethnicity==ethnicityCodes$Meaning[r])
  BMIvEthnicity$Ethnicity[r]=ethnicityCodes$Meaning[r]
  BMIvEthnicity$Number[r]=nrow(these)
  BMIvEthnicity$Mean[r]=sprintf("%.2f (%.2f)",mean(these$BMI),sd(these$BMI)/sqrt(nrow(these)))
}
print(BMIvEthnicity)
write.table(BMIvEthnicity,"BMIvEthnicity.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

toPlot=merge(eastTable,PCs,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$east!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="east")
g=ggplot(dfm, aes(x=east,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("East coordinate")+ylab("PC value")
# png("PCs.v.east.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()

toPlot=merge(northTable,PCs,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$north!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="north")
g=ggplot(dfm, aes(x=north,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("North coordinate")+ylab("PC value")
# png("PCs.v.north.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()

toPlot=merge(townsendTable,PCs,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$townsend!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="townsend")
g=ggplot(dfm, aes(x=townsend,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("Townsend index")+ylab("PC value")
# png("PCs.v.townsend.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()

toPlot=merge(YOBTable,PCs,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$YOB!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="YOB")
g=ggplot(dfm, aes(x=YOB,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("Year of birth")+ylab("PC value")
# png("PCs.v.YOB.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()

toPlot=merge(BMITable,PCs,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$BMI!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
dfm=melt(toPlot,id.vars="BMI")
g=ggplot(dfm, aes(x=BMI,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("BMI")+ylab("PC value")
# png("PCs.v.BMI.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()


toPlot=merge(BMITable,YOBTable,by="IID")
toPlot=merge(toPlot,northTable,by="IID")
toPlot=merge(toPlot,eastTable,by="IID")
toPlot=merge(toPlot,townsendTable,by="IID")
toPlot=toPlot[complete.cases(toPlot),]
toPlot=toPlot[toPlot$east!=0,]
toPlot=toPlot[toPlot$north!=0,]
drops <- c("IID","FID")
toPlot=toPlot[,!(names(toPlot) %in% drops)]
for (c in 2:5) { toPlot[,c]=(toPlot[,c]-mean(toPlot[,c]))/sd(toPlot[,c])}
dfm=melt(toPlot,id.vars="BMI")
g=ggplot(dfm, aes(x=BMI,y=value, color=variable)) +
  geom_smooth(se=FALSE) 
g=g+xlab("BMI")+ylab("Normalised values")
# png("BMI.v.others.png",width=12*ppi, height=6*ppi, res=ppi)
# print(g)
# dev.off()

RAB35=data.frame(read.table("UKBB.BMI.all.RAB35.sco",header=FALSE))
colnames(RAB35)=c("IID","BMI2","score")
toPlot=merge(toPlot,RAB35,by="IID")
g=ggplot(toPlot,aes_q(x=toPlot$BMI,y=toPlot$score,color=factor(toPlot$Ethnicity)))+geom_point(size=1)+ theme_bw() 
g+scheme + guides(colour = guide_legend(override.aes = list(size=5)))
fit =aov(y~ethnicity, data=ancestry)

