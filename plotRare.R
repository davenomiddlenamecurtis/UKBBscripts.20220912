# show geographical distribution of rare variants

wd="C:/Users/dave/OneDrive/sharedseq/UKBB/rare"
ppi=600

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)

# just some varied colours
colours=c(
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
"magenta" # Other ethnic group
)

setwd(wd)
rareVarsTable=data.frame(read.table("UKBB.rare.22.first100.txt",header=TRUE,sep=" "))
northTable=data.frame(read.table("../UKBB.north.txt",header=TRUE,sep="\t"))
eastTable=data.frame(read.table("../UKBB.east.txt",header=TRUE,sep="\t"))

nVars=length(colours)
varNames=vector(mode="character", length=nVars)

all=merge(na.omit(northTable[northTable$north!=0,]),rareVarsTable,by="IID")
all=merge(na.omit(eastTable[eastTable$east!=0,]),all,by="IID")
offset=0
toPlot=data.frame(IID=numeric(),east=numeric(),north=numeric(),var=numeric())
for (v in 1:nVars)
{
  toAdd=na.omit(all[all[,v+7+offset]!=0,(1:3)])
  varNames[v]=colnames(all)[v+7+offset]
  toAdd$var=colnames(all)[v+7+offset]
  toPlot=rbind(toPlot,toAdd)
  print(nrow(toAdd))
}

scheme=scale_color_manual(breaks=vars,values=colours)
g=ggplot(toPlot,aes_q(x=toPlot$east,y=toPlot$north,color=factor(toPlot$var)))+geom_point(size=1)+ theme_bw()
g=g+xlab("east")+ylab("north")

g


