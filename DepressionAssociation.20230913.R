#!/share/apps/R-3.6.1/bin/Rscript

# Case control study of 270K exome sample for 
# 2090-0.0	Seen doctor (GP) for nerves, anxiety, tension or depression
# Cases being HERC1 LOF carriers

CarrierListFile="HERC1.LOF.carriers.new.20230906.txt"
NewIDListFile="/SAN/ugi/UGIbiobank/data/downloaded/UKBB.newIn2023.exome.IDs.txt"
SeenGPFile="UKBB.SawGPNerves.txt"

Results=data.frame(matrix(ncol=2,nrow=2))
rownames(Results)=c("No","Yes")
colnames(Results)=c("NotLOFCarrier","LOFCarrier")
PResults=data.frame(matrix(ncol=2,nrow=2))
rownames(PResults)=c("No","Yes")
colnames(PResults)=c("NotLOFCarrier","LOFCarrier")
Table=data.frame(matrix(ncol=2,nrow=2))
rownames(Table)=c("No","Yes")
colnames(Table)=c("NotLOFCarrier","LOFCarrier")


NewIDs=data.frame(read.table(NewIDListFile,header=FALSE,stringsAsFactors=FALSE))
colnames(NewIDs)=c("IID")

Carriers=data.frame(read.table(CarrierListFile,header=FALSE,stringsAsFactors=FALSE))
colnames(Carriers)=c("IID")

SeenGP=data.frame(read.table(SeenGPFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

SeenGP=merge(SeenGP,NewIDS,by="IID",all=FALSE)
Results[1,1]=sum(SeenGP$SawGPNerves==0,na.rm=TRUE)
Results[2,1]=sum(SeenGP$SawGPNerves==1,na.rm=TRUE)
Carriers=merge(Carriers,SeenGP,by="IID",all=FALSE)
Results[1,2]=sum(Carriers$SawGPNerves==0,na.rm=TRUE)
Results[2,2]=sum(Carriers$SawGPNerves==1,na.rm=TRUE)
Results[1,1]=Results[1,1]-Results[1,2]
Results[2,1]=Results[2,1]-Results[2,2]
print(Results)

for (c in 1:2) {
	PResults[1,c]=Results[1,c]/(Results[1,c]+Results[2,c])*100
	PResults[2,c]=100-PResults[1,c]
	for (r in 1:2) {
	Table[r,c]=sprintf("%d (%.1f%%)",Results[r,c],PResults[r,c])
	}
}

print(Table)



