
##################
#Setting up
####################

ProjectSTRs=list()
#ProjectSTRs$ProjectRoot="~/..."

library(BioGeoBEARS)
library(foreach)
library(RMySQL)
library(doParallel)
library(survival)
library(ggplot2)
library(plyr)
library(stringr)
library(gplots)
library(devtools)
library(RColorBrewer)
library(R.matlab)
library(limma)
library(devtools)
library(pROC)
library(abind)
library(calibrate)
library(Rmisc)
library(pamr)
library(RCurl)
library(BioGeoBEARS)
library(HGNChelper)


#RscriptsPath="path_to_rscripts/"
#source(paste(RscriptsPath,"script",sep=""))

pcname <- system('uname -n',intern=T)
ProjectSTRs=list()
ProjectSTRs$ProjectRoot="~/Documents/Projects/TCGA_HNSC/"

######
#
#########

#Load DNA methylation and gene expression array data 
#METcancer refers to preprocessed DNA methylation data for tumor samples only (CpG in rows, patients in columns)
#METnormal refers to preprocessed DNA methylation data for normal adjacent tissue samples only (CpG in rows, patients in columns)
#MA cancer refers to gene expression (microarray or RNA-Seq) data for tumor samples only (Genes in rows, patients in columns)

#######################
#Script specific to NSD1 analysis 
#############################

#Use hierarchical clustering to cluster individual CpG sites into CpG cluster
Clusters_MET_KIRP= TCGA_GENERIC_MET_ClusterProbes_with_hclust(METcancer,METnormal,Parallel = T)

Clusters_METcancer <- Clusters_MET_KIRP[[1]]
Clusters_METnormal <- Clusters_MET_KIRP[[2]]
#Clean up the sample barcodes names
Clusters_METcancer <- TCGA_GENERIC_CleanUpSampleNames(Clusters_METcancer, 12)
Clusters_METnormal <- TCGA_GENERIC_CleanUpSampleNames(Clusters_METnormal, 12)

#Get the data that maps CpG probes to CpG clusters
ProbeMapping=Clusters_MET_KIRP[[3]]
colnames(ProbeMapping)=c("probe","MethylMixGene")
ProbeMapping=as.data.frame(ProbeMapping)

####################
#Run MethylMix to get abnormally methylated genes (Same for both HNSC and LUSC)
#####################

#Run MethylMix just for HNSC
library(MethylMix)
HNSC.MethylMix <- MethylMix("univ.beta", Clusters_METcancer, Clusters_METnormal, MAcancer, Parallel = T)
#Get methylation state data (The DM value matrix)
MethylationStates=HNSC.MethylMix$MethylationStates

################
#Run MethylMix and consensus clustering for all cancers with 10+ NSD1 mutations
########################
cancer=c("HNSC", "UCEC", "STAD","LUSC")

pathMeth='path to DNA methylation data'
pathExp='path to gene expression data'

for(i in 1:length(cancers)){
  cancer=cancers[i]
  load(paste("path/MET_ProcessedData_450k.RData", sep=""))
  METcancer=ProcessedData$MET_Data_Cancer
  METnormal=ProcessedData$MET_Data_Normal
  
  #Make CpG clusters
  Clusters_MET=TCGA_GENERIC_MET_ClusterProbes_with_hclust(METcancer,METnormal,Parallel = T)
  Clusters_METcancer <- Clusters_MET[[1]]
  Clusters_METnormal <- Clusters_MET[[2]]
  ProbeMapping=Clusters_MET[[3]]
  
  Clusters_MET_Results=list()
  Clusters_MET_Results[[1]]=Clusters_METcancer
  Clusters_MET_Results[[2]]=Clusters_METnormal
  Clusters_MET_Results[[3]]=ProbeMapping
  names(Clusters_MET_Results)=c("Clusters_MET_Results", "Clusters_MET_Results", "ProbeMapping")
  save(Clusters_MET_Results, file=paste("Directory",cancer[i], "/data/","Clusters_MET.RData", sep=""))
}

#Run MethylMix for NSD1 related cancers

for(i in 1:length(cancer)){
  #load Methylation cluster data
  cancer=cancer[i]
  dat=load(paste("Directory",cancer, "/data/","Clusters_MET.RData", sep=""))
  names(Clusters_MET_Results)=c("Clusters_MET_Cancer","Clusters_MET_Normal","ProbeMapping")
  Clusters_METcancer=Clusters_MET_Results$Clusters_MET_Cancer
  Clusters_METnormal=Clusters_MET_Results$Clusters_MET_Normal
  colnames(Clusters_METcancer)=substr(colnames(Clusters_METcancer),1,12)
  colnames(Clusters_METnormal)=substr(colnames(Clusters_METnormal),1,12)
  
  dat=load(paste("Directory",cancer, "/data/","MA_",cancer,"_Processed.RData", sep=""))
  MAcancer=ProcessedData$MA_Data_Cancer
  colnames(MAcancer)=substr(colnames(MAcancer),1,12)
  
  #Do MethylMix
  Results.MethylMix <- MethylMix("univ.beta", Clusters_METcancer, Clusters_METnormal, MAcancer, Parallel = T)
  
  save(Results.MethylMix, file=paste("Directory",cancer, "/data/Results.MethylMix_450k.RData", sep=""))
}

##################################
#For each cancer, get number of abnormally (hypomethylated or hypermethylated) genes
#################################

NSD1_HypoN=array(NA,c(length(cancer), 2))
rownames(NSD1_HypoN)=cancer
colnames(NSD1_HypoN)=c('N_HypoGene_Diff',"P.value")

NSD1_HyperN=array(NA,c(length(cancer), 2))
rownames(NSD1_HyperN)=cancer
colnames(NSD1_HyperN)=c('N_HypoGene_Diff',"P.value")

i=1
for(i in 1:length(cancers)){
  dat=load(paste("~/Documents/Projects/TCGA_",cancer, "/data/Results.MethylMix_450k.RData", sep=""))
  MethylationStates=Results.MethylMix$MethylationStates
  
  HyperMethGenesAll=list()
  for(i in 1:ncol(MethylationStates)){
    HyperMethGenesAll[[i]]=gsub("---.*","", names(MethylationStates[,i][MethylationStates[,i]>0]))
  }
  names(HyperMethGenesAll)=colnames(MethylationStates)
  
  HypoMethGenesAll=list()
  for(i in 1:ncol(MethylationStates)){
    HypoMethGenesAll[[i]]=gsub("---.*","", names(MethylationStates[,i][MethylationStates[,i]<0]))
  }
  names(HypoMethGenesAll)=colnames(MethylationStates)
  
  #
  Abnormal_MethylMix_state_N=array(NA, c(length(HyperMethGenesAll),2))
  rownames(Abnormal_MethylMix_state_N)=names(HyperMethGenesAll)
  colnames(Abnormal_MethylMix_state_N)=c("N_Hypermethylated","N_Hypomethylated")
  
  for(i in 1:length(HyperMethGenesAll)){
    Abnormal_MethylMix_state_N[i,1]=length(HyperMethGenesAll[[i]])
    Abnormal_MethylMix_state_N[i,2]=length(HypoMethGenesAll[[i]])
  }
  Abnormal_MethylMix_state_N=as.data.frame(Abnormal_MethylMix_state_N)
  
  
#################
#Run consensus clustering on the DM value matrix for each NSD1 related cancer
#######################

ConsensusClusterPlus(MethylationStates,maxK=10,reps=1000,pItem=0.8,pFeature=as.numeric(1),clusterAlg='km',distance='euclidean',verbose=TRUE,title='MethylMix_HNSC_ConsensusClustering', writeTable=TRUE, plot="png")

#Get the cluster labels for sample with 5 clusters

cluster5=read.csv("path..to..ConsensusClusterPlus..output/MethylMix_HNSC_ConsensusClustering.k=5.consensusClass.csv", header=F)
rownames(cluster5)=cluster5[,1]
colnames(cluster5)=c("patient_barcode","cluster")

#########################
#Make heatmap of DM valuesfor HNSC, as in figure 1a
##############################

#using heatmap 3, an amendment to heatmap2
library("gplots")
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#order samples by DNA methylation clusters
cluster5_order=cluster5[order(cluster5$cluster),]
MethylationStates_clust5=MethylationStates[,match(cluster5_order$patient_barcode, colnames(MethylationStates))]
MethylationStates_clust5=t(MethylationStates_clust5)
Mcolors=as.matrix(revalue(as.factor(cluster5_order[rownames(cluster5_order) %in% rownames(MethylationStates_clust5),][,"cluster"]),
                          c("1"="lightgrey", "2"="red", "3"="darkgrey", "4"="lightgrey","5"="darkgrey")))
#Make column labels
clab=cbind(Mcolors)
colnames(clab)=c("")

geneslab=matrix('white', length(colnames(MethylationStates_clust5)),1)
geneslab[colnames(MethylationStates_clust5) %in% geneshypo]="black"
rlab=t(geneslab) 
BuRd=c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7" ,"#FDDBC7", "#F4A582", "#D6604D", "#B2182B")    

file=paste("path/filename",sep="")
png(file=paste(file,'.png',sep=''), units="in", width=11, height=6, res=600)
heatmap.3(t(MethylationStates_clust5),Rowv=TRUE, Colv=FALSE,ColSideColors=clab, RowSideColors=rlab,dendrogram = "none",ColSideColorsSize=4,RowSideColorsSize=2, col=BuRd, cexRow=2, labRow ="",labCol = "")
dev.off()

#Make equivalent heatmap for LUSC
ov=overlap(MethylationStates, cluster6)
MethylationStates_clust6=ov[[1]]
cluster6_MethylationStates=ov[[2]]

cluster6_order=cluster6_MethylationStates[order(cluster6_MethylationStates$cluster),]

MethylationStates_clust6=MethylationStates_clust6[,match(cluster6_order$patient_barcode, colnames(MethylationStates_clust6))]
MethylationStates_clust6=t(MethylationStates_clust6)

Mcolors=as.matrix(revalue(as.factor(cluster6_order[rownames(cluster6_order) %in% rownames(MethylationStates_clust6),][,"cluster"]),
                          c("1"="lightgrey", "2"="darkgrey", "3"="red", "4"="lightgrey","5"="darkgrey", "6"="lightgrey")))
clab=cbind(Mcolors)
colnames(clab)=c("")

geneslab=matrix('white', length(colnames(MethylationStates_clust6)),1)
geneslab[colnames(MethylationStates_clust6) %in% geneshypo]="black"
rlab=t(geneslab) 

file=paste("path/file",sep="")
png(file=paste(file,'.png',sep=''), units="in", width=11, height=6, res=600)
heatmap.3(t(MethylationStates_clust6),Rowv=TRUE, Colv=FALSE, dendrogram = "none",RowSideColors=rlab, ColSideColors=clab, RowSideColorsSize=4, ColSideColorsSize=4, col=BuRd,labRow = "",labCol = "")
dev.off()

#########################
#Get and process MutSig.2 point mutationn data
#############################

QvalueThreshold=0.01
path=... #specify path
MutSigDirectory=paste(path,'gdac.broadinstitute.org_HNSC-TP.MutSigNozzleReport2.0.Level_4.2014041600.0.0/',sep="")
MutSig_TCGA_HNSC=read.table(paste(MutSigDirectory,"HNSC-TP.sig_genes.txt", sep=''), header=T, sep="\t")[1:43,2]

#remove silent mutations
silent<- which(with( MUT_TCGA_HNSC, (MUT_TCGA_HNSC$is_silent=="1")))
MUT_TCGA_HNSC=MUT_TCGA_HNSC[-silent,] 
MUT_TCGA_HNSC$bcr_patient_barcode=substr(MUT_TCGA_HNSC$Tumor_Sample_Barcode, 1,12)

#Make a matrix of the significantly mutated genes, indicating number of mutations for each patient
epigenes=MutSig_TCGA_HNSC
patients=unique(MUT_TCGA_HNSC$bcr_patient_barcode)

MutSigMatrix=array(NA, c(length(epigenes),length(patients)))
rownames(MutSigMatrix)=epigenes
colnames(MutSigMatrix)=patients

for(j in 1:length(epigenes)){
  Patients_with_Mutation=MUT_TCGA_HNSC[MUT_TCGA_HNSC$Hugo_Symbol==as.character(epigenes[j]),"bcr_patient_barcode"]
  Patients_with_Mutation=Patients_with_Mutation[!is.na(Patients_with_Mutation)]
  for(i in 1:length(patients)){
    MutSigMatrix[j,i]=length(Patients_with_Mutation[Patients_with_Mutation==patients[i]])
  }
}

#Make binary version (mutated/not mutated)
MutSigMatrix_binary=MutSigMatrix
for(i in 1:nrow(MutSigMatrix_binary)){
  MutSigMatrix_binary[i,MutSigMatrix_binary[i,]!=0]=1  
}

###############
#Intersect preprocessed DNA methylation array data with DNA methylation clusters for HNSC and LUSC
#################

#HNSC
ov=overlap(METcancer, cluster5)
METcancer_cluster5=ov[[1]]
cluster5_METcancer=ov[[2]]

#get probes in MethylMix genes
HNSC_NSD1_AllProbes=ProbeMappingHNSC[ProbeMappingHNSC$MethylMixGene %in% rownames(MethylationStates),"probe"]
#10,818 CpGs
METcancer_cluster5_MethylMixGenes=METcancer_cluster5[HNSC_NSD1_AllProbes,]

#scale by row/gene
METcancer_cluster5_MethylMixGenes=t(apply(METcancer_cluster5_MethylMixGenes,1, scale))
colnames(METcancer_cluster5_MethylMixGenes)=colnames(METcancer_cluster5)

#LUSC
ov=overlap(METcancerLUSC, cluster6)
METcancerLUSC_cluster6=ov[[1]]
cluster6_METcancerLUSC=ov[[2]]

#Restrict HNSC and LUSC data to probes that overlap between two cancer types
overlapProbes=intersect(rownames(METcancer_cluster5_MethylMixGenes), rownames(METcancerLUSC_cluster6))
METcancer_cluster5_MethylMixGenes=METcancer_cluster5_MethylMixGenes[overlapProbes,]
METcancerLUSC_cluster6_MethylMixGenes=METcancerLUSC_cluster6[overlapProbes,]

#Scale LUSC data
METcancerLUSC_cluster6_MethylMixGenes=t(apply(METcancerLUSC_cluster6_MethylMixGenes,1, scale))
colnames(METcancerLUSC_cluster6_MethylMixGenes)=colnames(METcancerLUSC_cluster6)

#Make heatmap of correlation matrix between HNSC and LUSC
pwcarray=array(NA, c(ncol(METcancerLUSC_cluster6_MethylMixGenes),ncol(METcancer_cluster5_MethylMixGenes)))
rownames(pwcarray)=colnames(METcancerLUSC_cluster6_MethylMixGenes)
colnames(pwcarray)=colnames(METcancer_cluster5_MethylMixGenes)

for(i in 1:ncol(METcancerLUSC_cluster6_MethylMixGenes)){
  for(j in 1:ncol(METcancer_cluster5_MethylMixGenes)){
    pwcarray[i,j]=cor(METcancerLUSC_cluster6_MethylMixGenes[,i],METcancer_cluster5_MethylMixGenes[,j], method="pearson")
  }
}

cluster6_METcancerLUSC_order=cluster6_METcancerLUSC[order(cluster6_METcancerLUSC$cluster),]
cluster5_METcancer_order=cluster5_METcancer[order(cluster5_METcancer$cluster),]

pwcarray2=pwcarray[match(cluster6_METcancerLUSC_order$patient_barcode, rownames(pwcarray)),match(cluster5_METcancer_order$patient_barcode, colnames(pwcarray))]

#Make heatmap of the correlation matrix 
file=paste('directory/filename')
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
heatmap.3(pwcarray2,Rowv=FALSE, Colv=FALSE, dendrogram = "none", col=BuRd, labRow = "", RowSideColors=t(rlab), labCol = "", ColSideColors=clab,  ColSideColorsSize=6,  RowSideColorsSize=6)
dev.off()

#########################
#Get NSD1 copy number and intersect with NSD1 point mutations
#################

GisticDirectory="/Users/kbrennan/Documents/Projects/TCGA_HNSC/data/2015-05-21/GISTIC/"
gistic=TCGA_Load_GISTICdata(GisticDirectory)
#457 AMP genes and 3112 DEL genes
CGH_Data_Thresholded=gistic$CGH_Data_Thresholded

#Intersect mutation and copy number data 
OverlapSamples=intersect(colnames(CGH_Data_Thresholded), colnames(MutSigMatrix))
MutSigMatrix_CGH=MutSigMatrix[,OverlapSamples]
CGH_MutSigMatrix=CGH_Data_Thresholded[,OverlapSamples]

#286 overlapping samples
NSD1_MutCGH=merge(MutSigMatrix_CGH["NSD1",],
                  CGH_MutSigMatrix["NSD1",], by=0)
rownames(NSD1_MutCGH)=NSD1_MutCGH[,1]
NSD1_MutCGH=NSD1_MutCGH[,2:ncol(NSD1_MutCGH)]
colnames(NSD1_MutCGH)=c("Mutation", "CNV")

#Make plots indicating percentage of mutations and deletions in each DNA methylation subtype
NSD1_CGH=CGH_Data_Thresholded["NSD1",]

cluster5_NSD1MUTCNV=cluster5

cluster5_NSD1MUTCNV$CNV=rep(NA, nrow(cluster5_NSD1MUTCNV))
cluster5_NSD1MUTCNV$Mut=rep(NA, nrow(cluster5_NSD1MUTCNV))

for(i in 1:length(cluster5_NSD1MUTCNV$CNV)){
  cluster5_NSD1MUTCNV$CNV[i]=ifelse(length(intersect(names(NSD1_CGH), rownames(cluster5_NSD1MUTCNV)[i]))!=0, as.character(NSD1_CGH[rownames(cluster5_NSD1MUTCNV)[i]]), NA) 
} 

NSD1_Mut=NSD1_Mutations_Cancers$HNSC$Samples_with_NonSilent_NSD1_barcodes

#NSD1_Mut=MutSig_TCGA_HNSC$Samples_with_NonSilent_NSD1_barcodes
names(NSD1_Mut)=rownames(NSD1_Mutations_Cancers$HNSC)

for(i in 1:length(cluster5_NSD1MUTCNV$Mut)){
  cluster5_NSD1MUTCNV$Mut[i]=ifelse(length(intersect(names(NSD1_Mut), rownames(cluster5_NSD1MUTCNV)[i]))!=0, as.character(NSD1_Mut[rownames(cluster5_NSD1MUTCNV)[i]]), NA) 
} 
cluster5_NSD1MUTCNV$dels=revalue(as.factor(cluster5_NSD1MUTCNV$CNV), c("0"="0","1"="0","-1"="1","-2"="1","2"="0"))

tab=array(NA,c(2,5))
rownames(tab)=c("Mut","Del")
colnames(tab)=c(1:5)

for(i in 1:5){
  tab["Mut",i]=length(which(cluster5_NSD1MUTCNV[cluster5_NSD1MUTCNV$cluster==i & !is.na(cluster5_NSD1MUTCNV$Mut),"Mut"]==1))/length(cluster5_NSD1MUTCNV[cluster5_NSD1MUTCNV$cluster==i & !is.na(cluster5_NSD1MUTCNV$Mut),"Mut"])
  tab["Del",i]=length(which(cluster5_NSD1MUTCNV[cluster5_NSD1MUTCNV$cluster==i & !is.na(cluster5_NSD1MUTCNV$dels),"dels"]==1))/length(cluster5_NSD1MUTCNV[cluster5_NSD1MUTCNV$cluster==i & !is.na(cluster5_NSD1MUTCNV$dels),"dels"])
}

#file=paste('Directory/file')
#png(file=paste(file,'.png',sep=''), units="in", width=5, height=3, res=600)
barplot(tab, beside=T, col=c("grey","grey","red","red","grey","grey","grey","grey","grey","grey"), density=c(10,100), cex.axis=1, xaxt="n")
#dev.off()

###################
#Overlap between cancer signature and Sotos syndrome (hypomethylation signature)
#####################

#Get downloaded list of hypomethylated CpGs from the Choufani et al paper
sotos=read.table("Directory/Sotos_probes.txt", header=T, sep="\t", fill=T, quote="", comment.char="")
sotosHypo=sotos[grep("Loss",sotos$DNA.methylation.effect),"Illumina.ID"]

#For each patient, get a list of genes that are hypomethylated in cancer, from the DM value matrix (MethylationStates)
HypoMethGenes_Cluster=list()
for(i in 1:ncol(MethylationStates)){
  HypoMethGenes_Cluster[[i]]=names(MethylationStates[,i][MethylationStates[,i]<0])
}
names(HypoMethGenes_Cluster)=colnames(MethylationStates)

#Get the list of CpGs within each gene 
Patient_HypocgsLUSC=list()
for(i in 1:length(HypoMethGenes_Cluster)){
  Patient_HypocgsLUSC[[i]]=ProbeMappingLUSC[ProbeMappingLUSC[,2] %in% HypoMethGenes_Cluster[[i]],1]
}
names(Patient_HypocgsLUSC)=colnames(MethylationStates)

#Make Sotos syndrome overlap index and test association with NSD1 mutations or deletions
GeneSig=sotosHypo
#get the CpGs that are measured on the array
arraygenes=rownames(METcancer)
GeneSample=Patient_Hypocgs

tab=make_hypergeometic(GeneSample=GeneSample, GeneSig=GeneSig, arraygenes=arraygenes)
hyper_1=tab
tab=hyper_1

#Make boxplot displaying levels of Sotos syndrome overlap index in each DNA methylation subtype
p1=make_ggplot2_boxplot(tab1=tab, tab2=cluster5)

file="Directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p1[[1]] + geom_boxplot(aes(fill = factor(cluster)), outlier.shape=NA) + geom_jitter() + theme(axis.ticks = element_blank(), legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.20.text) + scale_fill_manual(name = "NSD1\nmutations (n)", values = c("grey", "red","grey","grey","grey"), labels = c("1" = "1", "2" = "2","3" = "3","4" = "4","5" = "5")) 
dev.off()

#Make boxplot displaying levels of Sotos syndrome overlap index with the number of mutations and deletions 
p1=makep_mutdel(tab1=tab, tab2=NSD1_MutCGH)

file="Directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p1 + geom_boxplot( aes(fill =  factor(Mut)),outlier.shape=NA, position=position_dodge(0.8), width=0.8)+ scale_x_discrete("CNV", breaks=factor(1:4), drop=FALSE) + facet_grid(. ~ CNV, labeller="label_both") + geom_jitter() + theme(legend.position='none', axis.text.y = black.bold.20.text, axis.title.x=element_blank(), legend.text=element_text(size=25, face="bold"),legend.title=element_text(size=25, face="bold"), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.justification=c(.5,.5)) + scale_fill_manual(name = "NSD1\nmutations (n)", values = c("yellow", "orange","red"), labels = c("0" = "0", "1" = "1","2"="2")) + theme(strip.text.x = element_text(size = 30, colour = "black", face="bold"))
dev.off()

ov=overlap(tab, NSD1_MutCGH)
tab_NSD1_MutCGH=ov[[1]]
NSD1_MutCGH_tab=ov[[2]]

##Test association between Sotos styndrome overlap index and the NSD1 lesion score using linear regression 
summary(glm(tab_NSD1_MutCGH$Percentage_Overlap~NSD1_MutCGH_tab$nlesions))
#p=1.09e-12

#Testing difference in levels of Sotos syndrome overlap index between the NSD1 subtype and each other subtype
nlesions=levels(as.factor(NSD1_MutCGH_tab$nlesions))
for(i in 1:length(nlesions)){
  print(wilcox.test(tab_NSD1_MutCGH$Percentage_Overlap[NSD1_MutCGH_tab$nlesions==0],
                    tab_NSD1_MutCGH$Percentage_Overlap[NSD1_MutCGH_tab$nlesions==nlesions[i]]))
}

t1=cbind(tab_NSD1_MutCGH$Percentage_Overlap, NSD1_MutCGH_tab$nlesions)
colnames(t1)=c("Meth","Mut")
t1=as.data.frame(t1)
p=ggplot(t1, aes(factor(Mut), Meth))  

file="Directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p + geom_boxplot( aes(fill =  factor(Mut)),outlier.shape=NA, position=position_dodge(0.8), width=0.8)+ scale_x_discrete("CNV", breaks=factor(1:4), drop=FALSE) + geom_jitter() + theme(legend.position='none', axis.text.y = black.bold.20.text, axis.title.x=element_blank(), legend.text=element_text(size=25, face="bold"),legend.title=element_text(size=25, face="bold"), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.justification=c(.5,.5)) + scale_fill_manual(name = "NSD1\nmutations (n)", values = cols, labels = c("-1"="-1","0" = "0", "1" = "1","2"="2")) 
dev.off()

#make a random list of CpGs and calculate overlap between the random signature and each cancer patient's hypomethylated CpG lies
#Perform this ten times and make a list of arrays
Patient_Hypocgs=Patient_HypocgsLUSC

HypoSampleHNSCRRandom=list()
samplength=min(unlist(lapply(Patient_Hypocgs, length)))

GeneSig=sample(rownames(METcancerHNSC), length(sotosHypo))

for(i in 1:10){
  GeneSample=lapply(Patient_Hypocgs, function(x) sample(x, samplength))
  tab=make_hypergeometic(GeneSample=GeneSample, GeneSig=GeneSig, arraygenes=arraygenes)
  HypoSampleHNSCRRandom[[i]]=tab
}
hyperlist=HypoSampleHNSCRRandom

#Take the mean overlap acorss the ten interations 
for(i in 1:ncol(hyperlist[[1]])){
  hypermean[,i]=rowMeans(cbind(hyperlist[[1]][,i],
                               hyperlist[[2]][,i],
                               hyperlist[[3]][,i],
                               hyperlist[[4]][,i],
                               hyperlist[[5]][,i],
                               hyperlist[[6]][,i],
                               hyperlist[[7]][,i],
                               hyperlist[[8]][,i],
                               hyperlist[[9]][,i],
                               hyperlist[[10]][,i]))
}
hypermean=as.data.frame(hypermean)

#Make boxplot for association of random overlap with methylation subtypes
p1=make_ggplot2_boxplot(tab1=hypermean, tab2=cluster5)

#Make boxplot using ggplot2
file="directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p1[[1]] + geom_boxplot(aes(fill = factor(cluster)), outlier.shape=NA) + geom_jitter() + theme(axis.ticks = element_blank(), legend.position="none", axis.text.x = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(), plot.title = element_text(lineheight=.8, face="bold", size=20), axis.title = element_text(lineheight=.8, face="bold", size=16)) + theme(axis.text.y = black.bold.20.text) + scale_fill_manual(name = "NSD1\nmutations (n)", values = c("grey", "red","grey","grey","grey"), labels = c("1" = "1", "2" = "2","3" = "3","4" = "4","5" = "5")) 
dev.off()

#Make boxplot for association of random overlap with NSD1 mutations or deletions
p1=makep_mutdel(tab1=hypermean, tab2=NSD1_MutCGH)

file="~/Directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p1 + geom_boxplot( aes(fill =  factor(Mut)),outlier.shape=NA, position=position_dodge(0.8), width=0.8)+ scale_x_discrete("CNV", breaks=factor(1:4), drop=FALSE) + facet_grid(. ~ CNV, labeller="label_both") + geom_jitter() + theme(legend.position='none', axis.text.y = black.bold.20.text, axis.title.x=element_blank(), legend.text=element_text(size=25, face="bold"),legend.title=element_text(size=25, face="bold"), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.justification=c(.5,.5)) + scale_fill_manual(name = "NSD1\nmutations (n)", values = c("yellow", "orange","red"), labels = c("0" = "0", "1" = "1","2"="2")) + theme(strip.text.x = element_text(size = 30, colour = "black", face="bold"))
dev.off()

##Test overlap with the NSD1 lesion score (number of lesions)
ov=overlap(hypermean, NSD1_MutCGH)
tab_NSD1_MutCGH=ov[[1]]
NSD1_MutCGH_tab=ov[[2]]

#Test association between Sotos styndrome overlap index and the NSD1 lesion score using linear regression 
summary(glm(tab_NSD1_MutCGH$Percentage_Overlap~NSD1_MutCGH_tab$nlesions))
#p=1.09e-12

t1=cbind(tab_NSD1_MutCGH$Percentage_Overlap, NSD1_MutCGH_tab$nlesions)
colnames(t1)=c("Meth","Mut")
t1=as.data.frame(t1)
p=ggplot(t1, aes(factor(Mut), Meth))  

file="Directory/file"
png(file=paste(file,'.png',sep=''), units="in", width=11, height=8.5, res=600)
p + geom_boxplot( aes(fill =  factor(Mut)),outlier.shape=NA, position=position_dodge(0.8), width=0.8)+ scale_x_discrete("CNV", breaks=factor(1:4), drop=FALSE) + geom_jitter() + theme(legend.position='none', axis.text.y = black.bold.20.text, axis.title.x=element_blank(), legend.text=element_text(size=25, face="bold"),legend.title=element_text(size=25, face="bold"), axis.title.y=element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.justification=c(.5,.5)) + scale_fill_manual(name = "NSD1\nmutations (n)", values = cols, labels = c("-1"="-1","0" = "0", "1" = "1","2"="2")) 
dev.off()

#################
#CIBERSORT analysis
########################

#Get CIBERSORT data indicating predicted levels of each immune cell type
LoadFile='Dirctlry/file.tsv'
cs=read.delim(LoadFile, header = TRUE)

#split into tumor and normal tissue data 
csT=cs[substr(cs$SampleID,14,15)=="01",]
csT$patient_barcode=gsub("[.]","-",substr(csT$SampleID,1,12))

csN=cs[substr(cs$SampleID,14,15)=="11",]
csN$patient_barcode=gsub("[.]","-",substr(csN$SampleID,1,12))
#

#Overlap with subtypes
OverlapSamples=intersect(rownames(cluster5), csT$patient_barcode)
cluster5_csT=cluster5[OverlapSamples,]
csT_cluster5=csT[csT$patient_barcode %in% OverlapSamples,]

rownames(csT_cluster5)=csT_cluster5$patient_barcode

OverlapSamples=intersect(rownames(cluster5), csN$patient_barcode)
cluster5_csN=cluster5[OverlapSamples,]
csN_cluster5=csN[csN$patient_barcode %in% OverlapSamples,]

rownames(csN_cluster5)=csN_cluster5$patient_barcode

#########################
#Make heatmap of levels of immune variables in each subtype (as in Figure 3a)
##############################

#Get T cell signature (13 known T cell signatures)
TCellTranscripts=c("CD8A", "CCL2", "CCL3",
                   "CCL4", "CXCL9", "CXCL10", "ICOS", "GZMK", "IRF1", "HLA-DMA", "HLA-DMB", "HLA-DOA",
                   "HLA-DOB")
TCellTranscripts=checkGeneSymbols(TCellTranscripts)$Suggested.Symbol

#make matrix with the levels of cell types in each subtype
cellTypes=c("TotalLeukocyte","Plasma.cells", "Macrophages.M1","T.cells.CD8","T.cells.CD4.memory.resting")
tab=array(NA, c(length(cellTypes),5))

for(i in 1:5){
  for(j in 1:length(cellTypes)){
    tab[j,i]=mean(csT_cluster5[cluster5_csT$cluster==i,cellTypes[j] ])
  }
}      
rownames(tab)=cellTypes  
colnames(tab)=c("Non-CIMP-Atypical","NSD1","CIMP-Atypical","HPV","Stem-like_Smoking")
#tab=as.data.frame(tab)

#Add levels of expression of immunotherapy-related genes 
#tab=rowScales(tab, add.stats = FALSE)
genes=c("PDCD1","CD274","PDCD1LG2")

tab2=array(NA, c(length(genes)+1,5))

for(i in 1:5){
  for(j in 1:length(genes)){
    tab2[j,i]=mean(MAcancer_cluster5[genes[j], cluster5_MAcancer$cluster==i])
  }
}      
rownames(tab2)=c("PDCD1", "CD274", "PDCD1LG2", "TIL signature")  
colnames(tab2)=c("Non-CIMP-Atypical","NSD1","CIMP-Atypical","HPV","Stem-like_Smoking")

for(i in 1:5){
  tab2[4,i]=mean(colMeans(MAcancer_cluster5[TCellTranscripts,cluster5_MAcancer$cluster==i]))
}

tab2=tab2[ nrow(tab2):1, ]
tab2=t(apply(tab2,1,scale))
colnames(tab2)=c("Non-CIMP-Atypical","NSD1","CIMP-Atypical","HPV","Stem-like_Smoking")
tab2.m=melt(tab2)
colnames(tab2.m)=c("cells","subtype","mean")

#Make heatmap
p2=ggplot(tab2.m, aes(subtype, cells)) + geom_tile(aes(fill = mean), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
p2 + theme_grey(base_size = base_size) + labs(x = "",  y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())
ggsave(paste("Directory/file.png", sep=""), width = 10, height = 5)

#############################################
#Test association of immunotherapy-related genes with levels of different immune cell types across cancer types
###########################################

#PDL_data is a list of matrices, each with the expression data the immunotherapy-related genes in a different TCGA cancer type

#Make heatmap 
genes=c("PDCD1","CD274","PDCD1LG2")

for(n in 1:length(genes)){
  gene=genes[n]
  
  tab=array(NA, c(length(PDL_data), length(cells)))
  rownames(tab)=names(PDL_data)
  colnames(tab)=cells
  
  for(i in 1:length(PDL_data)){
    
    MA=PDL_data[[i]] 
    
    OverlapSamples=intersect(colnames(MA), csTu$patient_barcode)
    MA_csT=MA[,OverlapSamples]
    csT_MA=csTu[csTu$patient_barcode %in% OverlapSamples,]
    
    for(j in 1:length(cells)){
      cell=cells[j]
      tab[i,j]=ifelse(sd(csT_MA[, cell])!=0 & sd(MA_csT[gene,])!=0,cor(MA_csT[gene,], csT_MA[, cell], use = "complete.obs"),NA)
      
    }
  }
  
  tab=tab[!is.na(rowMeans(tab)),]
  tab=tab[ nrow(tab):1, ]
  #tab=t(apply(tab,1,scale))
  #colnames(tab)=c("cancer","cell","correlation")
  
  tab.m=melt(tab)
  colnames(tab.m)=c("cancer","cell","correlation")
  tab.m$cell=as.factor(tab.m$cell)
  levels(tab.m$cell)=cells
  
  p=ggplot(tab.m, aes(cancer, cell)) + geom_tile(aes(fill = correlation), colour = "white") + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") 
  
  p + coord_flip() + theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste("Directory/Heatmap_CIBERSORT_",gene,"_cor_with_TALs.png", sep=""), width = 3.5, height = 6)
}

########################
#All dependency scripts (Sourced in scripts above)
##########################

#Function to get 450k array data annotation for CpG probes

TCGA_GENERIC_GetInfinium450k_annotation <-function() {
  # start server using: /usr/local/mysql/bin/mysqld
  #conn<- dbConnect(MySQL(), user="root", password="basket", dbname="StanfordData")     
  #conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)
  
  pcname <- system('uname -n',intern=T) 
  if (regexpr('bmir',pcname)) conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
  if (pcname=='isis-genome.stanford.edu') conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3306)          
  ProbeAnnotation=dbGetQuery(conn, "SELECT ILMNID,GENESYMBOL FROM TCGA_MET_450k WHERE GENESYMBOL NOT IN (\"\");")     
  dbDisconnect(conn)
  return(ProbeAnnotation)
}

#Function to get 27k array data

TCGA_GENERIC_GetInfinium27k_annotation <-function() {
  #conn<- dbConnect(MySQL(), user="root", password="basket", dbname="StanfordData")
  #conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
  pcname <- system('uname -n',intern=T) 
  if (regexpr('bmir',pcname)) conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
  if (pcname=='isis-genome.stanford.edu') conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3306)               
  ProbeAnnotation=dbGetQuery(conn, "SELECT Illumina_ID,SYMBOL FROM TCGA_MET_27k WHERE SYMBOL NOT IN (\"\");")     
  dbDisconnect(conn)
  return(ProbeAnnotation)
}

#Fuction to get probes with SNPs

TCGA_GENERIC_GetInfinium450k_GetProbeswithSNPs <-function() {
  pcname <- system('uname -n',intern=T) 
  if (regexpr('bmir',pcname)) conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3307)     
  if (pcname=='isis-genome.stanford.edu') conn<- dbConnect(MySQL(), user="ogevaert", password="",host="127.0.0.1", dbname="StanfordData",port=3306)                    
  
  ProbeAnnotation=dbGetQuery(conn, "SELECT IlmnID,Probe_SNPs,Probe_SNPs_10 FROM TCGA_MET_450k_Extensive WHERE Probe_SNPs NOT IN (\"\") OR Probe_SNPs_10 NOT IN (\"\") ;")     
  #ProbeAnnotation=dbGetQuery(conn, "SELECT IlmnID,Probe_SNPs,Probe_SNPs_10,UCSC_RefGene_Name FROM TCGA_MET_450k_Extensive WHERE Probe_SNPs NOT IN (\"\") OR Probe_SNPs_10 NOT IN (\"\") ;")     
  
  dbDisconnect(conn)
  return(ProbeAnnotation[,1])
}


##Cluster the CpGs into clusters

TCGA_GENERIC_MET_ClusterProbes_with_hclust=function(MET_Cancer,MET_Normal,Parallel=FALSE) {
  # overlapping cancer & normal probes
  OverlapProbes=intersect(rownames(MET_Cancer),rownames(MET_Normal))
  MET_Cancer=MET_Cancer[OverlapProbes,]
  MET_Normal=MET_Normal[OverlapProbes,]
  
  #Get probe information
  #if (length(rownames(MET_Cancer))<30000) {
  #ProbeAnnotation_Old=TCGA_GENERIC_GetInfinium27k_annotation() # this contains the date genes errors !!!!          
  #MappingDataFile="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAovary/data/2011-04-15/METADATA/JHU_USC__HumanMethylation27/ProbeMapping.txt"
  #ProbeAnnotation <- read.table(MappingDataFile, header=FALSE, sep='\t')
  #ProbeAnnotation[,1]=as.character(ProbeAnnotation[,1])
  #ProbeAnnotation[,2]=as.character(ProbeAnnotation[,2])
  #} else {
  ProbeAnnotation=TCGA_GENERIC_GetInfinium450k_annotation()
  #}     
  
  # remove probes that are not annotated ??? (not necessary because I iterate over genes
  #MatchedProbes=match(rownames(MET_Cancer),ProbeAnnotation[,1])
  #MET_Cancer[is.na(MatchedProbes),]
  
  # remove probes with SNPs
  SNPprobes=TCGA_GENERIC_GetInfinium450k_GetProbeswithSNPs()     
  GoodProbes=setdiff(rownames(MET_Cancer),SNPprobes)
  NrProbesToRemove=length(rownames(MET_Cancer))-length(GoodProbes)
  cat("Removing",NrProbesToRemove,"probes with SNPs.\n")
  MET_Cancer=MET_Cancer[GoodProbes,]
  MET_Normal=MET_Normal[GoodProbes,]
  
  ###### only iterating over genes that have probes present
  # Getting the positions relative to probe annotation of the probes present in this data set. 
  PresentProbes=match(ProbeAnnotation[,1],rownames(MET_Cancer))
  UniqueGenes=sort(unique(ProbeAnnotation[!is.na(PresentProbes),2]))
  UniqueGenes=UniqueGenes[which(UniqueGenes != "")] # there is one empty one
  
  # create large matrix, delete zeros in the end.      
  MET_Cancer_C=matrix(0,length(rownames(MET_Cancer)),length(colnames(MET_Cancer)))
  MET_Normal_C=matrix(0,length(rownames(MET_Cancer)),length(colnames(MET_Normal)))
  colnames(MET_Cancer_C)=colnames(MET_Cancer)
  colnames(MET_Normal_C)=colnames(MET_Normal)
  ProbeMapping=matrix(0,length(rownames(MET_Cancer)),2)
  ProbeCounter=1
  METmatrixCounter=1
  ClusteredRownames=c()
  
  CorThreshold=0.4
  
  if (Parallel) {
    # alternative to do it in parallel
    NrCores=4
    cl <- makeCluster(NrCores)
    registerDoParallel(cl)          
    cat("Clustering",length(rownames(MET_Cancer)),"probes in CpG site clusters. Running in parallel mode with",NrCores,"cores.\n")
    # cluster needs to know the function because it is defined only here, so it does not know other functions.
    tmpClusterResults=foreach(i=1:length(UniqueGenes),.export='TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust') %dopar% {     #length(UniqueGenes)                              
      TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust(UniqueGenes[i],ProbeAnnotation,MET_Cancer,MET_Normal,CorThreshold)
    }
    #stopImplicitCluster(cl)   
    on.exit(stopCluster(cl))          
    for ( i in 1:length(UniqueGenes) ) {   #length(UniqueGenes)   
      # THIS IS VERY SLOW< rbind needs to reassign a matrix  every time. 
      #  MET_Cancer_C=rbind(MET_Cancer_C,tmpClusterResults[[i]][[1]])
      #  MET_Normal_C=rbind(MET_Normal_C,tmpClusterResults[[i]][[2]])              
      #  ProbeMapping=rbind(ProbeMapping,tmpClusterResults[[i]][[3]])                
      MET_Cancer_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[i]][[1]])),]=tmpClusterResults[[i]][[1]]               
      MET_Normal_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[i]][[2]])),]=tmpClusterResults[[i]][[2]]
      METmatrixCounter=METmatrixCounter+nrow(tmpClusterResults[[i]][[1]])
      ClusteredRownames=c(ClusteredRownames,rownames(tmpClusterResults[[i]][[1]]))
      ProbeMapping[ProbeCounter:(ProbeCounter-1+nrow(tmpClusterResults[[i]][[3]])),]=tmpClusterResults[[i]][[3]]
      ProbeCounter=ProbeCounter+nrow(tmpClusterResults[[i]][[3]])
    }
    # remove excessively large matrix
    MET_Cancer_C=MET_Cancer_C[1:length(ClusteredRownames),]
    MET_Normal_C=MET_Normal_C[1:length(ClusteredRownames),]
    rownames(MET_Cancer_C)=ClusteredRownames
    rownames(MET_Normal_C)=ClusteredRownames
    
  } else {   
    cat("Clustering",length(rownames(MET_Cancer)),"probes in CpG site clusters.\n")
    pb=txtProgressBar(1,length(UniqueGenes),style=3)
    for(i in 1:length(UniqueGenes)) {    #length(UniqueGenes)              
      setTxtProgressBar(pb,i)
      tmpClusterResults=TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust(UniqueGenes[i],ProbeAnnotation,MET_Cancer,MET_Normal,CorThreshold)
      # MET_Cancer_C=rbind(MET_Cancer_C,tmpClusterResults[[1]])
      # MET_Normal_C=rbind(MET_Normal_C,tmpClusterResults[[2]])
      # ProbeMapping=rbind(ProbeMapping,tmpClusterResults[[3]])
      
      MET_Cancer_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[1]])),]=tmpClusterResults[[1]]               
      MET_Normal_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[2]])),]=tmpClusterResults[[2]]
      METmatrixCounter=METmatrixCounter+nrow(tmpClusterResults[[1]])
      ClusteredRownames=c(ClusteredRownames,rownames(tmpClusterResults[[1]]))
      ProbeMapping[ProbeCounter:(ProbeCounter-1+nrow(tmpClusterResults[[3]])),]=tmpClusterResults[[3]]
      ProbeCounter=ProbeCounter+nrow(tmpClusterResults[[3]])
    }
    # remove excessively large matrix
    MET_Cancer_C=MET_Cancer_C[1:length(ClusteredRownames),]
    MET_Normal_C=MET_Normal_C[1:length(ClusteredRownames),]
    rownames(MET_Cancer_C)=ClusteredRownames
    rownames(MET_Normal_C)=ClusteredRownames
  }
  cat("\nFound",length(rownames(MET_Cancer_C)),"CpG site clusters.\n")
  return(list(MET_Cancer_Clustered=MET_Cancer_C,MET_Normal_Clustered=MET_Normal_C,ProbeMapping=ProbeMapping))
}

#Hierarchical clustering code

TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust <- function(Gene,ProbeAnnotation,MET_Cancer,MET_Normal=null,CorThreshold=0.3) {
  
  ####### first lookup the probes matching a single gene. DO NOT USE grep, it does not do exact matching, but looks for the pattern anywhere !!!
  Probes=ProbeAnnotation[which(ProbeAnnotation[,2] == Gene),1]
  Probes=Probes[which(Probes %in% rownames(MET_Cancer))]
  
  METcancer_Clustered=matrix(0,0,length(colnames(MET_Cancer)))
  if (!is.null(MET_Normal)) METnormal_Clustered=matrix(0,0,length(colnames(MET_Normal)))
  Clusternames=c()
  InverseCorrelationThreshold=1-CorThreshold
  GeneClustersForProbeMapping=array(dim=length(Probes))
  if (length(Probes)>1) {               
    ProbeCorrelation=cor(t(MET_Cancer[Probes,]),method='pearson')
    ClusterResults=hclust(as.dist(1-ProbeCorrelation), method = "complete", members = NULL)
    #plot(ClusterResults)
    Clusters=cutree(ClusterResults,h=InverseCorrelationThreshold)
    for (i in 1:length(unique(Clusters)) ) {
      tmpGeneProbes=Probes[Clusters==i]
      if ( length(tmpGeneProbes)>1 ) {
        
        tmpAveragedProfile=colMeans(MET_Cancer[tmpGeneProbes,])
        METcancer_Clustered=rbind(METcancer_Clustered,tmpAveragedProfile)          
        # Same for normal
        if (!is.null(MET_Normal)) tmpAveragedProfile=colMeans(MET_Normal[tmpGeneProbes,])
        if (!is.null(MET_Normal)) METnormal_Clustered=rbind(METnormal_Clustered,tmpAveragedProfile)
      } else {          
        METcancer_Clustered=rbind(METcancer_Clustered,MET_Cancer[tmpGeneProbes,])
        if (!is.null(MET_Normal)) METnormal_Clustered=rbind(METnormal_Clustered,MET_Normal[tmpGeneProbes,])
      }
      Clusternames=c(Clusternames,paste(Gene,'---Cluster',i,sep=""))               
      pos=which(Probes %in% tmpGeneProbes)
      GeneClustersForProbeMapping[pos]=paste(Gene,"---Cluster",i,sep="")
    }
    rownames(METcancer_Clustered)=Clusternames
    if (!is.null(MET_Normal)) rownames(METnormal_Clustered)=Clusternames         
    
  } else {
    METcancer_Clustered=MET_Cancer[Probes,,drop=FALSE]
    METnormal_Clustered=MET_Normal[Probes,,drop=FALSE]
    Clusternames=Gene
    GeneClustersForProbeMapping=Gene
    rownames(METcancer_Clustered)=Clusternames
    if (!is.null(MET_Normal)) rownames(METnormal_Clustered)=Clusternames
  }
  
  ProbeMapping=t(rbind(Probes,GeneClustersForProbeMapping))
  
  if (!is.null(MET_Normal)) {
    return(list(METcancer_Clustered,METnormal_Clustered,ProbeMapping))
  } else {
    return(list(METcancer_Clustered,ProbeMapping))
  }
  
}

#Script to fix TCGA patient sample barcodes

TCGA_GENERIC_CleanUpSampleNames <-function(GEN_Data,IDlength=12) {     
  SampleNames=colnames(GEN_Data)
  SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
  if (length(SampleNamesShort)!=length(unique(SampleNamesShort))) {
    # remove the doubles           
    Counts=table(SampleNamesShort)
    Doubles=rownames(Counts)[which(Counts>1)]
    
    cat("Removing doubles for",length(Doubles),"samples.\n")
    for(i in 1:length(Doubles)) {                         
      CurrentDouble=Doubles[i]          
      pos=grep(CurrentDouble,SampleNames)
      #GEN_Data[1:10,pos]
      #cor(GEN_Data[,pos])
      GEN_Data=GEN_Data[,-pos[2:length(pos)]]     
      SampleNames=colnames(GEN_Data) # need to update samplenames because pos is relative to this
    }
    SampleNames=colnames(GEN_Data)
    SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
    
    # now set the samplenames
    colnames(GEN_Data)=SampleNamesShort
  } else {
    colnames(GEN_Data)=SampleNamesShort     
  }     
  return(GEN_Data)
}

#Interect mutation data and DNA methylation cluster/subtype data

OverlapSamples=intersect(rownames(cluster5),colnames(MutSigMatrix))
cluster5_MutSigMatrix=cluster5[OverlapSamples,]
#289 samples
MutSigMatrix_cluster5=MutSigMatrix[,OverlapSamples]
#289 samples

#Run chi suqared test to identify signficantly mutated genes that are differentially mutated between DNA methylation clusters/subtypes
chisq_MutSigMatrix_SigGenes=array(NA, c(nrow(MutSigMatrix_cluster5), 1))
colnames(chisq_MutSigMatrix_SigGenes)=c("P.value")
rownames(chisq_MutSigMatrix_SigGenes)=rownames(MutSigMatrix_cluster5)

for(i in 1:nrow(MutSigMatrix_cluster5)){
  chisq_MutSigMatrix_SigGenes[i,1]=chisq.test(table(cluster5_MutSig$cluster,MutSigMatrix_cluster5[i,]), correct=T)$p.value
}
chisq_MutSigMatrix_SigGenes=as.data.frame(chisq_MutSigMatrix_SigGenes, drop=FALSE)
chisq_MutSigMatrix_SigGenes$Q.value=p.adjust(chisq_MutSigMatrix_SigGenes$P.value, method="fdr")
chisq_MutSigMatrix_SigGenes_rank<-chisq_MutSigMatrix_SigGenes[with(chisq_MutSigMatrix_SigGenes, order(chisq_MutSigMatrix_SigGenes$P.value)),]

#Functions to make boxplots of Sotos syndrome signature associated with DNA methylation subtypes
make_ggplot2_boxplot=function(tab1,tab2){
  ov=overlap(tab1, tab2)
  phyper_list_cluster5=ov[[1]]
  cluster5_phyper_list=ov[[2]]
  t1=cbind(phyper_list_cluster5$Percentage_Overlap, cluster5_phyper_list$cluster)
  colnames(t1)=c("Percent_overlap","cluster")
  t1=as.data.frame(t1)
  p1=ggplot(t1, aes(factor(cluster), Percent_overlap)) 
  t1=cbind(phyper_list_cluster5$Number_overlap, cluster5_phyper_list$cluster)
  colnames(t1)=c("Number_overlap","cluster")
  t1=as.data.frame(t1)
  p2=ggplot(t1, aes(factor(cluster), Number_overlap)) 
  t1=cbind(-log10(phyper_list_cluster5$P.value), cluster5_phyper_list$cluster)
  colnames(t1)=c("P.value","cluster")
  t1=as.data.frame(t1)
  p3=ggplot(t1, aes(factor(cluster), P.value))
  plots=list(p1,p2,p3)
  return(plots)
}

#Functions to make boxplots of Sotos syndrome signature associated with NSD1 mutations and deletions
makep_mutdel=function(tab1,tab2){
  ov=overlap(tab1, tab2)
  phyper_list_NSD1_MutCNV_LUSC=ov[[1]]
  NSD1_MutCNV_LUSC_phyper_list=ov[[2]]
  t1=cbind(phyper_list_NSD1_MutCNV_LUSC$Percentage_Overlap, NSD1_MutCNV_LUSC_phyper_list$Mut, NSD1_MutCNV_LUSC_phyper_list$CNV)
  colnames(t1)=c("Meth","Mut","CNV")
  t1=as.data.frame(t1)
  t1$CNV=as.factor(t1$CNV)
  t1$CNV=factor(t1$CNV,levels(t1$CNV)[c(4,3,2,1)])
  p=ggplot(t1, aes(factor(Mut), Meth))  
  return(p)
}

#Function to make matrix indicating the overlap between the hypomethylated CpG signatures of each patient and the Soto syndrome hypomethylated signature (or any other signature)
#It also get run the hypogeometric test for enrichment of the signature for each patient. 
make_hypergeometic=function(GeneSample, GeneSig, arraygenes){
  phyper_list=array(NA, c(length(GeneSample),3))
  rownames(phyper_list)=names(GeneSample)
  colnames(phyper_list)=c("P.value","Number_overlap","Percentage_Overlap")
  
  for(i in 1:length(GeneSample)){
    phyper_list[i,1]=phyper(length(intersect(GeneSample[[i]],GeneSig))-1,length(intersect(GeneSig, arraygenes)),length(setdiff(arraygenes,GeneSig)),length(GeneSample[[i]]),lower.tail=F, log.p = FALSE)
    phyper_list[i,2]=length(intersect(GeneSample[[i]],GeneSig))
    phyper_list[i,3]=length(intersect(GeneSample[[i]],GeneSig))/length(GeneSample[[i]])
  }
  phyper_list=as.data.frame(phyper_list)
  phyper_list$Q.value=p.adjust(phyper_list$P.value, method="fdr")
  return(phyper_list)
}

#Load GISTIC copy number data from data specified directory
TCGA_Load_GISTICdata <- function (GisticDirectory) {
  
  # loading the data in R
  GenesFile=paste(GisticDirectory,'all_data_by_genes.txt',sep="")
  CGH_Data=read.csv(GenesFile,sep="\t",row.names=1,head=TRUE,na.strings="NA")
  
  # Fix the sample names
  SampleNames=colnames(CGH_Data)
  SampleNames=gsub('\\.','-',SampleNames) # which way do we want to go . or -???
  colnames(CGH_Data)=SampleNames
  
  # TODO Process sample groups 
  # ADD CODE
  Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(CGH_Data))
  
  # CLeaning up sample names
  CGH_Data=TCGA_GENERIC_CleanUpSampleNames(CGH_Data,12)
  
  # remove the first two columns
  CGH_Data=CGH_Data[,-1]
  CGH_Data=CGH_Data[,-1]
  CGH_Data=as.matrix(CGH_Data)
  
  # loading the thresholded data. 
  GenesFile=paste(GisticDirectory,'all_thresholded.by_genes.txt',sep="")
  CGH_Data_Thresholded=read.csv(GenesFile,sep="\t",row.names=1,head=TRUE,na.strings="NA")
  
  # Fix the sample names
  SampleNames=colnames(CGH_Data_Thresholded)
  SampleNames=gsub('\\.','-',SampleNames) # which way do we want to go . or -???
  colnames(CGH_Data_Thresholded)=SampleNames
  
  CGH_Data_Thresholded=TCGA_GENERIC_CleanUpSampleNames(CGH_Data_Thresholded,12)
  CGH_Data_Thresholded=CGH_Data_Thresholded[,-1]
  CGH_Data_Thresholded=CGH_Data_Thresholded[,-1]
  CGH_Data_Thresholded=as.matrix(CGH_Data_Thresholded)
  
  # load the amp and del genes
  AMPfile=paste(GisticDirectory,'amp_genes.conf_99.txt',sep="")
  AMPtable=read.csv(AMPfile,sep="\t",head=TRUE,na.strings="NA")     
  AMPgenes=c()
  for (i in 2:ncol(AMPtable)) {
    tmpGenes=as.vector(AMPtable[4:nrow(AMPtable),i])
    tmpGenes[tmpGenes==""]=NA
    tmpGenes=na.omit(tmpGenes)
    AMPgenes=c(AMPgenes,tmpGenes)         
  }
  AMPgenes=unique(AMPgenes)
  AMPgenes=AMPgenes[order(AMPgenes)]
  
  # Same for DEL genes
  DELfile=paste(GisticDirectory,'del_genes.conf_99.txt',sep="")
  DELtable=read.csv(DELfile,sep="\t",head=TRUE,na.strings="NA")     
  DELgenes=c()
  for (i in 2:ncol(DELtable)) {
    tmpGenes=as.vector(DELtable[4:nrow(DELtable),i])
    tmpGenes[tmpGenes==""]=NA
    tmpGenes=na.omit(tmpGenes)
    DELgenes=c(DELgenes,tmpGenes)         
  }
  DELgenes=unique(DELgenes)
  DELgenes=DELgenes[order(DELgenes)]
  
  cat("There are",length(AMPgenes),"AMP genes and",length(DELgenes),"DEL genes.\n")
  
  return(list(CGH_Data_Segmented=CGH_Data,CGH_Data_Thresholded=CGH_Data_Thresholded,AMPgenes=AMPgenes,DELgenes=DELgenes))
}
}
