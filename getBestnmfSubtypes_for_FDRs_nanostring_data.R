rm(list=ls()) 

require(ConsensusClusterPlus)
require(NMF)
require(sva)
require(devtools)
#install_github("syspremed/exploBATCH")
require(exploBATCH)
############## Creating Folder structure

mainDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/"
dataDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/data/Patient/"
scrDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/scripts/immune_pdac_scripts/"
iDissectOut = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/output_sam_significant_genes_2017_18_08/microDissectOutput/"

outDir <- paste0(mainDir,"/output_removed_4_samples_nanostirng_patient_95_NMF_2017_24_08/")
dir.create(outDir, showWarnings = FALSE)

nmfDir <- paste0(outDir,"/output_95_NMF_2017_24_08/")
dir.create(nmfDir, showWarnings = FALSE)

#Reading File of immune specific genes from FDR cutoff - 80,85,90 and 95

imm_specFile = paste0(iDissectOut,"2017-08-18_significant_specific_genes_fVM_95_sd.txt")
cutoff = strsplit(imm_specFile,"_")[[1]][14]
#cutoff = as.numeric(cutoff)*100

imm_spec_data = read.delim(imm_specFile,header = TRUE,row.names = 1);dim(imm_spec_data)
imm_genes=dim(imm_spec_data)[1]
cutoff_genes = paste0(cutoff,"_",imm_genes)


# Reading nanostring patient file  -64 samples
inputPATFile1 = paste0(dataDir,"20170818_PT-72(8 removed)_final_RemovedHKS_4_normalised_log2_ready.csv")
nanostring_DataPAT1 = read.delim(inputPATFile1,header = TRUE,sep=",",row.names = 1);dim(nanostring_DataPAT1)
PAT1_samples = dim(nanostring_DataPAT1)[2]
# PCA plots of patient file - 64 sample
pdf(paste0(outDir,Sys.Date(),'_nanostring_PCA_PAT_',PAT1_samples,"_",cutoff_genes,'_log2_normalization.pdf'), width = 12, height = 9)
plot(prcomp(t(apply(nanostring_DataPAT1,2,as.numeric)),scale=TRUE)$x[,1:2],type = 'n')
text(prcomp(t(apply(nanostring_DataPAT1,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(nanostring_DataPAT1))
dev.off()

# Reading nanostring patient file - 32 samples
inputPATFile2 = paste0(dataDir,"20170802_PDAC_immune_NewRLF_final_removed_3290T_32PT_RemovedHKS_3_normalised_log2_ready.csv")
nanostring_DataPAT2 = read.delim(inputPATFile2,header = TRUE,sep=",",row.names = 1);dim(nanostring_DataPAT2)
nanostring_DataPAT2[, c('X985T','X1061T','X1128T','X1433T')] <- list(NULL);dim(nanostring_DataPAT2)
PAT2_samples = dim(nanostring_DataPAT2)[2]

# PCA plots of patient file - 32-4 = 28 sample
pdf(paste0(outDir,Sys.Date(),'_nanostring_PCA_PAT_',PAT2_samples,"_",cutoff_genes,'_log2_normalization.pdf'), width = 12, height = 9)
plot(prcomp(t(apply(nanostring_DataPAT2,2,as.numeric)),scale=TRUE)$x[,1:2],type = 'n')
text(prcomp(t(apply(nanostring_DataPAT2,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(nanostring_DataPAT2))
dev.off()

## combine the dataset and do combat/exploBATCH anlysis
m<-match(toupper(rownames(nanostring_DataPAT1)),toupper(rownames(nanostring_DataPAT2)))
w<-which(!is.na(m))
nanostring_DataPAT1m<-nanostring_DataPAT1[w,]; dim(nanostring_DataPAT1m)
nanostring_DataPAT2m<-nanostring_DataPAT2[m[w],]; dim(nanostring_DataPAT2m)
nanostringPATData <-cbind(nanostring_DataPAT1m,nanostring_DataPAT2m);dim(nanostringPATData)
sampleNo = dim(nanostringPATData)[2]
combined_data_nanostring = paste0(outDir,Sys.Date(),"_combined_data_nanostring_",sampleNo,".txt")
write.table(nanostringPATData,combined_data_nanostring,sep="\t",quote=FALSE,row.names = FALSE)

# matching immune specific genes from different FDRs- 80,85,90 and 95
m<-match(toupper(rownames(imm_spec_data)),toupper(rownames(nanostringPATData)))
w<-which(!is.na(m))
imm_from_nanostringPATData<-nanostringPATData[m[w],]; dim(imm_from_nanostringPATData)


# Running combat
bat<-rep(c(1,2),c(dim(nanostring_DataPAT1m)[2],dim(nanostring_DataPAT2m)[2]))
length(bat)
combat<-ComBat(dat=nanostringPATData,batch=bat,mod=NULL, par.prior=TRUE, prior.plots=FALSE)

pdf(paste0(outDir,Sys.Date(),'_pca_after_combat_nanostringPAT_all_',sampleNo,'.pdf'), width = 12, height = 9)
cl = rep(1:ncol(combat))
plot(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],col=cl,type = 'n')
text(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(combat),col=cl)
dev.off()

#pdf(paste0(outDir,Sys.Date(),'_pca_after_combat_nanostringPAT_4_samples',sampleNo,'.pdf'), width = 12, height = 9)
#rm = c('X985T','X1061T','X1128T','X1433T')
#lc=colnames(combat)%in%rm
#ind=rep(.001,ncol(combat))
#ind[lc]=2
#cl = rep(1:ncol(combat))
#plot(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],col=cl,type = 'n')
#text(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(combat),col=cl,cex = ind)
#dev.off()

# Running exploBATCH
expBdat = t(as.data.frame(nanostringPATData))
dim(expBdat)
setwd(outDir)
expBATCH(D = expBdat,batchCL = bat,mindim = 2,maxdim = 9,method = "ppcca",SDselect = 0)

# Reading exploBATCH corrected after ppcca method file
expBFile = paste0(outDir,"/exploBATCHresults/correctBATCH/2017-08-24_ppccaCorrectedData.txt")
ppccaCorrectData = read.table(expBFile,header=TRUE,sep="\t",row.names = 1)
dim(ppccaCorrectData)
tr_ppccaCorrectData = t(ppccaCorrectData)
dim(tr_ppccaCorrectData)

tr_ppccaCorrectDataFile = paste0(outDir,Sys.Date(),"2017-08-24_transpose_ppccaCorrectedData_SD_0.txt",sep="")
tr_ppccaCorrectData_genes=cbind(rownames(tr_ppccaCorrectData),tr_ppccaCorrectData)
colnames(tr_ppccaCorrectData_genes) <- c("Genes",colnames(tr_ppccaCorrectData))
write.table(tr_ppccaCorrectData_genes,tr_ppccaCorrectDataFile,sep="\t",quote=FALSE,row.names = FALSE)

# matching immune specific genes from different FDRs- 80,85,90 and 95
m<-match(toupper(rownames(imm_spec_data)),toupper(rownames(tr_ppccaCorrectData)))
w<-which(!is.na(m))
imm_spec_from_expBatch_nanostringPATData<-tr_ppccaCorrectData[m[w],]; dim(imm_spec_from_expBatch_nanostringPATData)

# median centering the explobatch nanostring data
row_med_cent_imm_spec_from_PAT = apply(imm_spec_from_expBatch_nanostringPATData,1,median)
med_cent_imm_spec_from_PAT = sweep(imm_spec_from_expBatch_nanostringPATData,1,row_med_cent_imm_spec_from_PAT,"-");dim(med_cent_imm_spec_from_PAT)

med_cent_imm_spec_from_PAT_file = paste0(outDir,Sys.Date(),"_after_explobatch_median_centering_with_PAT_samples_",cutoff_genes,"_imm_sp_genes.txt")
write.table(med_cent_imm_spec_from_PAT,med_cent_imm_spec_from_PAT_file,sep="\t",quote=FALSE)

############# NMF Analysis ###### 
setwd(nmfDir)
source(paste0(scrDir,"/nmfconsensus.R",sep = ""))
cat("Performing NMF analysis...")
nmfconsensus(input.ds=med_cent_imm_spec_from_PAT_file,k.init=2,k.final=10,num.clusterings=20,maxniter=500, error.function="euclidean")

nmfFile = paste0(nmfDir,"consensus.k.4.gct")
subtypes = read.delim(nmfFile,header = TRUE)

order_subtypes=subtypes[order(subtypes$membership.ordered),]

colnames(order_subtypes) = c("samples","cophenetic_subtypes")
subtype_out = paste0(outDir,Sys.Date(),"_4_order_subtypes_file_","with_nanostring_PAT_samples_",cutoff_genes,".txt",sep="")
write.table(order_subtypes,subtype_out,sep="\t",row.names = FALSE,quote=FALSE)


t_order_samples = t(imm_spec_from_expBatch_nanostringPATData)[ order(match(rownames(t(imm_spec_from_expBatch_nanostringPATData)), order_subtypes$Name)), ]
subtypes_order_samples=t(t_order_samples);dim(subtypes_order_samples)

order_sampleFile = paste0(outDir,Sys.Date(),"_4_subtypes_data_order_","with_nanostring_PAT_samples_",cutoff_genes,"_",sampleNo,"_imm_sp_genes.txt",sep="")
subtypes_order_samples_genes=cbind(rownames(subtypes_order_samples),subtypes_order_samples)
colnames(subtypes_order_samples_genes) <- c("Genes",colnames(subtypes_order_samples))
write.table(subtypes_order_samples_genes,order_sampleFile,sep="\t",quote=FALSE,row.names = FALSE)



#####

cat("Performing NMF.......")
setwd(outDir)
mdat=med_cent_imm_spec_from_PAT
datapos <- nneg(apply(mdat,2,as.numeric), method = 'pmax')
rownames(datapos) <- rownames(med_cent_imm_spec_from_PAT)

set.seed(1)
log2_nmfccp <- NMF::nmf(datapos, 2:10, method = 'brunet', nrun = 100)
pdf(paste0(Sys.Date(),"_with_nanostring_PAT_samples_",cutoff_genes,"_imm_sp_genes",'_human_log2_nmf_plots.pdf'), width = 9, height = 6)
print(plot(log2_nmfccp))
consensusmap(log2_nmfccp)
dev.off()
log2_nmfccp$measures$cophenetic

