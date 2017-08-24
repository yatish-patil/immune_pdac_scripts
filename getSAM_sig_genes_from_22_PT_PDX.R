rm(list=ls())

# Load libraries
require(siggenes)
require(devtools)
require(roxygen2)
install_github("yatish-patil/iDissect",force = TRUE)  

### Creating folder structure

mainDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/"
dataDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/data/PAT_PDX/"
scrDir = "/Users/ypatil/Dropbox (ICR)/Yatish/Analysis/PDA/immune_pdac/new_again/scripts/immue_pdac_scripts/"

outDir <- paste0(mainDir,"/output_sam_significant_genes_2017_18_08/")
dir.create(outDir, showWarnings = FALSE)


### Reading 22 Patient and PDX file

infile = paste0(dataDir,"20170818_PDAC_22PT_22Xeno_immune_final_removedHKS_13_normalised_log2_ready.csv")
rawData = read.delim(infile,header = TRUE,row.names=1,sep = ","); dim(rawData)

### PCA analysis

pdf(paste0(outDir,Sys.Date(),'_pca_PDAC_22PT_22Xeno_immune_final_removedHKS_13_normalised_log2_ready.pdf'), width = 12, height = 9)
cl = rep(1:ncol(rawData))
plot(prcomp(t(apply(rawData,2,as.numeric)),scale=TRUE)$x[,1:2],col=cl,type = 'n')
text(prcomp(t(apply(rawData,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(rawData),col=cl)
dev.off()

# Creating a data for xeno samples
rawpdxPos = grep("XENO",colnames(rawData))
rpdxData = rawData[,rawpdxPos]; dim(rpdxData)
# Creating a data for patient samples
rawpatPos = grep("XENO",colnames(rawData), invert = TRUE)
rpatData = rawData[,rawpatPos]; dim(rpatData)
# combine the two datasets
data=cbind(rpdxData,rpatData)
group=rep(c(0,1),c(22,22))

### Writing first pdx 22 samples and second patient 22 samples to file

dataFile = paste0(outDir,Sys.Date(),"_order_samples_20170818_PDAC_22PT_22Xeno_immune_final_removedHKS_13_normalised_log2_ready.txt",sep="")
data_genes=cbind(rownames(data),data)
colnames(data_genes) <- c("Genes",colnames(data))
write.table(data_genes,dataFile,sep = "\t",quote = FALSE,row.names=FALSE)


## sam analysis
sam_pat_pdx = sam(data,factor(group),rand=123,gene.names = rownames(data))

summary_sig_de_genes=summary(sam_pat_pdx,6.4)@row.sig.genes
sam_sig_genes_data = data[summary_sig_de_genes, ]; dim(sam_sig_genes_data)
samDataFile = paste0(outDir,Sys.Date(),"_289_sam_signficant_genes_PDAC_22PT_22Xeno_immune_final_removedHKS_13_normalised_log2_ready.txt",sep="")
sam_data_genes=cbind(rownames(sam_sig_genes_data),sam_sig_genes_data)
colnames(sam_data_genes) <- c("Genes",colnames(sam_sig_genes_data))
write.table(sam_data_genes,samDataFile,sep = "\t",quote = FALSE,row.names=FALSE)

### Run inSillico microdissection on 289 sam significant genes
### iDissect will give us immune specific genes where it is highly expressed in patient and lowly expressed in pdx
### with low variability in pdx condition
setwd(outDir)
samData = read.delim(samDataFile,header = TRUE,sep="\t",row.names = 1);dim(samData)
iDissect::iDissect(vData = samData,method = "sd",cond1 = 22,cond2 = 22,fdrCutOff = 80)
iDissect::iDissect(vData = samData,method = "sd",cond1 = 22,cond2 = 22,fdrCutOff = 85)
iDissect::iDissect(vData = samData,method = "sd",cond1 = 22,cond2 = 22,fdrCutOff = 90)
iDissect::iDissect(vData = samData,method = "sd",cond1 = 22,cond2 = 22,fdrCutOff = 95)
