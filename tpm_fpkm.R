#script for outputting normalized gene values (tpm and fpkm)

setwd("/Users/jleibowitz/Desktop/SK-3TW7/transcripts_quant")
library(tximport)
library(openxlsx)
require(biomaRt)
ensMart<-useMart("ensembl")
ensembl_ms_mart<-useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl_df<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id","ensembl_transcript_id_version"), mart = ensembl_ms_mart)

tx2gene<-data.frame()
tx2gene[1:nrow(ensembl_df),1]<-ensembl_df$ensembl_transcript_id_version
tx2gene[1:nrow(ensembl_df),2]<-ensembl_df$ensembl_gene_id_version
colnames(tx2gene)<-c("Transcript","Gene")

dir<-dir()
sampleList<-dir[1:96] #can subset
files<-vector()
for (i in 1:length(sampleList)){
  files[i]<-paste0(sampleList[i],"/quant.sf")
}
names(files)<-sampleList

txi<-tximport(files=files,type="salmon",tx2gene = tx2gene)
######fpkm
perMillion<-vector()
fpkm<-data.frame()
for (i in 1:ncol(txi$counts)){
  perMillion[i]<-sum(txi$counts[,i])/1000000
  
  
  fpkm[1:nrow(txi$counts),i]<-(txi$counts[,i]/perMillion[i])/(txi$length[,i]/1000)
  
}


#for(t in 1:nrow(txi$counts)){
 # fpkm[t,i]<-(txi$counts[t]/perMillion[i])/txi$length[t,i]
#}

######

ensMart<-useMart("ensembl")
ensembl_ms_mart<-useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl_df<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end','transcript_length'), mart = ensembl_ms_mart)
my_genes=rownames(txi$abundance)
my_genes_ann=ensembl_df[match(my_genes,ensembl_df$ensembl_gene_id_version),]
final<-cbind(my_genes_ann$mgi_symbol,txi$abundance)

write.xlsx(final,file="annotatedgenes_TPM.xlsx")

fpkm<-cbind(final[,1],fpkm)
colnames(fpkm)<-colnames(final)
write.xlsx(fpkm,file="annotatedgenes_FPKM.xlsx")
