setwd("/Users/jleibowitz/") #Set to directory containing metadata

library('openxlsx')
library('tximport')
library('DESeq2')

require(biomaRt)

ensMart<-useMart("ensembl") #connecting to ensembl to pull transcriptome information
ensembl_ms_mart<-useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl_df<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id","ensembl_transcript_id_version", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end','transcript_length'), mart = ensembl_ms_mart)

tx2gene<-data.frame() #create a transcript-to-gene data.frame for summarizing transcript level information
tx2gene[1:nrow(ensembl_df),1]<-ensembl_df$ensembl_transcript_id_version
tx2gene[1:nrow(ensembl_df),2]<-ensembl_df$ensembl_gene_id_version
colnames(tx2gene)<-c("Transcript","Gene")

sampleTable<-read.xlsx("sampleTable.xlsx") #reading metadata
colData<-sampleTable
colData[,3]<-as.factor(colData[,3])
colData[,4]<-as.factor(colData[,5])


setwd("/Users/jleibowitz/") #Set to directory containing quant.sf ouput from salmon (or other alignment softwares)

filenames<-c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11","A12",
             "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", 
             "C1","C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
             "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
             "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
             "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
             "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
             "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"
)

filenames<-paste0(filenames,"/quant.sf")

txi <- tximport(files=filenames, type = "salmon",tx2gene = tx2gene) #import aligned data and summarize to gene level

dds<-DESeqDataSetFromTximport(txi, colData, ~Disease) #Deseq2 - model must be colname from colData (metadata)
nrow(dds)
dds <- estimateSizeFactors(dds)
dds<-dds[rowSums(counts(dds, normalized = TRUE) > 10) >= 5] #filter out lowly express genes (optional)
nrow(dds)

vsd <- vst(dds, blind = FALSE) #normalization


library("dplyr")
library("ggplot2")

library("pheatmap")
library("RColorBrewer")

plotPCA(vsd, intgroup = c("Disease")) #visualizing differences between samples on PCA plots
plotPCA(vsd, intgroup = c("Digestion"))



dds$group <-factor(paste0(dds$Disease,dds$Digestion), level="Group1","Group2") #creating new groups combining all factors into one for LRT test, and leveling for base group
#dds$group <- factor(paste0(dds$Stimulation)) #stimulation effect
#dds$group <-factor(paste0(dds$Gender)) #gender effect
#dds$group <-factor(paste0(dds$Genotype)) #genotype effect
design(dds) <- ~ group #setting new group
dds <- DESeq(dds)


dds_lrt <- DESeq(dds,test="LRT",reduced = ~1) #LRT
res_LRT <- results(dds_lrt,alpha=0.05) #results

res_LRT_df<-as.data.frame(res_LRT)
res_LRT_df<-res_LRT_df[order(res_LRT_df$padj),] #ordering by padj (FDR)
orderedRes<-res_LRT_df


genes<-ensembl_df[match(rownames(orderedRes),ensembl_df$ensembl_gene_id_version),3] #matching ensembl ids and pulling gene symbols

final<-cbind(genes,orderedRes) #statistical data.frame

genesForHeatmap<-unique(final$genes[which(final$padj<0.05)]) #pulling genes for heatmap
