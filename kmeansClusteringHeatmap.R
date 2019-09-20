#set working directory
setwd("/Users/jleibowitz/")

final<-read.xlsx("annotategenes_TPM.xlsx") #load normalized gene sheet (can also pull from deseq2 object)
final<-final[match(genesForHeatmap,final$genes),] #keep only the identified sig genes
rownames(final)<-final$genes
final<-final[,-1]

forMatching<-c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11","A12",
               "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", 
               "C1","C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
               "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
               "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
               "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
               "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
               "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"
)
final<-final[,match(forMatching,colnames(final))] #order columns

colnames(final)<-paste0(colData[,3],colData[,5],colData[,1]) #rename colums to group


# calculate the z-score and generate heatmap 
library(gplots)
centered_data <- t(scale(t(data.matrix(final)), scale = T, center = T))


#determine cluster number
for(p in 1:ncol(final)){final[,p]<-as.numeric(final[,p])}
pca<-(prcomp(t(final), scale.=TRUE))
screeplot(pca) 

# perform k-means
set.seed(12345)
centers <- 6
final <- data.frame(centered_data)

# Assign genes to a cluster and sort clusters by 1-10
km <- kmeans(final, centers = centers)
m.kmeans<- cbind(final, cluster = km$cluster)
order <- m.kmeans[order(m.kmeans$cluster),]


# arrange clusters based on size (highest --> low )
len <- c()
for(i in 1:centers){
  len <- c(len, length(which(order$cluster == i)))
}

len2 <- order(-len)

sort.final <- data.frame()
for ( x in len2){
  sort.final <- rbind(sort.final, order[which(order$cluster == x),])
}

# assign colors to clusters
colours <-  colorRampPalette(brewer.pal(11,"Spectral"))(centers)[sort.final$cluster]

# generate heatmap with clusters 
sort.final<-data.matrix(sort.final)

#optional script for including interesting genes in row labels of heatmap
incl_gene_names<-vector(mode="character")
incl_gene_names<-c("^P2ry12$","^Apoe$","^Spp1$","^Gpnmb$","^Bhlhe40$","^Entpd1$","^Fcrls$","^Tmem119$","^Olfml3$","^Hexb$","^Irf8$","^Tgfbr1$",
                   "^Jun$","^Ccl2$","^C1qa$","^Bin1$","^Psen2$") #INPUT GENE NAMES TO INCLUDE ON HEAT MAP
rowlabels<-vector(mode="character", length=nrow(sort.final))

for (i in 1:length(incl_gene_names)){
  rowlabels[grep(incl_gene_names[i],rownames(sort.final))]<-incl_gene_names[i]
}

#order by cluster
sort.final<-sort.final[order(sort.final[,ncol(sort.final)]),]

library(RColorBrewer)

#create cluster annotation row (color annotation)
annotate_row<-as.data.frame(as.factor(sort.final[1:nrow(sort.final),ncol(sort.final)]))
rownames(annotate_row)<-rownames(sort.final)
colnames(annotate_row)<-"cluster"

sort.final2<-sort.final[,order(colnames(sort.final)[1:35])] #remove cluster column

#create column annotations (color annotation)
annotate_col<-as.data.frame(c(rep("Unstimulated",24),rep("TGFb",24),rep("IFNg",24),rep("TGFb_IFNg",24)))
annotate_col<-as.data.frame(c(rep("Female",48),rep("Male",48)))
annotate_col<-as.data.frame(c(rep("Healthy",10),rep("Mouse",6),rep("Progressive",9),rep("Relapsing",10)))
colnames(annotate_col)<-"Samples"
rownames(annotate_col)<-colnames(sort.final2)
Digestion<-as.factor(c(rep("No",5),rep("Yes",5),rep("No",2),rep("Yes",4),rep("No",4),rep("Yes",5),rep("No",5),rep("Yes",5)))
annotate_col$Digestion<-Digestion

breaksList<-seq(-2, 2, by = 0.25)

#generate heatmap
pheatmap(sort.final2,cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,annotation_row=annotate_row,labels_row = rownames(sort.final2),fontsize_row = 10,annotation_col = annotate_col, show_colnames = TRUE,show_rownames = FALSE)


# save the clusters and genes z-scores 
write.xlsx(sort.final2, 'clustered_genes.xlsx',
           row.names = TRUE)
