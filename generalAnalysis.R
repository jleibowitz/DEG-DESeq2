setwd("/users/jleibowitz/Lien tests")
sampleTable <- data.frame(type = c("neg","pos","neg","pos","neg","neg","neg","neg","neg","pos","neg","pos","pos","pos","pos","pos"),condition=c("unfl","fl","unfl","fl","fl","unfl","fl","unfl","fl","unfl","fl","unfl","fl","unfl","fl","unfl"))

files<-list.files(pattern="*.results")
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)



for (i in 1:nrow(txi.rsem$length)){for (t in 1:ncol(txi.rsem$length)){if (txi.rsem$length[i,t]==0){txi.rsem$length[i,t]=1}}}
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~type*condition)
dds$group<-factor(paste0(dds$type,dds$condition))
design(dds)<-~group
dds<-(DESeq(dds))
resultsNames(dds)

restest<-results(dds,contrast=c("group","posunfl","posfl"))