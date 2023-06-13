print("hello")
library(VirFinder)
help(VirFinder)
library('biomaRt')
library(tidyr)
library(dplyr)
install.packages("readxl")
library(readxl)
library(DESeq2)
install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
BiocManager::install("apeglm")


install.packages("ggrepel")



#Tez calismasi data processing

#Tez calismasi data eklenmesi
data_input1 = read.table("/Users/batuhanyolver/Downloads/Galaxy89-[featureCounts_on_data_36_and_data_88__Counts].tabular", header = TRUE, sep = "\t")
colnames(data_input1)[2] <- "ERR2902112"

data_input2 <- read.table("/Users/batuhanyolver/Downloads/Galaxy9-featureCounts_on_data_8_and_data_7__Counts.tabular", header = TRUE, sep = "\t")
colnames(data_input2)[2] <- "ERR2902102"

data_input3 <-  read.table("/Users/batuhanyolver/Downloads/Galaxy8-[featureCounts_ERR2902113_6__Counts].tabular", header = TRUE, sep = "\t")
colnames(data_input3)[2] <- "ERR2902113"

data_input4 <- read.table("/Users/batuhanyolver/Downloads/Galaxy108-[featureCounts_ERR2902108__Counts].tabular", header = TRUE, sep = "\t")
colnames(data_input4)[2] <- "ERR2902108"

data_input5 <- read.table("/Users/batuhanyolver/Downloads/Galaxy96-[featureCounts_on_data_36_and_data_95__Counts].tabular", header = TRUE, sep = "\t")
colnames(data_input5)[2] <- "SRR24804850_control"

data_input6 <- read.table("/Users/batuhanyolver/Downloads/SRR24804851_control.tabular", header = TRUE, sep = "\t")
colnames(data_input6)[2] <- "SRR24804851_control"

datname_tez = getBM(
  values =data_input1,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
  mart = useMart(biomart = "ensembl", dataset= "hsapiens_gene_ensembl")
)

names(datname_tez)[1] <- "Geneid"
join_of_tez <- full_join(datname_tez, data_input1, by="Geneid")
join_of_tez1 <- full_join(join_of_tez, data_input2, by="Geneid")
join_of_tez1 <- full_join(join_of_tez1, data_input3, by= "Geneid")
join_of_tez1 <- full_join(join_of_tez1, data_input4, by="Geneid")
join_of_tez1 <- full_join(join_of_tez1, data_input5, by="Geneid")
join_of_tez1 <- full_join(join_of_tez1, data_input6, by="Geneid")



statu_table <- read_excel("/Users/batuhanyolver/StatusOfData.xlsx", sheet = 1)
statu_table2 <- read_excel("/Users/batuhanyolver/StatusOfData.xlsx", sheet = 1)
statu_table <- statu_table[2]
rownames(statu_table) <- statu_table2$SampleID
statu_table$Status <- factor(statu_table$Status)
###
count_data <- join_of_tez1[,c(1,5,6,7,8,9,10)]
rownames(count_data) <- join_of_tez1$Geneid
#deseq analizi

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = statu_table, design = ~Status)

write.csv(count_data, "Count_verisi_7June.csv", row.names= FALSE)

a1 <- read.csv("Count_verisi_7June.csv", header= TRUE, row.names=1)

write.csv(statu_table, "STATU_Datasi.csv",row.names=FALSE)

a2 <- read.csv("STATU_Datasi.csv", header= TRUE, row.names=1)

a2$Status <- factor(a2$Status, levels = c("asd","control"))
a2$Status <- relevel(a2$Status, ref = "control")
dds1 <- DESeqDataSetFromMatrix(countData = a1, colData = a2, design = ~ Status)


dds1$Status <- factor(dds1$Status , levels = c("asd", "control"))

keep <- rowSums(counts(dds1)) >=5
dds2 <- dds1[keep,]

dds3 <- DESeq(dds2)

deseq_results <- results(dds3)

deseq_results <- as.data.frame(deseq_results)

deseq_results_ordered <- deseq_results[order(deseq_results$pvalue),]


filtered1 <-  deseq_results %>% filter(deseq_results$padj < 0.05)

filtered1 <- filtered1 %>% filter(abs(filtered1$log2FoldChange) > 1)


dim(deseq_results)

dim(filtered1)


plotDispEsts(dds3)


#PCA
vsd <- vst(dds3, blind=FALSE)

plotPCA(vsd, intgroup = c("Status"))

#generate the distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

#set color scheme
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(225)

#generate heatmap 
pheatmap(sampleDistMatrix, clustering_distance_rows= sampleDists,
         clustering_distance_cols = sampleDists, col = colors)

#top 10 hits, heatmap of log transformed normalized counts. We will use top 10 genes

top_hits <- deseq_results[order(deseq_results$padj), ][1:15,]
top_hits <- row.names(top_hits)

rld <- rlog(dds3, blind=FALSE)

pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE, show_rownames =T, cluster_cols = F)
pheatmap(assay(rld)[top_hits,],)


#MA Plot

plotMA(dds3, ylim=c(-2,2))

#remove the noise

resLFC <- lfcShrink(dds3, coef = "Status_control_vs_asd", type= "apeglm")

plotMA(resLFC, ylim=c(-2,2))

#volcano Plot

resLFC <- as.data.frame(resLFC)

resLFC$diffexpressed <- "NO"

resLFC$diffexpressed[resLFC$log2FoldChange > 0.1 & resLFC$padj <0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < 0.1 & resLFC$padj <0.05] <- "DOWN"

resLFC$delabel <- NA

ggplot(data= resLFC, aes(x=log2FoldChange, y= -log10(pvalue), col= diffexpressed, label= delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c("blue","black","red"))+
  theme(text= element_text(size=20))+
  geom_vline(xintercept = -10, color = "black", linetype = "solid", size = 1) +
  geom_vline(xintercept = 10, color = "black", linetype = "solid", size = 1) +
  geom_hline(yintercept = 150, color = "black", linetype = "solid", size = 1)
  
  



## Down ve up olmus gen isimleri
sonuc1 <- resLFC
sonuc1$ensemble_id <- rownames(sonuc1)

colnames(sonuc1)[8] <- "ensembl_gene_id"

datname_tez = getBM(
  values =sonuc1,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
  mart = useMart(biomart = "ensembl", dataset= "hsapiens_gene_ensembl")
)


sonuc2 <- sonuc1[sonuc1$log2FoldChange < -10 | sonuc1$log2FoldChange >10,]
sonuc3 <- sonuc2[sonuc2$pvalue > 1e-150,]

write.csv(sonuc3, "up-down-genler.csv", row.names= F)




