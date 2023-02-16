library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)

# Reading RNA-seq from MCF7 cell lines
Genes_TAMOH <- read.csv("../data/Unnormalized_Cell_Lines_Gene_Counts_TAM.csv")

# Fixing gene names and removing low counts reads
Genes_TAMOH <- Genes_TAMOH[-which(is.na(Genes_TAMOH$Name)),]
rownames(Genes_TAMOH) <- Genes_TAMOH[,2]
Genes_TAMOH_nz <- Genes_TAMOH[rowSums(Genes_TAMOH[,3:8]) >= 20 ,] 

################ DESeq2  ######################
condition <- data.frame(Cell_line= c("MCF7_OH_1", "MCF7_OH_2", "MCF7_OH_3",
                                     "MCF7_TAMR_1", "MCF7_TAMR_2", "MCF7_TAMR_3"),
                        Condition = as.factor(c("OH","OH","OH", "TAMR",  "TAMR",  "TAMR")))

deseq_TAMOH <- SummarizedExperiment(assays =  as.matrix(Genes_TAMOH_nz[,3:8]), #
                                    rowData = Genes_TAMOH_nz[,2],
                                    colData= condition) 

dds_TAMOH <- DESeqDataSet(deseq_TAMOH, design = ~ Condition)
dds_TAMOH <- DESeq(dds_TAMOH)
norm_counts_TAMOH <- as.data.frame(t(counts(dds_TAMOH, normalized = T)))
res_TAMOH <- results(dds_TAMOH,  contrast = c("Condition","TAMR", "OH")) 
res_TAMOH$SYMBOL<- rownames(res_TAMOH)
res_TAMOH <- as.data.frame(res_TAMOH)

## TCGA tamoxifen treated patients - Effective treatment vs resistant

#Reading response file
Response_TAM <- read.csv("../data/Response_TAM.csv")
cols <- c("Effective_treatment" = "green", "Resistant" = "red")
Response_TAM$Patient <-  gsub(x=Response_TAM$Patient, "\\.", "\\-")

#Change label for response
Response_TAM$Response <- gsub(x=Response_TAM$Response, "Good", "Effective_treatment")
Response_TAM$Response <- gsub(x=Response_TAM$Response, "Bad", "Resistant")

#Reading RNA-seq from TCGA TAM 
Genes_TAM <- read.csv("../data/Unnormalized_TCGA_Gene_Counts_TAM.csv")
colnames(Genes_TAM) <-  gsub(x=colnames(Genes_TAM), "\\.", "\\-")

# Fixing gene names and removing low counts reads
Genes_TAM <- Genes_TAM[,c(1,which(colnames(Genes_TAM) %in% Response_TAM$Patient))]
rownames(Genes_TAM) <- Genes_TAM[,1]
Genes_TAM <- Genes_TAM[,-1]
Genes_TAM_nz <- Genes_TAM[rowSums(Genes_TAM)>20,]

# Keep only patients pre-selected in the classification phase
Ordered_Patients <- data.frame(Patient=colnames(Genes_TAM_nz))

condition_TAM <- inner_join(Ordered_Patients, Response_TAM[,c(1,10)])
condition_TAM$Response <- as.factor(condition_TAM$Response)

################ DESeq2  ######################
deseq_TAM <- SummarizedExperiment(assays =  as.matrix(Genes_TAM_nz[,1:length(Genes_TAM_nz)]),
                                  rowData = rownames(Genes_TAM_nz),
                                  colData = condition_TAM)

dds_TAM <- DESeqDataSet(deseq_TAM, design = ~ Response)
dds_TAM <- DESeq(dds_TAM)
norm_counts_TAM <- as.data.frame(t(counts(dds_TAM, normalized = T)))
res_TAM <- results(dds_TAM,  contrast = c("Response","Resistant", "Effective_treatment"))
res_TAM$SYMBOL<- rownames(res_TAM)
res_TAM <- as.data.frame(res_TAM)


#### Results

res_join <- inner_join(res_TAMOH, res_TAM, by="SYMBOL")
res_join_plot <- res_join[,c(7,2,5,6,9,12,13)]
colnames(res_join_plot) <- c("SYMBOL", "log2FC_MCF7", "p_value_MCF7", "p_adj_MCF7", "log2FC_TCGA", "p_value_TCGA", "p_adj_TCGA")

## Apply filters to find common genes expressed enough in the same direction
res_join_plot$Significant <- ifelse(res_join_plot$p_adj_MCF7 < 0.1 & res_join_plot$p_adj_TCGA < 0.1,
                                    "Sig. Both", "Sig. Just One")

res_join_plot$DE <- ifelse(abs(res_join_plot$log2FC_MCF7) < 0.5 | abs(res_join_plot$log2FC_TCGA) < 0.5,
                           "Not enough DE", "DE")
res_join_plot$Sign <- ifelse((res_join_plot$log2FC_MCF7 * res_join_plot$log2FC_TCGA > 0),
                             "Same", "Opposite")
res_join_plot$DEandSign <- ifelse(res_join_plot$DE=="DE" & res_join_plot$Sign=="Same",
                                  "DE OK", "DE Fail")

res_join_plot$Finalcut <- ifelse(res_join_plot$Significant=="Sig. Both" & res_join_plot$DEandSign=="DE OK",
                                 "Signature Candidate", "Remove")

res_join_plot$SYMBOL <- ifelse(res_join_plot$Significant=="Sig. Both" & res_join_plot$Finalcut=="Signature Candidate",
                               res_join_plot$SYMBOL, "")
# Extract final canditate genes
Final_gene_selection <- res_join_plot$SYMBOL[which(res_join_plot$Finalcut=="Signature Candidate")]
Final_gene_selection


############ Plotting ######

### Heatmap plots (use TAM instead of TAMOH for TCGA plot)
vsd = vst(dds_TAMOH, blind=FALSE)
vst_counts = as.matrix(assay(vsd))

candidate_genes_sig <- res_TAMOH %>% 
  filter(padj < 0.1, abs(log2FoldChange) > 0.5) %>%    # filter table
  pull(SYMBOL) %>%          # extract the gene column as a vector
  unique()                   # retain only unique values


hclust_matrix_sig <-norm_counts_TAMOH %>% 
  as.matrix() 
hclust_matrix_sig <- t(hclust_matrix_sig[,candidate_genes_sig ])


#Renaming for TAM cells (comment and uncomment next line for TCGA)
colnames(hclust_matrix_sig) <- c("MCF7c_1", "MCF7c_2", "MCF7c_3", "TamR_1","TamR_2", "TamR_3")
#colnames(hclust_matrix_sig) <- gsub("TCGA-", "", colnames(hclust_matrix_sig))

hclust_matrix_sig <- hclust_matrix_sig %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

gene_dist_sig <- dist(hclust_matrix_sig)
gene_hclust_sig <- hclust(gene_dist_sig, method = "cen")
my.legend <-as.expression(bquote('log'['2']*' FC '))

#Plot heatmap 
Heatmap(hclust_matrix_sig, show_row_names = F, show_row_dend = F,
        column_names_rot=90,
        heatmap_legend_param = list(
          title = "Z-score", at = c(-4, 0, 5)))


### Gene Cloud plots in Supp. Material (add different sections +" for progressive plots)
# Plotting options for gene cloud plots
annotations <- data.frame(
  xpos = c(-8,-8,8,8),
  ypos =  c(-5, 5,-5,5),
  annotateText = c("Same Sign","Opposite Sign"
                   ,"Opposite Sign", "Same Sign"),
  hjustvar = c(0,0,1,1) ,
  vjustvar = c(0,1,0,1)) #<- adjust

ggplot(res_join_plot, aes(x=log2FC_MCF7, y= log2FC_TCGA, label=SYMBOL)) +
  #geom_point(colour="darkgreen", size =2)  +
  
  # geom_point(aes(colour=DE, alpha=DE)) +
  # scale_color_manual(values = c("High DE" = "darkgreen", "Low DE" = "grey50"), labels=NULL) +
  # scale_alpha_manual(values =c("Low DE"=0.2, "High DE" = 1) , labels=NULL) +
  # scale_size_manual(values=  c("Low DE" = 2, "High DE" = 3)) +
  # 
  # geom_point(aes(colour=DEandSign, alpha=DEandSign)) +
  #  scale_color_manual(values = c("In" = "darkgreen", "Out" = "grey50"), labels=NULL) +
  #  scale_alpha_manual(values =c("In"= 1, "Out" = 0.2 ), labels=NULL) +
  #  scale_size_manual(values=  c("Out" = 2, "In" = 3)) +
#  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText, fontface = "bold"), size=6) +

geom_point(aes(colour=Finalcut, alpha=Finalcut, size=Finalcut)) +
  scale_color_manual(values = c("Signature Candidate" = "darkgreen", "Remove" = "grey50"), labels=NULL) +
  scale_alpha_manual(values =c("Remove"=0.2, "Signature Candidate" = 1 , labels=NULL)) +
  scale_size_manual(values=  c("Signature Candidate" = 3, "Remove" = 2)) +
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText, fontface = "bold"), size=6) +
  
  
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-8,8)) +
  scale_y_continuous(limits = c(-5,5)) +
  geom_hline(yintercept = 0.5,  linetype="dotted") + 
  geom_vline(xintercept = 0.5, linetype="dotted") +
  geom_hline(yintercept = -0.5,  linetype="dotted") + 
  geom_vline(xintercept = -0.5, linetype="dotted") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     axis.title = element_text(face="bold", size = 20),
                     axis.text = element_text(face="bold", size = 20), 
                     legend.position = "none") + 
  geom_label_repel(aes(label=SYMBOL),  max.overlaps = Inf, force=12.6, fontface = "bold") 



  
  
  
  
  geom_point(aes(colour=Significant)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) #Add clour legend for significance

