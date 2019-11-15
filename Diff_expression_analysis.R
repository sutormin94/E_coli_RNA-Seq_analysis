##############
#Script for the analysis of differential expression in RNA-Seq data obtained for E. coli DY330.
#Dmitry Sutormin, 2019.
#############

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
install.packages("xlsx")

library(edgeR)
require(xlsx)

###
#Identify dif exp genes among all the 3 conditions.
###
#Quick start

Ident_dif_exp <- function(exp_dataframe, group_factor, test_to_use, coef_to_use) {
  y <- DGEList(counts=exp_dataframe,group=group_factor)
  plotMDS(y, col=rep(1:3, each=3))
  
  message ("\nInitial samples description table:")
  print(y$samples)
  
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  message ("\nSamples table after low expressed genes are discarded:")
  print(y$samples)
  
  y <- calcNormFactors(y)
  message ("\nSamples table after samples are normalized on total reads number:")
  print(y$samples)
  plotMDS(y, col=rep(1:3, each=3))
  
  design <- model.matrix(~group_factor)
  y <- estimateDisp(y,design)
  
  if(test_to_use=="F-test"){
    #To perform quasi-likelihood F-tests:
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=coef_to_use)
    result <- qlf
    }
    
  
  if(test_to_use=="LRT"){
    #To perform likelihood ratio tests:
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=coef_to_use)
    result <- lrt
    }
  
  return(result)
}


#Read table into dataframe
x <- read.xlsx("F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis_correct/Genes_del_cor/DY330_RNA-Seq_data_merged_fragment_counts_columns_for_edgeR.xlsx", row.names="accession", sheetName="Sheet1")
group <- factor(c(1,1,1,2,2,2,3,3,3))
All_vs_all <- Ident_dif_exp(x, group, 'F-test', 2:3)
print(All_vs_all)
DE_All_vs_all_significance <- All_vs_all$table$PValue
hist(DE_All_vs_all_significance, breaks=seq(0,1, by=0.01), xlab="PValue", ylab="Number of genes", main="DE all vs all")

#Create subsets corresponding to pairs of culturing conditions
E_ES_x <- x[, c("Frag_count_EP1", "Frag_count_EP2", "Frag_count_EP3", "Frag_count_ESP1", "Frag_count_ESP2", "Frag_count_ESP3")]
E_S_x <- x[ , c("Frag_count_EP1", "Frag_count_EP2", "Frag_count_EP3", "Frag_count_SP1", "Frag_count_SP2", "Frag_count_SP3")]
ES_S_x <- x[ , c("Frag_count_ESP1", "Frag_count_ESP2", "Frag_count_ESP3", "Frag_count_SP1", "Frag_count_SP2", "Frag_count_SP3")]

group_pairs <- factor(c(1,1,1,2,2,2))
E_vs_ES <- Ident_dif_exp(E_ES_x, group_pairs, 'F-test', 2)
write.table(E_vs_ES$table, 'F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis_correct/Genes_del_cor/EdgeR_analysis_dif_expression/E_vs_ES_DE.txt', row.names=TRUE, col.names=TRUE, sep='\t')
DE_E_vs_ES_significance <- E_vs_ES$table$PValue
hist(DE_E_vs_ES_significance, breaks=seq(0,1, by=0.01), xlab="PValue", ylab="Number of genes", main="DE E vs ES")

E_vs_S <- Ident_dif_exp(E_S_x, group_pairs, 'F-test', 2)
write.table(E_vs_S$table, 'F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis_correct/Genes_del_cor/EdgeR_analysis_dif_expression/E_vs_S_DE.txt', row.names=TRUE, col.names=TRUE, sep='\t')
DE_E_vs_S_significance <- E_vs_S$table$PValue
hist(DE_E_vs_S_significance, breaks=seq(0,1, by=0.01), xlab="PValue", ylab="Number of genes", main="DE E vs S")

ES_vs_S <- Ident_dif_exp(ES_S_x, group_pairs, 'F-test', 2)
write.table(ES_vs_S$table, 'F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis_correct/Genes_del_cor/EdgeR_analysis_dif_expression/ES_vs_S_DE.txt', row.names=TRUE, col.names=TRUE, sep='\t')
DE_ES_vs_S_significance <- ES_vs_S$table$PValue
hist(DE_ES_vs_S_significance, breaks=seq(0,1, by=0.01), xlab="PValue", ylab="Number of genes", main="DE ES vs S")

print(topTags(E_vs_ES))
print(topTags(E_vs_S))
print(topTags(ES_vs_S))

###
###Remove rRNA genes before DE analysis?
###

par(mfrow=c(3,3))
hist(x$Frag_count_E1, main="E1", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_E2, main="E2", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_E3, main="E3", xlab="Expression level", ylab="Number of genes", breaks=10)

hist(x$Frag_count_ES1, main="ES1", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_ES2, main="ES2", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_ES3, main="ES3", xlab="Expression level", ylab="Number of genes", breaks=10)

hist(x$Frag_count_S1, main="S1", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_S2, main="S2", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(x$Frag_count_S3, main="S3", xlab="Expression level", ylab="Number of genes", breaks=10)
dev.off()


sort_return_top <- function(vec_in, thr) {
  vec_in_sorted <- sort(vec_in, decreasing = TRUE)
  vec_in__top <- vec_in_sorted[which(vec_in_sorted > thr)]
  return(vec_in__top)
}

par(mfrow=c(3,3))
hist(sort_return_top(x$Frag_count_E1, 100000), main="E1 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_E2, 100000), main="E2 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_E3, 100000), main="E3 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)

hist(sort_return_top(x$Frag_count_ES1, 100000), main="ES1 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_ES2, 100000), main="ES2 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_ES3, 100000), main="ES3 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)

hist(sort_return_top(x$Frag_count_S1, 100000), main="S1 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_S1, 100000), main="S2 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
hist(sort_return_top(x$Frag_count_S1, 100000), main="S3 top > 100000", xlab="Expression level", ylab="Number of genes", breaks=10)
dev.off()


###
###Fraction of rRNA.
###

exp_in_all_samples <- apply(x, 1, sum)
exp_in_all_samples_sorted <- sort(exp_in_all_samples, decreasing=TRUE)
write.table(exp_in_all_samples_sorted, 'F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis/Total_expression_of_genes.txt', row.names=TRUE, col.names=FALSE, sep='\t')

fraction_of_rRNA <- sum(exp_in_all_samples_sorted[which(exp_in_all_samples_sorted > 1800000)])/sum(exp_in_all_samples_sorted)
message('Mean fraction of rRNA in E. coli is:')
print(fraction_of_rRNA)

rRNA_names <- read.table('F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis/rRNA_genes_list.txt', sep="")
rRNA_names <- levels(unlist(rRNA_names))
rRNA_genes_expression <- x[rRNA_names, ]
rRNA_in_each_dataset <- apply(rRNA_genes_expression, 2, sum)
rRNA_in_each_dataset
RNA_in_each_dataset <- apply(x, 2, sum)
RNA_in_each_dataset

boxplot(as.vector(rRNA_in_each_dataset)/as.vector(RNA_in_each_dataset), ylab="Fraction of rRNA", main="Fraction of rRNA in RNA-Seq data")


###
###Remove the top expressed genes (rRNA genes).
###

not_rRNA_genes_expression <- x[which(!(row.names(x) %in% rRNA_names)), ]
All_vs_all_nr <- Ident_dif_exp(not_rRNA_genes_expression, group, 'F-test', 2:3)
print(All_vs_all_nr)
DE_All_vs_all_nr_significance <- All_vs_all_nr$table$PValue
hist(DE_All_vs_all_nr_significance, breaks=seq(0,1, by=0.01), xlab="PValue", ylab="Number of genes", main="DE all vs all no rRNA")


###
###Old expression data added.
###

x_ext <- read.delim("F:/E_coli_RNA-Seq/E_coli_DY330_RNA-Seq/FPKM_analysis/All_RNA-Seq_data_merged_RPKM_columns_plus_old_data_for_edgeR.txt",row.names="accession_trimmed")
group_ext <- factor(c(1,1,1,2,2,2,3,3,3,4))
All_vs_all_ext <- Ident_dif_exp(x_ext, group_ext, 'F-test', 2:4)
print(All_vs_all_ext)
