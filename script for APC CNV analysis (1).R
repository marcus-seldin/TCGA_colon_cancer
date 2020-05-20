library(dplyr)
library(tidyr)
library(data.table)
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(pheatmap)
library(colormap)
library(WGCNA)
library(ggrepel)
library(ggplot2)
library(reshape2)
setwd('D:/My Drive/lab files/lab files/selma masri/APC snp focused analysis/')

#Import Gene Annotations
gene_annots = read.delim('HGNC_annotations.txt')

#Read in the first type of file - CNV data
dt = fread("TCGA-COAD.gistic.tsv.gz")
dt = as.data.frame(dt)
dt$gene = factor(dt$`Gene Symbol`)
row.names(dt) = dt$gene
dt$`Gene Symbol` = NULL
dt$gene = NULL

#melt the matrix to examine systematically
melted_cnv = reshape2::melt(as.matrix(dt))
head(melted_cnv)
colnames(melted_cnv) = c('ens', 'TCGA_ID', 'value')

#Add the gene names
melted_cnv$ens_sub = gsub("\\..*","",melted_cnv$ens)
melted_cnv$gene_symbol = gene_annots$Approved.symbol[match(melted_cnv$ens_sub, gene_annots$Ensembl.ID.supplied.by.Ensembl.)]

#Lets look at CNVs for APC
APC_muts = melted_cnv[melted_cnv$gene_symbol=='APC',]
length(unique(APC_muts$TCGA_ID))
table(APC_muts$value)

#There are some CNVs for APC, we can also inspect the somatic mutation data in teh same way
#Read in the somatic mutation data
snps = fread("TCGA-COAD.somaticsniper_snv.tsv.gz")
snps[1:5,1:5]
length(unique(snps$gene))

#It looks like there are plently of disruptive SNPS in the APC gene
clock_genes = c('ARNTL', 'CRY1', 'CRY2', 'PER1', 'PER2', 'CLOCK')
clocksnps = snps[snps$gene %in% clock_genes,]
effects = c('stop_gained', 'missense_variant')

APC_snps = clocksnps[clocksnps$effect %in% effects,]  
#APC_snps = snps[snps$gene=='APC',]
table(APC_snps$effect)
length(unique(APC_snps$Sample_ID))

#Contrast these mutations to anotehr gene
snps$gene[1:10]
tp73 = snps[snps$gene=='TP73',]
table(tp73$effect)


apc_loss = APC_snps
#Make sure it contains what we want
table(apc_loss$effect)

#We have the individuals with LOF APC, these data will be used for all other analyses.  First the sequencing data is imported - should look familiar

dt = fread("Colon_cancer_Seq.gz")
row.names(dt) = dt$Ensembl_ID
dt$Ensembl_ID = NULL
tt = dt
tt = as.data.frame(t(as.matrix(tt)))
colnames(tt) = row.names(dt)
colnames(tt)[1:10]
row.names(tt)[1:10]
lcols = row.names(tt)

new_df = unlist(strsplit(lcols, '-'))
new_df = matrix(new_df, ncol=4, byrow=T)
new_df = as.data.frame(new_df)
head(new_df)


colnames(new_df) = c('TCGA', 'Center', 'Patient_ID', 'Tumor_class')
row.names(new_df) = lcols
new_df$c = new_df$Tumor_class
table(new_df$c)
new_df = new_df[grepl('A', new_df$c),]
new_df$cc = gsub('A', '', new_df$c)
row.names(new_df)[1:10]
tt$Patient_ID = new_df$Patient_ID[match(row.names(tt), row.names(new_df))]
tt$Tumor_class = as.numeric(as.character(new_df$cc[match(row.names(tt), row.names(new_df))]))
table(tt$Tumor_class)
tt$tumor_type = tt$Tumor_class > 10
tumors = tt[tt$tumor_type=='FALSE',]
tumors$tumor_type = NULL
tumors$Tumor_class = NULL
tumors$Patient_ID = NULL
tumors_melted = reshape2::melt(as.matrix(tumors))
colnames(tumors_melted) = c('Sample_ID', 'Ens', 'value')
tumors_melted$ens_sub = gsub("\\..*","",tumors_melted$Ens)
tumors_melted$gene_symbol = gene_annots$Approved.symbol[match(tumors_melted$ens_sub, gene_annots$Ensembl.ID.supplied.by.Ensembl.)]

#Initially we can inspect if APC expression changes with the diruptive snps
apc_loss$Sample_ID[1:10]

tumors_melted$apc_status = match(tumors_melted$Sample_ID, apc_loss$Sample_ID,nomatch = 0)
tumors_melted$apc_status = ifelse(tumors_melted$apc_status >0, 'APC_LOF', 'APC_WT')


traits = fread("colon_cancer_phenotypes.gz")
set_diagnosis = as.data.frame(table(traits$tumor_stage.diagnoses))
set_diagnosis$Var1
set_diagnosis$stage_score = paste0(1:15)
set_diagnosis$stage_score = ifelse(set_diagnosis$stage_score=='1', 'NA', paste0(set_diagnosis$stage_score))
set_diagnosis$stage_score = as.numeric(as.character(set_diagnosis$stage_score))
traits$stage_score = set_diagnosis$stage_score[match(traits$tumor_stage.diagnoses, set_diagnosis$Var1)]

head(tumors_melted)
traits$submitter_id.samples[1:10]
stage_plotting = tumors_melted[!duplicated(tumors_melted$Sample_ID),]
stage_plotting$stage_score = traits$stage_score[match(stage_plotting$Sample_ID,traits$submitter_id.samples,)]

table(stage_plotting$stage_score)
stage_plotting = na.omit(stage_plotting)


#It looks as though diagnosis age seems somewhat linear with APC CNV
stage_plotting$age = traits$age_at_diagnosis.diagnoses[match(stage_plotting$Sample_ID,traits$submitter_id.samples,)]

stage_plotting$bmi = traits$bmi.exposures[match(stage_plotting$Sample_ID,traits$submitter_id.samples,)]

stage_plotting = na.omit(stage_plotting)
ggplot(stage_plotting, aes(x=apc_status, y=stage_score, fill=apc_status)) + geom_boxplot(fill=c('darkorchid3', 'darkorange')) + theme_minimal() + ylab('stage score') + ggtitle('Stage Score vs circadian status')

ggplot(stage_plotting, aes(x=apc_status, y=age, fill=apc_status)) + geom_boxplot(fill=c('darkorchid3', 'darkorange')) + theme_minimal() + ylab('diagnosis age') + ggtitle('Diagnosis age vs circadian status')

stage_plotting = stage_plotting[stage_plotting$bmi < 100,]
ggplot(stage_plotting, aes(x=apc_status, y=bmi, fill=apc_status)) + geom_boxplot(fill=c('darkorchid3', 'darkorange')) + theme_minimal() + ylab('BMI') + ggtitle('BMI vs circadian status')


#We use a logistic regression to estimate the significance.  While we cannot predict age of diagnosis using CNV, the specific CNV is more likely than the other
stats = glm(formula= apc_status ~ bmi, data=stage_plotting, family = binomial)
summary(stats)

#Since we know that APC mutations are important (a significant number of mutations are found in colon cancer in general), we can find out which genes might mediate these effecs by performing DE.
#Most of these were performed in the 3-29-20 analysis to generate a counts matrix and sample table.  The only difference is that we will add CNV and LOF mutations as variable to the table
summary(tumors_melted$value)
tumors_melted$scaled_value = scales:::rescale(tumors_melted$value, to = c(0, 3e6))
mode(tumors_melted$scaled_value) = 'integer'
summary(tumors_melted$scaled_value)
tt1 = na.omit(tumors_melted)
counts_matrix = reshape2::dcast(tt1, Sample_ID ~ gene_symbol, fun.aggregate = mean, value.var = 'scaled_value')


counts_matrix = na.omit(counts_matrix)
row.names(counts_matrix) = counts_matrix$Sample_ID
counts_matrix[1:5,1:5]
sample_table = as.data.frame(row.names(counts_matrix))
colnames(sample_table) = c('TCGA_ID')
#sample_table$cnv = APC_muts$value[match(sample_table$TCGA_ID, APC_muts$TCGA_ID)]
sample_table$apc_LOF = match(sample_table$TCGA_ID, apc_loss$Sample_ID,nomatch = 0)
sample_table$apc_LOF = ifelse(sample_table$apc_LOF >0, 'Circadian_Missense', 'Circadian_WT')
row.names(sample_table) = sample_table$TCGA_ID
cm = as.matrix(t(counts_matrix))
mode(cm) = 'integer'
cm =cm[apply(cm, 1, function(x) sum( is.na(x) ))==0,]
#First we will contrast to LOF mutations

dds <- DESeqDataSetFromMatrix(cm, sample_table, formula(~ apc_LOF))
#Its a good idea to filter counts matrices for genes which have low or no counts.  We couly use the following to filter for gene which ahve an fpkm of 5 in at least 10 samples

dim(assay(dds))
keep <- rowSums(counts(dds) >=5) >= 10
dds <- dds[keep,]

#WE perform PCA
res.pca <- PCA(t(assay(dds)),  graph = F)
fviz_pca_ind(res.pca,addEllipses =T, col.ind = dds$apc_LOF, repel = F) + ggtitle('Raw counts PCA')
#Repeat PCA using a variance-stablizing approach
vsd <- vst(dds, blind = FALSE)
res.pca <- PCA(t(assay(vsd)),  graph = F)
fviz_pca_ind(res.pca,addEllipses =T, labels=F,col.ind = dds$apc_LOF, repel = F) + ggtitle('Vst transformation counts PCA')
dds1 = DESeq(dds)
#Check the order of levels.  This means that the log2FC will be up in tumors and down in surrounding
table(dds1$apc_LOF)
res <- results(dds1)
res2 = as.data.frame(res)
#look at numebr of genes
summary(res)
lof_DE = as.data.frame(res)
lof_DE = lof_DE[order(lof_DE$padj, decreasing = F),]

go_terms = read.delim('Human_all_GO_terms.tab.gz')
lof_DE$gene_symbol = row.names(lof_DE)
write.table(lof_DE, file = 'DEGs circ missense.txt', sep = '\t', row.names = F)


#go_terms = go_terms[!grepl('TLE7', go_terms$Gene.names),]

#wnt = go_terms[grepl('Catenin', go_terms$Gene.ontology..biological.process.) | grepl('Wnt', go_terms$Gene.ontology..biological.process.) | grepl('Lrp', go_terms$Gene.ontology..biological.process.),]

#wnt = go_terms[grepl('glycolysis', go_terms$Gene.ontology..biological.process.) | grepl('glucose', go_terms$Gene.ontology..biological.process.) | grepl('warburg', go_terms$Gene.ontology..biological.process.),]

#wnt = go_terms[grepl('pentose', go_terms$Gene.ontology..biological.process.) | grepl('NADP', go_terms$Gene.ontology..biological.process.) | grepl('ribose', go_terms$Gene.ontology..biological.process.),]

#wnt = go_terms[grepl('circadian', go_terms$Gene.ontology..biological.process.) | grepl('rhythm', go_terms$Gene.ontology..biological.process.) ,]

#res1 = res2[row.names(res2) %in% wnt$Gene.names,]
res1 = lof_DE
hist(-log10(res1$pvalue), n = 200)
res1 = na.omit(res1)
res1$label = ifelse(-log10(res1$pvalue) > 42, paste0(row.names(res1)), "")

res1$label_col = ifelse(res1$log2FoldChange > 0, 'dodgerblue4', 'firebrick3')
res1$label_col = ifelse(res1$label =='', '', paste0(res1$label_col))
table(res1$label_col)
ggplot(res1, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point() +
  geom_label_repel(
    aes(x=log2FoldChange, y=-log10(pvalue), label = res1$label),
    color = res1$label_col,
    box.padding = 0.5, point.padding = 0.3,
    segment.color = 'grey50'
  ) + theme_classic() + ggtitle('Volcano Plot Circadian Variants')







lof_DE$gene_symbol = row.names(lof_DE)
write.table(lof_DE, file = 'DEGs circ missense.txt', sep = '\t', row.names = F)


#go_terms = go_terms[!grepl('TLE7', go_terms$Gene.names),]

wnt = go_terms[grepl('Catenin', go_terms$Gene.ontology..biological.process.) | grepl('Wnt', go_terms$Gene.ontology..biological.process.) | grepl('Lrp', go_terms$Gene.ontology..biological.process.),]

wnt2 = go_terms[grepl('glycolysis', go_terms$Gene.ontology..biological.process.) | grepl('glucose', go_terms$Gene.ontology..biological.process.) | grepl('warburg', go_terms$Gene.ontology..biological.process.),]

wnt3 = go_terms[grepl('pentose', go_terms$Gene.ontology..biological.process.) | grepl('NADP', go_terms$Gene.ontology..biological.process.) | grepl('ribose', go_terms$Gene.ontology..biological.process.),]

#wnt = go_terms[grepl('circadian', go_terms$Gene.ontology..biological.process.) | grepl('rhythm', go_terms$Gene.ontology..biological.process.) ,]

#res1 = res2[row.names(res2) %in% wnt$Gene.names,]
circ_genes = read.delim('new_circs.txt')

pathways = as.data.frame(rbind(wnt, wnt2, wnt3))

gene_ids = unlist(strsplit(pathways$Gene.names, ' '))
s <- strsplit(pathways$Gene.names, split = " ")
gene_id = data.frame(Category = rep(pathways$Gene.names, sapply(s, length)), Gene = unlist(s))

new_gene_set = paste0(c(circ_genes$Gene, 'FGF21', 'GAPDHS', 'SLC25A27', 'AXIN2', 'WNT5B', 'APC', 'MYC', 'PRDM1'))

res1 = lof_DE
res1 = res1[res1$gene_symbol %in% new_gene_set,] 
#res1 = res1[abs((res1$log2FoldChange)) < 10,]
#write.table(res1, file = 'pathways-specific DEGs circadian.txt', sep = '\t', row.names = F)
hist(-log10(res1$pvalue), n = 200)
hist(abs(res1$log2FoldChange), n = 200)
res1 = na.omit(res1)
res1 = res1[abs(res1$log2FoldChange) <10,]
res1$label = ifelse(-log10(res1$pvalue) > 1 | res1$log2FoldChange > 0.5 | res1$log2FoldChange < -0.5, paste0(row.names(res1)), "")

res1$label_col = ifelse(res1$log2FoldChange > 0, 'dodgerblue4', 'firebrick3')
res1$label_col = ifelse(res1$label =='', '', paste0(res1$label_col))
table(res1$label_col)
ggplot(res1, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point() +
  geom_label_repel(
    aes(x=log2FoldChange, y=-log10(pvalue), label = res1$label),
    color = res1$label_col,
    box.padding = 0.5, point.padding = 0.3,
    segment.color = 'grey50'
  ) + theme_classic() + ggtitle('Volcano Plot Circadian Variants')


sig_genes_circ = res1$label
length(unique(sig_genes_circ))
spec_genes = tumors_melted[tumors_melted$gene_symbol %in% sig_genes_circ,]

gene_cor = reshape2::dcast(spec_genes, Sample_ID ~ gene_symbol, value.var = 'value', fun.aggregate = mean)
trait_cor = traits %>% select('Sample_ID' = submitter_id.samples, 'BMI' = bmi.exposures, stage_score, age_at_initial_pathologic_diagnosis)

gene_cor = gene_cor[gene_cor$Sample_ID %in% trait_cor$Sample_ID,]
trait_cor = trait_cor[trait_cor$Sample_ID %in% gene_cor$Sample_ID,]
gene_cor = gene_cor[order(as.character(gene_cor$Sample_ID)),]
trait_cor = trait_cor[order(as.character(trait_cor$Sample_ID)),]
gene_cor$Sample_ID[1:10]
trait_cor$Sample_ID[1:10]

trait_cor$Sample_ID = NULL
gene_cor$Sample_ID = NULL

pp1 = bicorAndPvalue(gene_cor, trait_cor, use = 'p')
cc2 = pp1$bicor
tt4 = pp1$p
tt3 = ifelse(tt4 < 0.05,"*","")

pheatmap(cc2, fontsize_number = 20, display_numbers = tt3, cluster_cols = F, number_color = "black", color = colormap(colormap = colormaps$picnic, nshades = 50), main='Correlation between clock DEGS and traits', fontsize_row = 5, fontsize_col = 5)

