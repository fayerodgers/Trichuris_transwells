#!/usr/bin/env Rscript
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(scales)
library(stringr)
library(goseq)
library("clusterProfiler")
library("pathview")
library("fgsea")
library("GO.db")
library("msigdbr")
library("reshape2")

### paths and metadata ###
main_dir<-'/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection/March2021'
mouse_kallisto_dir<-file.path(main_dir,'mouse_kallisto')
trichuris_kallisto_dir<-file.path(main_dir,'trichuris_kallisto')
metadata<-read.table(file.path(main_dir,'metadata.txt'),header=T,stringsAsFactors = TRUE)
out_dir<-file.path(main_dir,'results')
tx2gene <- read.table('~/git_repos/Trichuris_transwells/transcripts2genes.E98.tsv', header = T)


### QC: import the table of read numbers and counts and plot ####

mouse_counts <- read.table(file.path(main_dir,'mouse_kallisto_counts.txt'),header=FALSE)
names(mouse_counts) <- c("sample","total_reads","n_pseudoaligned_mouse","n_unique_mouse")
trichuris_counts <- read.table(file.path(main_dir,'trichuris_kallisto_counts.txt'),header=FALSE)
names(trichuris_counts) <- c("sample","total_reads","n_pseudoaligned_trichuris","n_unique_trichuris")
all_counts<-merge(mouse_counts,trichuris_counts)
all_counts<-merge(all_counts,metadata)
all_counts<-all_counts[,c("sample","total_reads","n_pseudoaligned_mouse","n_pseudoaligned_trichuris","condition")]
read_counts<-melt(all_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=10), axis.text.x=element_text(angle=45,hjust=1,size=8)) +
  facet_grid(~condition,scales="free",space="free_x")

pdf(file.path(out_dir,'counts.pdf'),30,10)
print(p)
dev.off()

### Analysis ###
# PCA of all samples
design_formula<-"~condition"
dir.all<-file.path(main_dir,'all_data')
dir.create(dir.all)

# write a table so we know which samples were used in the analysis
write.table(metadata,file.path(dir.all,'metadata.txt'),sep="\t")

# setup dds.whole object
dds.all.mouse<-m$setup_dds_from_kallisto(metadata,file.path(main_dir,"mouse"),as.formula(design_formula),tx2gene)

# plot PCAs
pca<-m$plot_pca(dds.all.mouse,'condition',dir.all)
pca<-m$plot_pca(dds.all.mouse,'plate',dir.all)

# select samples to analyse
# remove samples with no/low material detected
samples_to_remove <- c("4672STDY9495748","4672STDY9495753","4672STDY9495856", # D1 control
                       "4672STDY9495721","4672STDY9495722","4672STDY9495723","4672STDY9495724","4672STDY9495727","4672STDY9495728","4672STDY9495729","4672STDY9495732","4672STDY9495734","4672STDY9495735","4672STDY9495736", "4672STDY9495819", "4672STDY9495820", "4672STDY9495821","4672STDY9495823","4672STDY9495824","4672STDY9495825","4672STDY9495826","4672STDY9495828","4672STDY9495830","4672STDY9495817","4672STDY9495834","4672STDY9495835",  # D1 infected
                       "4672STDY9495794","4672STDY9495795","4672STDY9495796","4672STDY9495798","4672STDY9495799","4672STDY9495808","4672STDY9495809","4672STDY9495890","4672STDY9495892","4672STDY9495898","4672STDY9495900","4672STDY9495904", # D3 control
                       "4672STDY9495776","4672STDY9495868","4672STDY9495872","4672STDY9495873","4672STDY9495875","4672STDY9495877","4672STDY9495880","4672STDY9495885","4672STDY9495886","4672STDY9495887","4672STDY9495888","4672STDY9495770","4672STDY9495772","4672STDY9495778","4672STDY9495782","4672STDY9495786","4672STDY9495867","4672STDY9495869","4672STDY9495876","4672STDY9495881","4672STDY9495883") # D3 infected
metadata.qc<-metadata[which(! metadata$sample %in% samples_to_remove),]
metadata.qc<-metadata.qc[which(metadata.qc$larvae_only_control == FALSE & metadata.qc$technical_control == FALSE),]

# replot PCA
design_formula<-"~condition"
dir.qc<-file.path(main_dir,'all_data_QC_pass')
dir.create(dir.qc)

# write a table so we know which samples were used in the analysis
write.table(metadata.qc,file.path(dir.qc,'metadata.txt'),sep="\t")

# setup dds object
dds.all.mouse.qc<-m$setup_dds_from_kallisto(metadata.qc,file.path(main_dir,"mouse"),as.formula(design_formula),tx2gene)

# plot PCAs
pca<-m$plot_pca(dds.all.mouse.qc,'condition',dir.qc)
pca<-m$plot_pca(dds.all.mouse.qc,'plate',dir.qc)

all_counts_for_qc <- all_counts[,c("n_pseudoaligned_mouse","n_pseudoaligned_trichuris","condition","sample")]
colnames(all_counts_for_qc)<-c("M. musculus","T. muris","condition","sample")


#all_counts_for_qc$condition[all_counts_for_qc$condition == "D1.infected"] <- "D1 Infected (n=17)"
#all_counts_for_qc$condition[all_counts_for_qc$condition == "D3.control"] <- "D3 Control (n=28)"
#all_counts_for_qc$condition[all_counts_for_qc$condition == "D1.control"] <- "D1 Control (n=37)"
#all_counts_for_qc$condition[all_counts_for_qc$condition == "D3.infected"] <- "D3 Infected (n=27)"



read_counts.qc <- melt(all_counts_for_qc,id.vars = c("condition","sample"))
read_counts.qc <- read_counts.qc[which(read_counts.qc$sample %in% metadata.qc$sample),]
colnames(read_counts.qc)<- c("condition","sample","Genome","value")

# new facet label names
new.labs<-c("D1 Control (n=37)", "D1 Infected (n=17)", "D3 Control (n=28)", "D3 Infected (n=27)" )
names(new.labs) <- c("D1.control","D1.infected", "D3.control", "D3.infected")

p<-ggplot(data=read_counts.qc, aes(x=Genome,y=value, fill=Genome)) +
  geom_boxplot() +
#  geom_dotplot(binaxis = "y", position="jitter") +
  geom_jitter(position = position_jitter(height = .2, width = .2)) +
  ylab('Aligned reads') +
  xlab('') +
  theme_minimal()+
  theme(axis.title = element_text(size=12), axis.text  = element_text(size=12), axis.text.x=element_text(size=0), legend.text = element_text(face = "italic")) +
  facet_grid(~condition,scales="free",space="free_x", labeller = labeller(condition = new.labs))

pdf(file.path(dir.qc,'counts.pdf'),8,3)
print(p)
dev.off()

# D1 samples
metadata.d1 <- metadata.qc[which(metadata.qc$time == "D1"),]
dir.d1 <- file.path(main_dir,"day1")
dir.create(dir.d1)

# write a table so we know which samples were used in the analysis
write.table(metadata.d1,file.path(dir.d1,'metadata.txt'),sep="\t")

# replot read counts
read_counts.d1 <- read_counts[which(read_counts$sample %in% metadata.d1$sample),]

p<-ggplot(data=read_counts.d1, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=10), axis.text.x=element_text(angle=45,hjust=1,size=8)) +
  facet_grid(~condition,scales="free",space="free_x")

pdf(file.path(dir.d1,'counts.pdf'),30,10)
print(p)
dev.off()

# setup dds object
design_formula<-"~plate + infection_status"
reduced_formula<-"~plate"
dds.d1.mouse<-m$setup_dds_from_kallisto(metadata.d1,file.path(main_dir,"mouse"),as.formula(design_formula),tx2gene)

# plot PCAs
pca<-m$plot_pca(dds.d1.mouse,'infection_status',dir.d1)
pca<-m$plot_pca(dds.d1.mouse,'plate',dir.d1)
pca<-m$plot_pca(dds.d1.mouse,'sample',dir.d1)

# results
dds.d1.mouse<-DESeq(dds.d1.mouse,test="LRT",reduced=as.formula(reduced_formula))
results_d1<-results(dds.d1.mouse,contrast=c("infection_status","infected","control"))
results_d1<-m$mart(results_d1,0.05)
results_d1<-results_d1[order(results_d1$padj),]
results_d1<-m$write_results(results_d1,dir.d1,"LRT")
volcano<-m$volcano_plot(results_d1,dir.d1,'volcano')

# D1 GO
res.d1<-results_d1
tx.lengths<-rowMedians(dds.d1.mouse@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.d1.mouse)

# remove genes that are not expressed
res.d1.expressed<-res.d1[!(res.d1$baseMean == 0),]

sig.d1<-as.integer(res.d1.expressed$padj<0.05 & !is.na(res.d1.expressed$padj))
names(sig.d1)<-res.d1.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.d1) ]

pwf <- nullp(sig.d1, "mm10", "ensGene",bias.data=tx.lengths)

go.d1 <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
go.d1$p.adjust <- p.adjust(go.d1$over_represented_pvalue, method = "BH")
go.d1<-m$write_results(go.d1,dir.d1,"GO.BP")

# d1 by plate


# d3 samples
metadata.d3 <- metadata.qc[which(metadata.qc$time == "D3"),]
dir.d3 <- file.path(main_dir,"day3")
dir.create(dir.d3)

# write a table so we know which samples were used in the analysis
write.table(metadata.d3,file.path(dir.d3,'metadata.txt'),sep="\t")

# replot read counts
read_counts.d3 <- read_counts[which(read_counts$sample %in% metadata.d3$sample),]

p<-ggplot(data=read_counts.d3, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=10), axis.text.x=element_text(angle=45,hjust=1,size=8)) +
  facet_grid(~condition,scales="free",space="free_x")

pdf(file.path(dir.d3,'counts.pdf'),30,10)
print(p)
dev.off()

# setup dds object
design_formula<-"~plate+infection_status"
reduced_formula<-"~plate"
dds.d3.mouse<-m$setup_dds_from_kallisto(metadata.d3,file.path(main_dir,"mouse"),as.formula(design_formula),tx2gene)

# plot PCAs
pca<-m$plot_pca(dds.d3.mouse,'infection_status',dir.d3)
pca<-m$plot_pca(dds.d3.mouse,'plate',dir.d3)
pca<-m$plot_pca(dds.d3.mouse,'sample',dir.d3)

# results
dds.d3.mouse<-DESeq(dds.d3.mouse,test="LRT",reduced=as.formula(reduced_formula))
results_d3<-results(dds.d3.mouse,contrast=c("infection_status","infected","control"))
results_d3<-m$mart(results_d3,0.05)
results_d3<-results_d3[order(results_d3$padj),]
results_d3<-m$write_results(results_d3,dir.d3,"LRT")
volcano<-m$volcano_plot(results_d3,dir.d3,'volcano')

# D3 GO
res.d3<-results_d3
tx.lengths<-rowMedians(dds.d3.mouse@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.d3.mouse)

# remove genes that are not expressed
res.d3.expressed<-res.d3[!(res.d3$baseMean == 0),]

sig.d3<-as.integer(res.d3.expressed$padj<0.05 & !is.na(res.d3.expressed$padj))
names(sig.d3)<-res.d3.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.d3) ]

pwf <- nullp(sig.d3, "mm10", "ensGene",bias.data=tx.lengths)

go.d3 <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
go.d3$p.adjust <- p.adjust(go.d3$over_represented_pvalue, method = "BH")
go.d3<-m$write_results(go.d3,dir.d3,"GO.BP")
