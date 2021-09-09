library('tximport')
library('DESeq2')
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("biomaRt")
library("gProfileR")
library("plyr")
library("gridExtra")
library("EnhancedVolcano")
library("dplyr")
library("scales")
library("cowplot")
library("stringr")
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

###### things to edit ######

#paths and metadata
maindir<-'/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection/run_34491/laser_microdissection/star'
# figures_dir<-'/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection/run_34491/laser_microdissection/star/figures'
metadata<-read.table(file.path(maindir,'../laser_microdissection_manifest.txt'),header=TRUE,stringsAsFactors = TRUE)
counts_dir<-file.path(maindir,'counts_trimmed_trichuris')

# deseq2 design + variables
# select which variables to use for sample selection

# nworms : bulk, single
nworms<-c("single")

# infection_status : infected, control
infection_status<-c("infected")

# fix : Clarks, Carnoys
fix<-c("Clarks")

# time : day.1, day.3
time<-c("day.3")

# technical_control : TRUE, FALSE
technical_control<-c("TRUE","FALSE")

# trichuris_detected : TRUE, FALSE
trichuris_detected<-c("TRUE")

design <- as.formula("~infection_status")
reduced <- as.formula("~1")

file_prefix <- "day3-clarks-single"

#### should be able to run everything below here ####

#select metadata

dir.create(file.path(maindir,file_prefix))

metadata<-metadata[which(metadata$nworms %in% nworms 
                      &  metadata$infection_status %in% infection_status
                      &  metadata$fix %in% fix
                      &  metadata$time %in% time
                      &  (metadata$trichuris_detected %in% trichuris_detected | metadata$infection_status == "control" )
                         ), ]

# add a new column to metadata for condition
metadata$condition<-paste(metadata$time, metadata$nworms,metadata$infection_status, metadata$fix, sep="_")

# metadata[which(metadata$technical_control == "TRUE"),]$condition <- "techincal_control"

# remove bad samples
metadata<-metadata[which(! metadata$sample %in% c("4672STDY8917864","4672STDY8917908","4672STDY8917909","4672STDY8917905","4672STDY8917906","4672STDY8917907","4672STDY8917939","4672STDY8917869","4672STDY8917881","4672STDY8917874")),]

###### QC: import the table of read numbers and counts and plot ######

import_counts_table<-function(sample,dir){
  count_table<-read.table(file.path(dir,paste0(sample,".counts.txt")),header=F,col.names=c('feature','count'))
  return(count_table)
}

select_features_only<-function(sample.df){
  features<-sample.df[grep("^ENSMUSG",rownames(sample.df)),]
  return(features)
}

summary_table<-function(sample.df){
  sample.df["n_features",]<-colSums(sample.df[grep("^ENSMUSG",rownames(sample.df)),])
  summary<-sample.df[ -(grep("^ENSMUSG",rownames(sample.df))),]
  return(summary)
}

# read all counts files
counts<-lapply(metadata$sample,import_counts_table,file.path(maindir,"counts_trimmed_trichuris"))
names(counts)<-metadata$sample

# merge count files
all_counts<-Reduce(function(x,y) merge( x= x, y = y, by = 'feature', all=T), counts)
colnames(all_counts)<-c('feature',names(counts))
row.names(all_counts)<-all_counts$feature
all_counts<-all_counts[-1]

# remove genes with <10 counts in one or more samples

#colnames(all_counts)<-c("Test_1","Test_2","Test_3","Control_1","Control_2","Test_4","Test_5","Test_6","Control_3","Control_4")
#temp<-all_counts[-c(1:5),c(4,5,9,1,2,3)]
# remove genes with <10 counts in one or more samples
#temp<-temp[apply(temp,1,min) > 10, ]

colnames(all_counts)<-c("Control_1","Control_2","Control_3")
temp<-all_counts[-c(1:5),]
write.table(temp,'laser_counts_worm.txt',sep="\t",quote=FALSE)



# get just gene features
features<-select_features_only(all_counts)
# for validating HT-seq- get tables with summed features and non-aligned reads
summary<-summary_table(all_counts)
summary$feature_type <- rownames(summary)
# reformat for ggplot
summary.melt <- melt(summary,id.vars="feature_type")
colnames(summary.melt)<-c("feature_type","sample","reads")
# merge with metadata
summary.melt <- merge(summary.melt,  metadata[,c("sample","condition")])

# plot summary table
p<-ggplot(data=summary.melt, aes(x=sample,y=reads,fill=feature_type)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,file_prefix,paste0(file_prefix,'.htseq_qc.pdf')),50,10)
print(p)
dev.off()


###### make dds object #####
setup_dds_from_htseq<-function(sample_table,dir,design){
  dds<-DESeqDataSetFromHTSeqCount(sample_table,dir,design)
  return(dds)
}
metadata$file<-paste0(metadata$sample,".counts.txt")
#reorder columns for deseq2
metadata<-subset(metadata, select=c(sample,file,time:condition,sample))

# remove technical controls
metadata<-metadata[which(metadata$technical_control == "FALSE"),]

#write a table so we know which samples were used in the analysis
write.table(metadata,file.path(maindir,file_prefix,'metadata.txt'),sep="\t")

dds<-DESeqDataSetFromHTSeqCount(metadata,counts_dir,design)

###### QC: PCA plots  ######

vars<-c(names(colData(dds)))
for (i in vars){
  pca<-m$plot_pca(dds,i,file.path(maindir,file_prefix))
}

##### DE ######

#relevel factors
dds$infection_status <- relevel(dds$infection_status, ref = "control")

#get results - Wald test
dds<-DESeq(dds, test ="LRT", reduced = reduced)
res <- results(dds)
res<-res[order(res$padj),]

#biomart
res<-m$mart(res,0.05)
res<-res[order(res$padj),]

#write and plot results
res<-m$write_results(res,file.path(maindir,file_prefix),"LRT")
volcano<-m$volcano_plot(res,file.path(maindir,file_prefix),"volcano")


