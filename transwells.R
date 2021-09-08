library(RUVSeq)
library(reshape2)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(biomaRt)
library(EnhancedVolcano)

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")

#things to edit
organoid_line<-c('WT3','WT4')
experiment<-factor(c('3'))
condition<-c('control','ES')
test<-"conditionES"
#

subdir<-paste0('expt.',paste(experiment,collapse ="."),'.',paste(organoid_line,collapse ="."),'.',paste(condition,collapse="."))
dir.create(file.path('/Users/fr7/git_repos/Trichuris_transwells', subdir))
dir<-paste0('/Users/fr7/git_repos/Trichuris_transwells/',subdir,'/')

meta_data<-read.table('/Users/fr7/git_repos/Trichuris_transwells/meta_data.txt',header=T)
meta_data<-meta_data[which(meta_data$organoid_line %in% organoid_line
                         & meta_data$experiment %in% experiment
                         & meta_data$condition %in% condition
                         ),]
meta_data$experiment<-factor(meta_data$experiment)
meta_data<-droplevels(meta_data)
meta_data$condition<-relevel(meta_data$condition, ref="control")
write.table(meta_data,paste0(dir,'/meta_data.txt'),sep="\t") #so we know which samples were used in the analysis

vars<-c()

for (var in c('organoid_line','experiment','condition')){
  if (length(eval(parse(text=var))) > 1){
    vars <- c(vars,var)
  }
}

vars<-paste0('~',paste(vars,collapse='+'))

colors<-brewer.pal(3,"Set2")

import_counts_table<-function(sample,dir){
  count_table<-read.table(paste0(dir,sample,'.count'),header=F,col.names=c('feature','count'))
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

make_plots<-function(dir,file_name,set,meta_data){
  
  pdf(paste0(dir, file_name, "_RLE.pdf"))
  plotRLE(set,outline=FALSE,col=colors[meta_data$condition])
  dev.off()
  
  pdf(paste0(dir,file_name,"_PCA_condition.pdf"))
  plotPCA(set,col=colors[meta_data$condition],k=4)
  dev.off()
  
  pdf(paste0(dir,file_name,"_PCA_line.pdf"))
  plotPCA(set,col=colors[meta_data$organoid_line],k=4)
  dev.off()
  
  pdf(paste0(dir,file_name,"_PCA_experiment.pdf"))
  plotPCA(set,col=colors[meta_data$experiment],k=4)
  dev.off()
}

run_edgeR<-function(set,group,design){
  y <- DGEList(counts=counts(set),group=group)
  y<-calcNormFactors(y,method="upperquartile")
  y<-estimateGLMCommonDisp(y,design)
  y<-estimateGLMTagwiseDisp(y, design)
  fit<-glmFit(y,design)
  return(fit)
}

export_results<-function(DGEGLM,dir,name){
  DGEGLM$table$FDR<-p.adjust(DGEGLM$table$PValue,method="BH")
  DGEGLM$table<-DGEGLM$table[order(DGEGLM$table$PValue),]
  write.table(DGEGLM$table,paste0(dir,name,"_results.tsv"),sep="\t")
  pdf(paste0(dir,name,"_MDPlot.pdf"))
  print(plotMD(DGEGLM))
  dev.off()
  return(DGEGLM)
}

#The action

#read all counts files
counts<-lapply(meta_data$sample,import_counts_table,'/Users/fr7/git_repos/Trichuris_transwells/counts/')
names(counts)<-meta_data$sample

#merge count files
all_counts<-Reduce(function(x,y) merge( x= x, y = y, by = 'feature', all=T), counts)
colnames(all_counts)<-c('feature',names(counts))
row.names(all_counts)<-all_counts$feature
all_counts<-all_counts[-1]

#get just gene features
features<-select_features_only(all_counts)
#for validating HT-seq- get tables with summed features and non-aligned reads
summary<-summary_table(all_counts)
row.names(meta_data)<-meta_data$sample
meta_data<-meta_data[-1]

#remove features that aren't expressed
keep <- filterByExpr(features)
features <- features[keep,]

set<-newSeqExpressionSet(as.matrix(features),phenoData = meta_data)
make_plots(dir,"raw_",set,meta_data)
set<-betweenLaneNormalization(set,which="upper")
make_plots(dir,"UQ_",set,meta_data)


#RUVg: estimate unwanted variation using negative control genes.
#first pass DE analysis to identify least differentially expressed genes.
design<-model.matrix(as.formula(vars), data=pData(set))
fit<-run_edgeR(set,meta_data$condition,design)
lrt<-glmLRT(fit,coef=test)
#also export the edgeR results
export_results(lrt,dir,"vanilla")
top<-topTags(lrt,n=nrow(set))$table
empirical<-rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
res<-residuals(fit,type="deviance") #get the residuals for RUVr
set2<-RUVg(set,empirical, k=1)
make_plots(dir,"RUVg",set2,meta_data)
design<-model.matrix(as.formula(paste0(vars, "+W_1"))  , data=pData(set2))
fit<-run_edgeR(set2,meta_data$condition,design)
##check coefficient
lrt_ruvg<-glmLRT(fit,coef=test) ##
lrt_ruvg<-export_results(lrt_ruvg,dir,"RUVg")

#RUVs: estimate unwanted variation using replicate samples
differences<-makeGroups(meta_data$condition)
genes<-rownames(counts(set))
set3<-RUVs(set,genes,k=1,differences)
make_plots(dir,"RUVs",set3,meta_data)
design<-model.matrix(as.formula(paste0(vars,"+W_1")) , data=pData(set3))
fit<-run_edgeR(set3,meta_data$condition,design)
lrt_ruvs<-glmLRT(fit,coef=test) ##
lrt_ruvs<-export_results(lrt_ruvs,dir,"RUVs")


#RUVr: estimate unwanted variation using residuals
set4<-RUVr(set,genes,k=1,res)
make_plots(dir,"RUVr",set4,meta_data)
design<-model.matrix(as.formula(paste0(vars,"+W_1")), data=pData(set4))
fit<-run_edgeR(set4,meta_data$condition,design)
lrt_ruvr<-glmLRT(fit,coef=test) ##
lrt_ruvr<-export_results(lrt_ruvs,dir,"RUVr")

#make volcano plot
filters <- ('ensembl_gene_id')
attributes <- c('ensembl_gene_id','external_gene_name','description')
values<-row.names(lrt_ruvg$table)
query <- getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
res.df<-as.data.frame(lrt_ruvg$table)
res.df<-merge(query,res.df,by.x = "ensembl_gene_id",by.y="row.names",all.x=TRUE)
pdf(paste0(dir,'/volcano.pdf'))
print(EnhancedVolcano(res.df,lab=res.df$external_gene_name,x='logFC',y='FDR',pCutoff = 0.05, gridlines.major = FALSE,gridlines.minor = FALSE,FCcutoff = 0.5))
dev.off()
write.table(res.df[which(res.df$FDR < 0.05),], paste0(dir,"/sig_results_ruvg.tsv"),sep="\t")

