import("DESeq2")
import("tximport")
import("pheatmap")
import("vsn")
import("RColorBrewer")
import("ggplot2")
import("reshape2")
import("biomaRt")
import("gProfileR")
import("EnhancedVolcano")
import("grDevices")
import("utils")
import("stats")
import("pathview")
import("fgsea")
import("clusterProfiler")

prepare_df_for_plotting<-function(df,facet,x,y){
  cols<-c(x,y,deparse(substitute(facet)))
  df<-df[which(eval(substitute(facet), df) %in% facet),cols]
  df<-melt(df)
  return(df)
}

bar_plot<-function(df,x,y,facet,labels){
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols<-c(x,y,deparse(substitute(facet)))
  df<-df[which(eval(substitute(facet), df) %in% facet),cols]
  df<-melt(df)
  p<-ggplot(df,aes(eval(parse(text=x)) ,value))+
    geom_col(aes(fill=variable))+
    facet_grid(~eval(substitute(facet), df),scales="free", space="free_x", labeller=labeller(labels))+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=20))+
    ylab('Number of reads')+
    xlab('Sample')+
    labs(fill='Mapped to T muris')
  return(labels)
}

setup_dds_from_kallisto <- function(meta_data,maindir,design,tx2gene){
  files<-file.path(maindir, 'kallisto_counts',meta_data$sample,'abundance.h5')
  names(files) <- meta_data$sample
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds <- DESeqDataSetFromTximport(txi,colData=meta_data, design)
  return(dds)
}


setup_dds_from_counts_table<-function(counts,meta_data,design){
  dds<-DESeqDataSetFromMatrix(counts,colData=meta_data,design)
  return(dds)
}

setup_dds_from_htseq<-function(sample_table,dir,design){
	dds<-DESeqDataSetFromHTSeqCount(sample_table,dir,design)
	return(dds)
}


plot_pca <- function(dds,group,dir){
  vsd<- vst(dds)
  pcaData <- plotPCA(vsd, intgroup=group, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  g<-(ggplot(pcaData, aes(PC1, PC2, color=group)) +
        geom_point(size=3) +
        labs(color = group) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        theme(text = element_text(size=10)) +
        coord_fixed())
  pdf(paste0(dir,'/PCA_',group,'.pdf'))
  print(g)
  dev.off()
  return(g)
}

deseq2_wald <- function(dds){
    dds<-DESeq(dds,test="Wald")
    return(dds)
}

deseq2_lrt <- function(dds,reduced){
  dds<-DESeq(dds,test="LRT",reduced=reduced)
  return(dds)
}

ordered_results<-function(dds){
  res<-results(dds,alpha=0.1)
  res<-res[order(res$padj),]
  return(res)
}

write_results<-function(results,dir,test){
  write.table(results,file=file.path(dir,paste0(test,"_results.tsv")),row.names=F,sep="\t",eol="\n")
  return(results)
}

write_results_text<-function(results,dir,file_name){
  write.table(results,file=file.path(dir,file_name),row.names=F,col.names=F,sep="\t",eol="\n")
  return(results)
}

mart<-function(results,cutoff=NULL){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  if (!is.null(cutoff)){
    values <- row.names(subset(results, padj < cutoff))
  }
  else {
    values <- row.names(results)
  }
  filters <- ('ensembl_gene_id')
  attributes <- c('ensembl_gene_id','external_gene_name','description')
  query <- getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
  results<-merge(query,results,by.x = "ensembl_gene_id",by.y="row.names",all.x=TRUE,all.y=TRUE)
  return(results)
}

ensembl_to_entrez <-function(ensembl_id_list){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  filters <- ('ensembl_gene_id')
  attributes <- c('ensembl_gene_id','entrezgene_id','entrezgene_accession','external_gene_name')
  query <- getBM(attributes = attributes, filters = filters, values = ensembl_id_list, mart = ensembl)
  return(query)
}

genename_to_ensembl <-function(external_gene_names){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  filters <- ('external_gene_name')
  attributes <- c('ensembl_gene_id','external_gene_name')
  query <- getBM(attributes = attributes, filters = filters, values = external_gene_names, mart = ensembl)
  return(query)
}

retrieve_gene_from_go<-function(go_id){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  filters<-('go')
  attributes <- c('ensembl_gene_id')
  values <- go_id
  query<-getBM(attributes = attributes,filters = filters, values = values, mart = ensembl)
  return(query)
}

retrieve_go_descriptions<-function(go_list){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  filters<-('go')
  attributes <- c('go_id','name_1006')
  query<-getBM(attributes = attributes, filters = filters, values = go_list, mart = ensembl)
  return(query)
}

draw_kegg_pathway<-function(pathway,genes,dir){
  setwd(dir)
  pathview(gene.data = genes, 
           pathway.id = pathway, 
           species = "mmu", 
           limit = list(gene=5, cpd=1),
           kegg.dir = dir)
}

plot_enrichment_score<-function(plot_title, msigdbr_list,gene_stats,dir){
  pathway<-msigdbr_list[[plot_title]]
  p<-plotEnrichment(pathway,gene_stats)
  q<-p + ggtitle(plot_title)
  pdf(file.path(dir,paste0(plot_title,".pdf")),10,5)
  print(q)
  dev.off()
  return(p)
}

volcano_plot<-function(results,dir,file_name){
  results<-results[complete.cases(results[,c('padj')]),]
  to_label<-results[which(results$padj<0.001 | abs(results$log2FoldChange) > 2 ),'external_gene_name']
  g<-EnhancedVolcano(results,lab=results$external_gene_name,x='log2FoldChange',y='padj',selectLab=to_label,gridlines.major = FALSE,gridlines.minor = FALSE,pCutoff = 0.05,FCcutoff = 0,xlim=c(-10,10),ylim = c(0, max(-log10(results[,'padj']), na.rm=TRUE) + 1))
 #g<-EnhancedVolcano(results,lab=results$external_gene_name,x='log2FoldChange',y='padj',selectLab=to_label,gridlines.major = FALSE,gridlines.minor = FALSE,pCutoff = 0.05,FCcutoff = 0,xlim=c(-10,10),ylim = c(0,3))
  pdf(paste0(dir,'/',file_name,'.pdf'))
  print(g)
  dev.off()
}

text_matrix <- function(data, table_text) {
  
  rbind(c(table_text, rep('', ncol(data)-1)), # title
        rep('', ncol(data)), # blank spacer row
        names(data), # column names
        unname(sapply(data, as.character))) # data
  
}