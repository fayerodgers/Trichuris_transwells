#!/usr/bin/env Rscript
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#things to edit
variables<-c('experiment','condition') #put the factor that should be removed in the LRT test LAST.
organoid_line<-c('WT1','WT3','WT4','None')
experiment<-factor(c('4','5','6','7','8','9','10','11'))
condition<-c('coculture','free_larvae','cells_naive')
maindir<-'/Users/fr7/git_repos/Trichuris_transwells/low_input_libraries'
stats<-read.table(file.path(maindir,'alignment_stats.txt'),header=TRUE)
stats$unmapped<-stats$total-stats$unique-stats$multi
stats$percent_mapped<-(stats$unique+stats$multi)*100/stats$total
#

subdir<-paste0('expt.',paste(experiment,collapse ="."),'.',paste(organoid_line,collapse ="."),'.',paste(condition,collapse="."))
dir.create(file.path(maindir, subdir))
dir<-file.path(maindir,subdir)

meta_data<-read.table('/Users/fr7/git_repos/Trichuris_transwells/low_input_libraries/metadata.txt',header=T)
meta_data<-meta_data[which(meta_data$organoid_line %in% organoid_line
                           & meta_data$experiment %in% experiment
                           & meta_data$condition %in% condition
),]

meta_data$experiment<-factor(meta_data$experiment)
meta_data<-droplevels(meta_data)
meta_data$condition<-relevel(meta_data$condition, ref="free_larvae")
write.table(meta_data,paste0(dir,'/meta_data.txt'),sep="\t") #so we know which samples were used in the analysis
stats<-merge(stats,meta_data,by='sample')

labels<-c(
  'free_larvae' = 'Free larvae',
  'cells_naive' = 'Uninfected transwells',
  'coculture' = 'Infected transwells'
)

y<-c('unique','multi','unmapped')
x<-'sample'
stats.1<-m$prepare_df_for_plotting(stats,condition,x,y)
p1<-ggplot(stats.1,aes(sample,value))+
  geom_col(aes(fill=variable))+
  facet_grid(~condition,scales="free", space="free_x", labeller=labeller(condition = labels))+
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=20))+
  ylab('Number of reads')+
  xlab('Sample')+
  labs(fill='Mapped to T muris')
pdf('tmuris_mapping.pdf',30,5)
print(p1)
dev.off()


condition<-c('cells_naive','coculture')
y<-c('percent_mapped')
stats.2<-m$prepare_df_for_plotting(stats,condition,x,y)
p2<-ggplot(stats.2,aes(sample,value))+
  geom_col(aes(fill=variable))+
  facet_grid(~condition,scales="free", space="free_x", labeller=labeller(condition = labels))+
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=20))+
  ylab('% mapped')+
  xlab('Sample')+
  labs(fill='Mapped to T muris')
pdf('tmuris_mapping2.pdf',30,5)
print(p2)
dev.off()

y<-c('unique','multi','unmapped')
stats.3<-m$prepare_df_for_plotting(stats,condition,x,y)
p3<-ggplot(stats.3,aes(sample,value))+
  geom_col(aes(fill=variable))+
  facet_grid(~condition,scales="free", space="free_x", labeller=labeller(condition = labels))+
  scale_fill_manual(values=cbPalette)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=20))+
  ylab('Number of reads')+
  xlab('Sample')+
  labs(fill='Mapped to T muris')+
  scale_y_log10()
pdf('tmuris_mapping3.pdf',30,5)
print(p3)
dev.off()


design<-c()
for (var in variables){
  if (length(eval(parse(text=var))) > 1){
    design <- c(design,var)
  }
}
reduced<-design[-length(design)]
design_formula<-paste0('~',paste(design,collapse='+'))
if (length(reduced) == 0){
  reduced_formula<-"~1"
}else{
  reduced_formula<-paste0('~',paste(reduced,collapse='+'))
}

meta_data$file<-paste0(meta_data$sample,".txt")
meta_data<-meta_data[,c("sample","file","experiment","condition","organoid_line")]
count_dir<-file.path(maindir,"htseq_counts")

dds<-m$setup_dds_from_htseq(meta_data,count_dir,as.formula(design_formula))

for (i in design){
  pca<-m$plot_pca(dds,i,dir)
}
dds_wald<-m$deseq2_wald(dds)
results_wald<-m$ordered_results(dds_wald)
results_wald<-m$write_results(results_wald,dir,"Wald")

dds_lrt<-m$deseq2_lrt(dds,as.formula(reduced_formula))
results_lrt<-m$ordered_results(dds_lrt)
#results_lrt<-m$mart(results_lrt,0.1)
results_lrt<-m$write_results(results_lrt,dir,"LRT")
volcano<-m$volcano_plot(results_lrt,dir)




