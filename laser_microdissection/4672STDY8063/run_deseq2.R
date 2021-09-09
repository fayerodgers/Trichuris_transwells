#!/usr/bin/env Rscript
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

#things to edit
variables<-c('infection_status','organoid_line','fixative') #put the factor that should be removed in the LRT test LAST.
organoid_line<-c('WT1','WT4')
infection_status<-c('control','larvae')
#plate<-factor(c('1','2','NA'))
fixative<-c('CAR','CLA')
#test<-"conditionlarvae"
tx2gene <- read.table('~/git_repos/Trichuris_transwells/transcripts2genes.E98.tsv', header = T)
#

maindir<-'/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection'
subdir<-paste0('fix.',paste(fixative,collapse ="."),'.',paste(organoid_line,collapse ="."),'.',paste(infection_status,collapse ="."))
dir.create(file.path(maindir, subdir))
dir<-file.path(maindir,subdir)

meta_data<-read.table('/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection/meta_data.txt',header=T)
meta_data<-meta_data[which(meta_data$organoid_line %in% organoid_line
                           & meta_data$fixative %in% fixative
                           & meta_data$infection_status %in% infection_status
                           ),]
#meta_data$plate<-factor(meta_data$plate)
meta_data<-droplevels(meta_data)
meta_data$infection_status<-relevel(meta_data$infection_status, ref="control")
write.table(meta_data,paste0(dir,'/meta_data.txt'),sep="\t") #so we know which samples were used in the analysis

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


dds<-m$setup_dds_from_kallisto(meta_data,maindir,as.formula(design_formula),tx2gene)
for (i in design){
  pca<-m$plot_pca(dds,i,dir)
}
dds_wald<-m$deseq2_wald(dds)
results_wald<-m$ordered_results(dds_wald)
results_wald<-m$write_results(results_wald,dir,"Wald")

dds_lrt<-m$deseq2_lrt(dds,as.formula(reduced_formula))
results_lrt<-m$ordered_results(dds_lrt)
results_lrt<-m$mart(results_lrt,0.1)
results_lrt<-m$write_results(results_lrt,dir,"LRT")
volcano<-m$volcano_plot(results_lrt,dir)




