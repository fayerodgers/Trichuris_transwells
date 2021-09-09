library("ggplot2")
library("reshape2")
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

###### things to edit ######
#paths and metadata
maindir<-'/Users/fr7/git_repos/Trichuris_transwells/laser_microdissection/run_34491/laser_microdissection/mapping_QC'
metadata<-read.table(file.path(maindir,'../laser_microdissection_manifest.txt'),header=TRUE,stringsAsFactors = TRUE)

#add a new column for condition
metadata$condition<-paste(metadata$time, metadata$nworms,metadata$infection_status, metadata$fix, sep="_")
metadata[which(metadata$technical_control == "TRUE"),]$condition <- "technical_control"

###### QC: import various tables of read counts and plot ######

# STAR ~ standard params #
star_stats_mouse <-read.table(file.path(maindir,"star_mouse","star_stats_mouse.txt"),header=FALSE)
colnames(star_stats_mouse) <- c('sample','total','mouse.unique','mouse.multi')
star_stats_trichuris <-read.table(file.path(maindir,"star_trichuris","star_stats_trichuris.txt"),header=FALSE)
colnames(star_stats_trichuris) <- c('sample','total','trichuris.unique','trichuris.multi')

# merge
star_stats <- merge(star_stats_mouse, star_stats_trichuris, all = TRUE)
read_counts <- star_stats

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'star_standard.pdf'),50,10)
print(p)
dev.off()

# kallisto ~ standard #
kallisto_stats_mouse <-read.table(file.path(maindir,"kallisto_mouse","kallisto_mouse_counts_2.txt"),header=FALSE)
colnames(kallisto_stats_mouse) <- c('sample','total','mouse.pseudoaligned','mouse.unique')
# new col for multiple pseudoalignments
kallisto_stats_mouse$multi.mouse<-kallisto_stats_mouse$mouse.pseudoaligned - kallisto_stats_mouse$mouse.unique

kallisto_stats_trichuris <-read.table(file.path(maindir,"kallisto_trichuris","kallisto_trichuris_counts_2.txt"),header=FALSE)
colnames(kallisto_stats_trichuris) <- c('sample','total','trichuris.pseudoaligned','trichuris.unique')
kallisto_stats_trichuris$multi.trichuris<-kallisto_stats_trichuris$trichuris.pseudoaligned - kallisto_stats_trichuris$trichuris.unique

# merge
kallisto_stats <- merge(kallisto_stats_mouse, kallisto_stats_trichuris, all = TRUE)
read_counts <- kallisto_stats

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'kallisto_standard.pdf'),50,10)
print(p)
dev.off()

# kallisto ~ trimmed #
kallisto_stats_mouse <-read.table(file.path(maindir,"kallisto_mouse_trimmed","kallisto_mouse_trimmed_counts.txt"),header=FALSE)
colnames(kallisto_stats_mouse) <- c('sample','total','mouse.pseudoaligned','mouse.unique')
# new col for multiple pseudoalignments
kallisto_stats_mouse$multi.mouse<-kallisto_stats_mouse$mouse.pseudoaligned - kallisto_stats_mouse$mouse.unique

kallisto_stats_trichuris <-read.table(file.path(maindir,"kallisto_trichuris_trimmed","kallisto_trichuris_trimmed_counts.txt"),header=FALSE)
colnames(kallisto_stats_trichuris) <- c('sample','total','trichuris.pseudoaligned','trichuris.unique')
kallisto_stats_trichuris$multi.trichuris<-kallisto_stats_trichuris$trichuris.pseudoaligned - kallisto_stats_trichuris$trichuris.unique

# merge
kallisto_stats <- merge(kallisto_stats_mouse, kallisto_stats_trichuris, all = TRUE)
read_counts <- kallisto_stats

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'kallisto_trimmed.pdf'),50,10)
print(p)
dev.off()

# STAR ~ SE #
star_stats_mouse_se <-read.table(file.path(maindir,"star_mouse_se","star_mouse_se_counts.txt"),header=FALSE)
colnames(star_stats_mouse_se) <- c('id','total','mouse.unique','mouse.multi')
temp <- data.frame(do.call('rbind', strsplit(as.character(star_stats_mouse_se$id),'_',fixed=TRUE)))
colnames(temp) <-c('sample','read')
star_stats <- cbind(star_stats_mouse_se,temp)

read_counts <- star_stats[ , -which(names(star_stats) %in% c("id"))]

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition","read"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'star_se.pdf'),50,10)
print(p)
dev.off()
head(read_counts)

# STAR ~ permissive #
star_stats_mouse_new_params <-read.table(file.path(maindir,"star_mouse_new_params","star_mouse_new_params_counts.txt"),header=FALSE)
colnames(star_stats_mouse_new_params) <- c('sample','total','mouse.unique','mouse.multi')

read_counts <- star_stats_mouse_new_params

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'star_new_params.pdf'),50,10)
print(p)
dev.off()

# STAR ~ trimmed #
star_stats_mouse <-read.table(file.path(maindir,"star_mouse_trimmed","star_stats.txt"),header=FALSE)
colnames(star_stats_mouse) <- c('sample','total','mouse.unique','mouse.multi')
star_stats_trichuris <-read.table(file.path(maindir,"star_trichuris_trimmed","star_stats.txt"),header=FALSE)
colnames(star_stats_trichuris) <- c('sample','total','trichuris.unique','trichuris.multi')

# merge
star_stats <- merge(star_stats_mouse, star_stats_trichuris, all = TRUE)
read_counts <- star_stats

# temp
metadata<-metadata[which(metadata$time == "day.3" & metadata$fix == "Carnoys") ,]

# just choose the samples in this experiment
read_counts <- read_counts[which(read_counts$sample %in% metadata$sample),]
read_counts <- merge(read_counts,  metadata[,c("sample","condition")])
read_counts<-melt(read_counts,id.vars = c("sample","condition"))

p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") +    
  xlab('Sample') + 
  ylab('Reads') +
  facet_grid(~condition,scales="free",space="free_x") +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90), strip.text.x = element_text(size = 20))

pdf(file.path(maindir,'star_trimmed_d3_carnoys.pdf'),20,10)
print(p)
dev.off()


