require(ggplot2)
require(data.table)
require(ggpubr)
require(rstatix)
require(tidyr)
require(ggfortify)
require(umap)
require(rlist)
require(dplyr)
require(Matrix)
require(lme4)
require(glmm)
require(ggstats)
require(GGally)
require(caret)
library(e1071)
require(pROC)
require(data.table)
library(forestploter)
require(glmnet)
require(glmnetSE)
require(mplot)
require(selectiveInference)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
require(progress)
library(scattermore)
require(pheatmap)
require(SCGLR)
require(car)
require(corrplot)
require(ggrepel)
require(ggvenn)

pub <- theme(text = element_text(size = 20))
pubrot <- theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

prot <- c(
  "HLA-A*11:01",
  "HLA-A*66:01",
  "HLA-A*32:01",
  "HLA-B*35:02",
  "HLA-B*46:01",
  "HLA-B*57:01",
  "HLA-B*07:02",
  "HLA-B*44:03",
  "HLA-B*57:03",
  "HLA-DRB1*03:02",
  "HLA-DRB1*04:03",
  "HLA-DRB1*08:04",
  "HLA-DRB1*11:01",
  "HLA-DRB1*11:04",
  "HLA-DRB1*13:03",
  "HLA-DRB1*14:01",
  "HLA-DRB1*15:01",
  "HLA-DRB1*15:03",
  "HLA-DRB1*16:02",
  "HLA-DQA1*01:03",
  "HLA-DQA1*02:01",
  "HLA-DQA1*04:01",
  "HLA-DQA1*05:05",
  "HLA-DQB1*03:01",
  "HLA-DQB1*04:02",
  "HLA-DQB1*05:01",
  "HLA-DQB1*05:03",
  "HLA-DQB1*06:01",
  "HLA-DQB1*06:02",
  "HLA-DPB1*04:01",
  "HLA-DPB1*04:02"
)

risk <- c(
  "HLA-A*02:01",
  "HLA-A*24:02",
  "HLA-A*33:03",
  "HLA-B*08:01",
  "HLA-B*15:10",
  "HLA-B*18:01",
  "HLA-B*39:06",
  "HLA-B*54:01",
  "HLA-C*03:03",
  "HLA-C*03:04",
  "HLA-C*15:02",
  "HLA-DRB1*03:01",
  "HLA-DRB1*04:01",
  "HLA-DRB1*04:02",
  "HLA-DRB1*04:04",
  "HLA-DRB1*04:05",
  "HLA-DRB1*09:01",
  "HLA-DPB1*02:02",
  "HLA-DPB1*03:01"
)

#########################################
#########       Setup`      #############
#########################################`
setwd("D://Lab/DirtyRepertoire/")

test <- read.table(".//allstudies.tsv",header = T,sep = "\t")
testsub <- test[!duplicated(test$subject_id),]
testsub <- within(testsub,studyshort <- unlist(lapply(strsplit(testsub$study_title," et"),function(x) x[[1]])))
table(testsub$studyshort)
testsub[testsub$studyshort == "Rawat",]$allele_designation

################################
#########Meta data##############
################################

#Since all of the metadata are in a bunch of different files let's concat and put them together

tab <- read.table(file = ".//ir_2025-06-07_1629_684468e13c815.tsv",header = T,sep = "\t")
head(tab)
unique(tab$study_id)

test <- read.csv(file = ".//Data/IReceptor/Sequences/IR-T1D-02/t1d.tsv",header = T,sep = "\t")
head(test)

and <- read.csv(file = ".//Data/IReceptor/Metadata/anderson-01.tsv",header = T,sep = "\t")
gom <- read.csv(file = ".//Data/IReceptor/Metadata/gomez-01.tsv",header = T,sep = "\t")
mitch <- read.csv(file = ".//Data/IReceptor/Metadata/Daisy-01.tsv",header = T,sep = "\t")
seay <- read.csv(file = ".//Data/IReceptor/Metadata/seay-01.tsv",header = T,sep = "\t")

I01 <- read.csv(file = ".//Data/IReceptor/Metadata/IR-T1D-01.tsv",header = T,sep = "\t")
I02 <- read.table(file = ".//Data/IReceptor/Metadata/IR-T1D-02.tsv",header = T,sep = "\t")
I03 <- read.table(file = ".//Data/IReceptor/Metadata/IR-T1D-03.tsv",header = T,sep = "\t")
tabl <- rbind(and,gom,seay,mitch,I01,I02,I03)

head(tabl)

unique(tabl$biomaterial_provider) #Not a whole lot of difference in environment

rownames(tabl) <- tabl$repertoire_id
case <- rep("control",nrow(tabl))
case[tabl$disease_diagnosis == "type 1 diabetes mellitus"] <- "case"
tabl <- within(tabl,Status <- case)

write.table(tabl,".//Data/IReceptor/Metadata/full.tsv",sep = "\t")
tabl <- read.csv(".//Data/IReceptor/Metadata/full.tsv",header = T,sep = "\t")


################################
#########stats##################
################################

#We need to first use hamming (or leven) distance to cluster/group the CDR3 seqs and then 
folders <- c("anderson_01","gomez_01","mitchell_01","seay_01","IR-T1D-01","IR-T1D-02","IR-T1D-03")
files <- paste(".//Data/IReceptor/Sequences/",folders,"/t1d.tsv",sep="")
allseq <- data.frame(NULL)
stats <- data.frame(NULL)

# for(x in files){
#   print(paste("grabbing",x))
#   tmp <- fread(file = x,header = T,sep = "\t")[,c(10,22,146,151)]
#   # print(dim(tmp))
#   # allseq <- rbind(allseq,tmp)
#   tmp.ls <- split(tmp,tmp$repertoire_id)
#   names <- names(tmp.ls)
#   unique <- unlist(lapply(tmp.ls,nrow))
#   total <- unlist(lapply(tmp.ls,function(y) sum(y$duplicate_count,na.rm = T)))
#   statmp <- data.frame(Subject = names,Unique = unique,Total = total,Study = rep(x,length(names)))
#   stats <- rbind(stats,statmp)
#   rm(statmp)
#   rm(tmp.ls)
#   rm(tmp)
#   gc()
# }

#save the stats
write.table(stats,".//Data/IReceptor/Sequences/All/stats.tsv",sep = "\t")


###############################
####Read back in###############
###############################

# seqs <- read.csv(file = ".//Data/IReceptor/Sequences/All/geneusage.tsv",header = T,sep = "\t",row.names = NULL)
stats <- read.table(".//Data/IReceptor/Sequences/All/stats.tsv",header = T,sep = "\t")
table <- read.csv(".//Data/IReceptor/Metadata/full.tsv",header = T,sep = "\t")


###################Plot general statistics
stats <- within(stats,Study_Name <- gsub("/t1d.tsv","",gsub(".//Data/IReceptor/Sequences/","",stats$Study)))
stats <- within(stats,Subject_ID <- table[rownames(stats),]$subject_id)
stats <- within(stats,Locus <- table[rownames(stats),]$pcr_target_locus)
stats <- within(stats, Study_Short <- table[rownames(stats),]$Study_Short)

###################Create a table with statistics per study
table <- within(table,Study_Short <- unlist(lapply(strsplit(table$study_title,split = " et"),function(x) x[[1]])))
norm <- data.frame(race = c("nr","Black","Asian","Mixed","Black","Hispanic","White","Hispanic","Hispanic","Asian",
                            "Mixed","nr","nr","White","White","White","White","White"),
                   row.names = names(table(table$race)))
table <- within(table,race_norm <- norm[table$race,])
meta <- table[!duplicated(table$subject_id),]
table(meta$Study_Short)
table(table$race_norm)

getStats <- function(df){
  total <- nrow(df)
  diab <- sum(df$Status == "case")
  nondiab <- sum(df$Status == "control")
  mf_ratio <- sum(df$sex == "male")/nrow(df)
  avg_age <- mean(df$age_max)
  sd_age <- sd(df$age_max)
  min_age <- min(df$age_max)
  max_age <- max(df$age_max)
  return(data.frame(total = total, diab = diab, nondiab = nondiab, sex_male = mf_ratio, 
                    age = avg_age, sd_age = sd_age, min_age = min_age,max_age = max_age))
}
studies <- unique(table$Study_Short)
study_stats <- rbind(within(getStats(meta),Study <- "all"),rbindlist(lapply(split(meta,meta$Study_Short),getStats),idcol = "Study"))
write.table(study_stats,file = ".//Results/IRepertoire/Tables/studystats.tsv",sep = "\t")

################Plot time points
order <- unique(table[order(table$Study_Short),]$subject_id)
ggplot(table,aes(x = as.numeric(table$age_max),y = factor(subject_id,levels = order),group = factor(subject_id,levels = order))) +
  geom_point(aes(color = Study_Short)) +
  geom_line(aes(color = Study_Short)) +
  ylab("Subject") +
  xlab("Age (years)") +
  ggtitle("Subject samples by age") +
  facet_wrap(~pcr_target_locus) +
  pub +
  theme(axis.text.y = element_blank())
ggsave(".//Results/IRepertoire/TCRfigs/fig_Agebysample.svg",width = 7,height = 15,units = "in",dpi = 600 )


######################Get tissue, TRA type and cell subset if present
tcode <- data.frame(Tissue = c("Islets","Blood","Spleen","Lymph Node","Lymph Node","Blood","Blood","Blood","Blood"),row.names = unique(table$tissue))
tissue <- tcode[table$tissue,]
cellsub <- table[,c("cell_phenotype","cell_subset")]
cells <- rep(NA,nrow(cellsub))
cells[unique(grep("CD4",cellsub$cell_phenotype),grep("CD4",cellsub$cell_subset))] <- "CD4"
cells[unique(grep("CD8",cellsub$cell_phenotype),grep("CD8",cellsub$cell_subset))] <- "CD8"
type <- data.frame(ID = table$subject_id, tissue = tissue, cells = cells, locus = table$pcr_target_locus,study = table$Study_Short,row.names = rownames(table)) 

ggplot(type,aes(x = tissue, fill = locus)) +
  geom_bar(position = "stack") +
  facet_grid(study~cells) +
  pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(".//Results/IRepertoire/TCRfigs/fig_Sampletypes.svg",width = 12,height = 10,units = "in",dpi = 600 )

########################
#We know that we have multiple files for each patient some representing time points others representing what is the variance per file?
ggplot(stats,aes(x = Subject_ID,y = Unique,color = Study_Name)) +
  geom_boxplot() +
  geom_point() +
  pub +
  theme(axis.text.x = element_blank())

ggplot(stats,aes(x = Total, y = Unique, color = Study_Name)) +
  geom_point() +
  pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################Plot number of unique TCRs per patient in each study
statbypt <- as.data.frame(stats %>%
                            group_by(Subject_ID) %>%
                            summarize(Unique_avg = mean(Unique),Unique_sum = sum(Unique),
                                      Total_avg = mean(Total),Total_sum = sum(Total),Study_Name = Study_Name[1],
                                      Locus = Locus[1],Study_Short = Study_Short[1], .groups = 'drop') )
rownames(statbypt) <- as.character(statbypt$Subject_ID)

ggplot(statbypt,aes(x = Study_Short,y = Unique_avg,color = Study_Short)) +
  geom_boxplot() +
  geom_jitter() +
  ylab("Unique TCR Sequences") +
  pub +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(".//Results/IRepertoire/TCRfigs/fig_TCRperstudy.svg",width = 10,height = 5,units = "in",dpi = 600 )

ggplot(statbypt,aes(x = Study_Name,y = Total_sum,color = Study_Name)) +
  geom_boxplot() +
  geom_jitter() +
  ylab("Total TCR Sequences") +
  pub +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(".//Results/IRepertoire/TCRfigs/fig_totalTCRperstudy.svg",width = 10,height = 5,units = "in",dpi = 600 )

statbyptlocus <- stats %>%
  group_by(Locus) %>%
  summarize(Unique = sum(Unique),Total = sum(Total),Study_Name = Study_Name[1],Locus = Locus[1], .groups = 'drop') 
statbyptlocus #Table with the stats per loci
write.table(statbyptlocus,file = ".//Results/IRepertoire/Tables/TRABstats.tsv",sep="\t")

ggplot(within(statbypt,age <- meta[match(statbypt$Subject_ID,meta$subject_id),]$age_max),aes(x = age,y = Unique_avg, color = Study_Name)) + #Plot by age
  geom_point() +
  scale_y_continuous(limits = c(0,3000000))
#Not really bad bias

#Plot HLA availablility

mhc <- lapply(as.list(meta$allele_designation),function(x) unlist(strsplit(x,", ")))
mhc_num <- as.numeric(unlist(lapply(mhc,length)))


names(mhc) <- meta$subject_id#metadat$subject_id#table$repertoire_id
mhc_ls <- unique(unlist(mhc))
mhc_ls <- mhc_ls[!is.na(mhc_ls)]
mhc_df <- matrix(F,nrow = length(mhc),ncol = length(mhc_ls))
colnames(mhc_df) <- mhc_ls
rownames(mhc_df) <- names(mhc)
for(i in 1:length(mhc)) {
  # for(i in 1:10){#length(mhc)) {
  rowid <- as.character(names(mhc)[i])
  sub <- mhc[[i]]
  for(x in 1:length(sub)){
    allele <- sub[x]
    mhc_df[rowid,colnames(mhc_df) %in% allele] <- T
  }
}
num_mhc <- matrix(as.numeric(mhc_df), ncol = ncol(mhc_df), nrow = nrow(mhc_df))
colnames(num_mhc) <- mhc_ls
rownames(num_mhc) <- names(mhc)
hist(rowSums(num_mhc))
as.numeric(mhc_df)

diab <- meta[meta$Status == "case",]$subject_id
non <- meta[meta$Status == "control",]$subject_id
# diab <- table[table$Status == "case",]$repertoire_id
# non <- table[table$Status == "control",]$repertoire_id
allelerep <- cbind(apply(mhc_df[diab,],2,sum),apply(mhc_df[non,],2,sum))
colnames(allelerep) <- c("Diabetic","Nondiabetic")
allelerep
write.table(allelerep,".//Results/IRepertoire/Tables/AllelicRepstats.tsv",sep = '\t')

# prot <- c("HLA-DRB1*02","HLA-DQB1*05:01","HLA-DRB1*07:01")
# risk <- c("HLA-B*08:01","HLA-B*15:01","HLA-DRB1*03:01","HLA-DRB1*04:01")

prot <- c("HLA-A*03:25","HLA-A*11:311","HLA-A*66:01","HLA-B*07:341","HLA-B*18:176","HLA-B*18:202" ,"HLA-B*39:05","HLA-B*48:01",
          "HLA-B*57:01","HLA-C*08:01","HLA-DOB*01:02","HLA-DOB*01:18N",
          "HLA-DPB1*18:01","HLA-DQA1*01:03","HLA-DQA1*01:04","HLA-DQA1*05:05","HLA-DQB1*03:19","HLA-DQB1*06:01",
          "HLA-DQB1*06:02","HLA-DRB1*03:02","HLA-DRB1*07:01","HLA-DRB1*11:01","HLA-DRB1*11:02","HLA-DRB1*11:04",
          "HLA-DRB1*13:03","HLA-DRB1*14:01","HLA-DRB1*14:54","HLA-DRB1*15:01","HLA-DRB1*15:03","HLA-DRB5*01:01",
          "HLA-W*03:01","MICB*003:01") 
prot <- prot[prot %in% colnames(mhc_df)]
risk <- c("HLA-DRB1*04:01","HLA-DRB1*04:02","HLA-DRB1*04:05","MICA*009:02")  
risk <- risk[risk %in% colnames(mhc_df)]

#plot number with protective or risk or both
mhc_df <- mhc_df[meta$subject_id,]
risk_factor <- rep("Neither",nrow(mhc_df))
risk_factor[apply(mhc_df[,prot],1,any)] <- "Protective"
risk_factor[apply(mhc_df[,risk],1,any)] <- "Risk"
risk_factor[apply(mhc_df[,prot],1,any) & apply(mhc_df[,risk],1,any)] <- "Both"

ggplot(within(meta, risk <- risk_factor),aes(x = Study_Short,fill = factor(risk,levels = c("Both","Risk","Protective","Neither")))) +
  geom_bar(stat = "count") +
  facet_wrap(~Status) +
  xlab("Study name") +
  ylab("Number of patients") +
  pubrot +
  theme(legend.title = element_blank())
ggsave(".//Results/IRepertoire/TCRfigs/fig_Studyrisk.svg",width = 10,height = 5,units = "in",dpi = 600 )



#########################################
#####,, Consolidate files ############
##########################################

#We should get one table for each patient and one table for each sample containing only important variables
sample <- cbind(table[,c("Study_Short","subject_id","sex","age_max","race_norm","cell_number","Status")],
                type[rownames(table),c(2:4)],
                stats[rownames(table),c("Unique","Total")])
metadat <- cbind(meta[,c("Study_Short","subject_id","sex","age_max","race_norm","cell_number","Status")],
                 data.frame(risk = risk_factor))
rownames(metadat) <- metadat$subject_id
alleles <- `dim<-`(as.numeric(mhc_df), dim(mhc_df))
rownames(alleles) <- rownames(mhc_df)
colnames(alleles) <- colnames(mhc_df)
alleles <- alleles[,!colnames(alleles) == "NA."]

# write.table(sample,file = ".//Data/IReceptor/samples.tsv",sep = "\t")
# write.table(metadat,file = ".//Data/IReceptor/meta.tsv",sep = "\t")
# write.table(alleles,file = ".//Data/IReceptor/alleles.tsv",sep = "\t")





