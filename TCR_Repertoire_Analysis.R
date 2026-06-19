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
require(vegan)
require(circlize)
require(stringr)
require(RColorBrewer)
require(networkD3)
require(igraph)
require(ggraph)
require(colormap)
library(performance)
library(see)
library(sjPlot)
require(ggheatmap)

##########################################################
########Setup############################################
##########################################################
setwd(dir = "D://Lab/DirtyRepertoire/")

pub <- theme(text = element_text(size = 20))
pubrot <- theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

#Prot and Risk taken from comparison of multipe papers and RNA prot/risk groups
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
riskt <- gsub("[^[:alnum:]]", ".", risk)
prott <- gsub("[^[:alnum:]]", ".", prot)

maprisk <- data.frame(Allele = c(prott,riskt),Color = c(rep("blue",length(prot)),rep("red",length(risk))))
rownames(maprisk) <- maprisk$Allele
removeNA <- function(x) {
  return(x[!is.na(x)])
}

################################################
###### Read in final files #####################
################################################

sample <- read.table(file = ".//Data/IReceptor/samples.tsv",sep = "\t",header = T)
metadat <- read.table(file = ".//Data/IReceptor/meta.tsv",sep = "\t",header = T)
alleles <- read.table(file = ".//Data/IReceptor/alleles.tsv",sep = "\t",header = T)

metadat$risk[metadat$risk == "Protective"] <- "Protected"
metadat$risk[metadat$risk == "Risk"] <- "At-risk"

ggplot(metadat[!metadat$Study_Short %in% c("Culina","Eugster"),],aes(x = Study_Short,fill = Status)) +
  geom_bar(stat = "count") +
  xlab("Study") +
  ylab("Num patients") +
  pubrot


rownames(metadat) <- metadat$subject_id

table(metadat$Status,metadat$risk)



#We should go ahead and drop Eugster patients because they are too little as well as nPOD69 
#so let's drop samples below 10,000 unique TCR detected
sample <- within(sample,sample_id <- rownames(sample))
dim(sample)
sample <- sample[sample$Unique > 10000,]
dim(sample) #We have 825 from 993, but there are only 770 sequencing files so there is probably overlap with dropped samples
#Next add the oldest age and classify as terminal measurement or not

sample <- as.data.frame(rbindlist(lapply(split(sample,sample$subject_id),function(x) within(x,oldest <- max(x$age_max)))))
rownames(sample) <- sample$sample_id
sample <- within(sample,Terminal <- sample$age_max == sample$oldest)
table(sample$Terminal,sample$Status) #So for the most part there is 
# And lastly add the risk category from the metadata
sample <- within(sample,risk <- metadat[sample$subject_id,]$risk)
term_meta <- sample[sample$Terminal == T,]
#Seems for the most part that Terminal samples ages are greater than non-terminal however there is 
ggplot(sample,aes(x = Status,y = age_max)) +
  geom_boxplot() +
  geom_jitter(aes(colour = Status),alpha = 0.7) +
  facet_wrap(~Terminal) +
  ylab("Age") +
  xlab("Status") +
  pub
#Plot the sequencing depth by the number of unique detected
ggplot(sample,aes(x = Total,y = Unique,colour = Study_Short)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, max(sample$Total), by = 10000000)) +
  pub

################################################
######### HLA Allele analysis ##################
################################################

#Allele usage based on study
newallele <- alleles[rowSums(alleles) != 0,]

ls <- split(newallele,paste(metadat[rownames(newallele),]$Study_Short,metadat[rownames(newallele),]$Status,sep="_"))
tbl <- rbindlist(lapply(ls,
                        function(x) as.data.frame(t(as.data.frame(colSums(x))))),idcol = "Comb")
tbl <- tbl[rowSums(tbl[,-1]) != 0,]
tblorder <- names(sort(colSums(tbl[,-1]),decreasing = T))
tbl <- tbl %>% pivot_longer(colnames(tbl)[-1], names_to = "Allele",values_to = "Number")
tbl <- within(tbl, Type <- unlist(lapply(strsplit(tbl$Allele,"\\."), function(x) paste(x[[1]],x[[2]],sep="-"))))
tbl <- within(tbl,Study <- unlist(lapply(strsplit(tbl$Comb,"_"),function(x) x[[1]])))
tbl <- within(tbl,Status <- unlist(lapply(strsplit(tbl$Comb,"_"),function(x) x[[2]])))

ggplot(tbl, aes(x = Number, y = factor(Allele,levels = rev(tblorder)))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("HLA Allele") +
  xlab("Count") +
  facet_grid(Type~Study,scales = "free",space = "free_y") +
  theme(axis.text.y = element_text(size = 10)) +
  pub
ggsave(".//Results/IRepertoire/TCRfigs/Allelebias.svg",width = 10,height = 20,units = "in",dpi = 600 )

agg <- aggregate(Number ~ Allele + Type + Status, data = tbl, FUN = sum)
riskcol <- rep(NA,nrow(agg))
riskcol[agg$Allele %in% riskt] <- "Risk"
riskcol[agg$Allele %in% prott] <- "Prot"
agg <- within(agg,Risk <- riskcol)
agg[agg$Status == "case",]$Number <- agg[agg$Status == "case",]$Number * -1 
agg <- agg[!is.na(agg$Risk),]
sortorder <- rev(order(agg$Risk))
agg <- agg[sortorder,]
agg <- within(agg,riskcolor <- maprisk[agg$Allele,]$Color)

ggplot(rbind(agg,data.frame(Allele = c("dummy","dummy"),Type = c("dummy","dummy"),Status = c("case","control"),
                            Number = c(-100,100),Risk = c("Risk","Risk"),riskcolor = c("red","red"))),
             aes(x = Number,y = factor(Allele,levels = unique(Allele)),fill = Status)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Status,scales="free_x")+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.spacing.x = unit(0, "mm")) +
  ylab("Alleles") +
  xlab ("# of Patients") +
  pub

#Now check alleles completedness
typing <- unlist(lapply(strsplit(colnames(alleles),"\\."), function(x) paste(x[[1]],x[[2]],sep="-")))
typerep <- t(do.call('rbind',lapply(split(as.data.frame(t(alleles)),typing), colSums)))
mean(rowSums(typerep))

##################Any bias for certain HLAs in T1D status?
#Enrichment with a fishers exact test
allele_rep <- data.frame(Allele = colnames(alleles), 
                         ncase = length(case_ids),
                         ncontrol = length(control_ids),
                         repcase = colSums(alleles[case_ids,]),
                         repcontrol = colSums(alleles[control_ids,]))

exacttest <- function(a,b,ap,bp){
  t <- fisher.test(matrix(c(a,ap-a,b,bp-b),2,2),alternative = "two.sided")
  return(data.frame(LFC = log2(a/b),OR = t$estimate,pval = t$p.value))
}
allele_enrich <- rbindlist(apply(allele_rep,1,function(x) exacttest(as.numeric(x[4]),as.numeric(x[5]),
                                         as.numeric(x[2]),as.numeric(x[3]))))

allele_total <- cbind(allele_rep,allele_enrich)
anno <- term_meta[,c(7,13)]
rownames(anno) <- term_meta$subject_id
pheatmap(alleles[order(metadat$Status),rev(order(allele_total$LFC))],
         cluster_rows = F,cluster_cols = T,annotation_row = anno,show_rownames = F,show_colnames = F)

riskv <- rep("None",nrow(allele_total))
riskv[allele_total$Allele %in% riskt] <- "Risk"
riskv[allele_total$Allele %in% prott] <- "Protective"

allele_total <- within(allele_total,risk <- riskv)
allele_tot_norm <- allele_total[!is.infinite(allele_total$OR) & allele_total$OR != 0,]
allele_tot_norm[rev(order(allele_tot_norm$OR)),]

ggplot(metadat,aes(x = risk, fill = Status)) +
  geom_bar() +
  ylab("# of patients") +
  xlab("Allelic risk level") +
  pubrot

###############################################
########Phenotypic summary#####################
###############################################

#Make a summarization table

summariseStudy <- function(x,hla) {
  num <- nrow(x)
  female <- table(x$sex)[1]
  male <- table(x$sex)[2]
  age <- mean(x$age_max)
  agesd <- sd(x$age_max)
  race <- table(x$race_norm)
  race <- t(data.frame(as.numeric(race),row.names = names(race)))
  status <- (table(x$status))
  status <- t(data.frame(as.numeric(status),row.names = names(status)))
  HLA_profiled <- mean(rowSums(hla))
  HLA_sd <- sd(rowSums(hla))
  return(data.frame(Study = x$Study_Short[1],n = num,male = male, female = female,
                    age = age, agesd = agesd,race,status,Num_HLA = HLA_profiled,HLA_sd = HLA_sd))
}

e <- summariseStudy(metadat,typerep)
f <- lapply(split(metadat,metadat$Study_Short),function(x) summariseStudy(metadat[rownames(x),],typerep[rownames(x),]))
write.table(t(as.data.frame(rbind(e,do.call('bind_rows',f)))),file = ".//Results/IRepertoire/Tables/TCRprofileSummary.tsv",sep = "\t")






###############################################
####### Read in clonal sparse matrix###########
###############################################

#On the server, we used the increased resources to reassemble the datasets into a sparse matrix with rows as clones and  
clonal <- readRDS(".//Data/IReceptor/Sequences/All/sparse_repertoire.rds")
col_ids <- colnames(clonal)
clonal <- clonal[nchar(rownames(clonal)) > 6,]
table(sample[col_ids,]$Study_Short)
table(sample$Study_Short)
head(nchar(rownames(clonal)))
clonsums <- colSums(clonal)

sample <- within(sample, ntotal <- clonsums[rownames(sample)])

#Get a representation of how many times a clone is detected in different individuals
clonrep <- clonal
clonrep@x <- clonal@x / clonal@x
clonrow <- rowSums(clonrep,na.rm = T)
repsums <- colSums(clonrep)


#Check the representation of clones across individuals (i.e. how many share a clone)
# clontap <- table(clonrow)
# clontap <- clontap/sum(clontap)
# plot(clontap)
# overrep <- as.numeric(format(c(clontap[1],clontap[2],clontap[3],clontap[4],sum(clontap[5:length(clontap)])),scientific = F))
# names(overrep) <- c("1","2","3","4",">4")
# overrep
# rm(overrep)
# gc()

#Compare the number of shared sequences between individuals
# comb <- combn(1:length(clonsums),2,simplify = T)[,]
# comb <- comb[,sample(ncol(comb),1000)]
# shared <- NULL
# for (x in 1:ncol(comb)){
#   i <- comb[1,x]
#   j <- comb[2,x]
#   shared <- c(shared,sum((clonrep[,i] + clonrep[,j]) == 2) / clonsums[i])
# }
# gghistogram(shared,y = "count")


#Update the sample and metadat dataframes to line up with what data we have from clonal
sample <- sample[col_ids,]
sample <- sample[!is.na(sample$Study_Short),]
dim(sample) #Okay now we are down to 662 samples

#Now save the ids of the case and control patients
casept <- rownames(sample[sample$Status == "case",])
controlpt <- rownames(sample[sample$Status == "control",])

case_term <- rownames(sample[sample$Status == "case" & sample$Terminal,])
case_nonterm <- rownames(sample[sample$Status == "case" & !sample$Terminal,])
ctrl_term <- rownames(sample[sample$Status == "control" & sample$Terminal,])
ctrl_nonterm <- rownames(sample[sample$Status == "control" & !sample$Terminal,])

case_risk <- rownames(sample[sample$Status == "case" & sample$risk == "Risk",])
ctrl_risk <- rownames(sample[sample$Status == "control" & sample$risk == "Risk",])
case_prot <- rownames(sample[sample$Status == "case" & sample$risk == "Protective",])
ctrl_prot <- rownames(sample[sample$Status == "control" & sample$risk == "Protective",])

compsamps <- list(term_risk = sample[sample$risk == "Risk" & sample$Terminal,],
                  term_prot = sample[sample$risk == "Protective" & sample$Terminal,],
                  nterm_risk = sample[sample$risk == "Risk" & !sample$Terminal,],
                  nterm_prot = sample[sample$risk == "Protective" & !sample$Terminal,])

compsamps <- lapply(compsamps,function(x) lapply(split(x,x$Status),rownames))

sub <- rownames(sample) %in% names(clonsums)
summarydf <- rbind(data.frame(Depth = mean(clonsums,na.rm = T), dsd = sd(clonsums,na.rm = T),
                 Coverage = mean(repsums,na.rm = T), csd = sd(repsums,na.rm = T)),
      rbindlist(lapply(split(rownames(sample)[sub],sample$Study_Short[sub]),
                       function(x) data.frame(Depth = mean(clonsums[x],na.rm = T), dsd = sd(clonsums[x],na.rm = T),
                                              Coverage = mean(repsums[x],na.rm = T), csd = sd(repsums[x],na.rm = T)))))
rownames(summarydf) <- c("All",names(split(rownames(sample)[sub],sample$Study_Short[sub])))
write.table(summarydf,".//Results/IRepertoire/Tables/clonesummary.tsv",sep = "\t")





#########################################
######Clonal basic analsyis#############
#########################################

#1. Compare the CDR3 lengths of case vs control 

getCDR3Lengths <- function(x){
  ls <- list()
  for(i in 1:ncol(x)){
    if(i %% 10 == 0) {print(i)}
    name <- colnames(x)[i]
    vec <- x[,i]
    ls[[name]] <- as.data.frame(table(nchar((rownames(x)[vec > 0]))))
    rm(vec)
  }
  return(ls)
  # return(as.data.frame(table(nchar(x))))
}
cdlen_ls <- getCDR3Lengths(clonrep)
gc()
#Normalize the frequency of each list to the total number of sequences detected
for(i in 1:length(cdlen_ls)){
  cdlen_ls[[i]][,2] <- cdlen_ls[[i]][,2]/repsums[i]
}

cdlen <- rbindlist(cdlen_ls,idcol="Sample")
colnames(cdlen) <- c("Sample","Length","Freq")
cdlen_meta <- cdlen[cdlen$Sample %in% rownames(sample),]
cdlen_meta <- within(cdlen_meta,status <- sample[cdlen_meta$Sample,]$Status)
cdlen_meta <- within(cdlen_meta,term <- sample[cdlen_meta$Sample,]$Terminal)
cdlen_meta <- within(cdlen_meta,study <- sample[cdlen_meta$Sample,]$Study_Short)

ggplot(cdlen_meta,aes(x = as.numeric(Length),y = Freq,colour = status)) +
  # stat_summary(fun = "mean",geom = "point",position = "dodge") +
  # stat_summary(fun = "mean",geom = "line",position = "dodge") +
  geom_smooth(alpha = 0.5)+
  # facet_grid(study~term) +
  facet_wrap(~term) +
  ggtitle("CDR3 Lengths") +
  ylab("% representation") +
  xlab("Lengths") +
  pub


#2. Compare the alpha diversity
#For loops through the abundances
calcDiversity <- function(x){
  shannon <- NULL
  for(i in 1:ncol(x)){
    if(i %% 10 == 0){print(i)}
    # rare <- rarefy()
    sum <- sum(x[,i])
    shannon <- c(shannon,vegan::diversity(x[,i]/sum,index = "shannon"))
  }
  return(shannon)
}
diverse <- calcDiversity(clonal)
diverse.df <- data.frame(Sample = colnames(clonal),Shannon = diverse,Status = sample[colnames(clonal),]$Status,
                         term = sample[colnames(clonal),]$Terminal,risk = sample[colnames(clonal),]$risk)
# plot(diverse,clonsums)
ggplot(diverse.df[!is.na(diverse.df$Status),],aes(x = term,y = Shannon,fill = Status)) +
  geom_boxplot() +
  ylab("Shannon Diversity")+
  xlab("Terminal measurment")+
  ggtitle("Shannon diversity of samples") +
  pubrot

ggplot(diverse.df[!is.na(diverse.df$Status),],aes(x = risk,y = Shannon,fill = Status)) +
  geom_boxplot() +
  ylab("Shannon Diversity")+
  xlab("Risk level")+
  ggtitle("Shannon diversity of samples") +
  pubrot

ggplot(diverse.df[!is.na(diverse.df$Status),],aes(x = risk,y = Shannon,fill = Status)) +
  geom_boxplot() +
  ylab("Shannon Diversity")+
  xlab("Risk level")+
  facet_wrap(~term) +
  stat_compare_means(comparisons = list(c("case","control")))+
  ggtitle("Shannon diversity of samples") +
  pubrot

ano <- aov(Shannon~risk * term * Status,data = diverse.df[!is.na(diverse.df$Status),])
summary(ano)
TukeyHSD(ano)



#3. Compare beta diversityChange over time (like for patients will multiple is there a difference in changes in diversity?)

#Since the datasets are very sparse, we will take the strategy in pairwise comparisons to rarify only common ones/normalizing to total common
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
rarefaction <- function(n,dat){
  rowdat <- data.frame(NULL)
  sums <- colSums(dat)
  t <- Sys.time()
  for(i in 1:ncol(dat)){
    print(i)
    name <- colnames(dat)[i]
    prob <- dat[,i]/sums[i]
    samp <- table(sample(rownames(dat),size = n,replace = T, prob = prob))
    cum <- cumsum(sort(samp))
    rows <- data.frame(Name = names(samp),Value = as.numeric(samp),cum = cum,running = 1:length(samp),j = i)
    rowdat <- rbind(rowdat,rows)
  }
  print(paste("Took ", Sys.time() - t," Seconds",sep=""))
  species <- unique(rowdat$Name)
  index <- data.frame(Name = species,index = 1:length(species))
  rownames(index) <- index$Name
  rowdat <- within(rowdat,i <- index[rowdat$Name,]$index)
  sm <- sparseMatrix(i = rowdat$i,j = rowdat$j,x = rowdat$Value)
  colnames(sm) <- colnames(dat)
  rownames(sm) <- species
  #make the rarefaction curve plot
  g <- ggplot(rowdat,aes(x = cum,y = running,group = j))+
    geom_line() +
    ylab("# of CDR3 detected") +
    xlab("# of samples")
  return(list(g = g, sm = sm))
}
bray <- function(x){
  # 1. rarefy
  tmp.r <- rrarefy(x = t(x),sample = min(colSums(x)))
  # 2. calc bray curtis dist
  v <- vegdist(tmp.r,method = "bray")
  #We need to return the value
  return(v)
} #Takes in the full vector and calculates dissimilarities with abundance in mind
getBetaDiversity <- function(x){
  lower <- lower.tri(matrix(NA,nrow = ncol(x),ncol = ncol(x)))
  mat <- matrix(NA,nrow = ncol(x),ncol = ncol(x))
  overlap <- mat
  for (i in 1:ncol(x)){
    tryCatch(print(i))
    for(j in 1:ncol(x)){
      if(lower[i,j]){
        match <- (x[,i] > 0 & x[,j] > 0)
        bray <- tryCatch(bray(x[match,c(i,j)]),error = function(e) {return(NA)})
        mat[i,j] <- bray[1] #jaccard(x[[i]],x[[j]])
        overlap[i,j] <- sum(match)
      }
    }
  }
  return(list(Bray = mat,Overlap = overlap))
}

comps <- unlist(lapply(split(rownames(sample),paste(sample$Study_Short,sample$Status,sample$Terminal)),function(x) sample(x,10,replace=T)))
# tmp <- clonal[(clonal[,1] > 0 & clonal[,2] > 0),c(1:2)]
# 
# v <- 
# rarefaction(10000000,clonal[,sample(1:ncol(clonal),10)])
# test <- lapply(as.list(c(10,100,1000,10000,100000)),function(x) rarefaction(x,clonal[,sample(1:ncol(clonal),3)])) #Test what is the best rarefication (might have to do opposite of rarification)

# rowlist <- list(NULL)
# for(i in 1:ncol(clonal)){
#   if(i %% 10 == 0){print(i)}
#   rowlist[[colnames(clonal)[i]]] <- rownames(clonal)[clonal[,i] > 0]
# }
# rowlist <- rowlist[unlist(lapply(rowlist,length)) > 10000]

beta <- getBetaDiversity(clonal[,comps])#clonal[,clonsums > 10000])

#Get lower distance matrix mirrored for full distance matrix
beta_full <- beta$Bray
beta_full[upper.tri(beta_full)] <- t(beta_full)[upper.tri(beta_full)]
diag(beta_full) <- 0
beta_full[is.na(beta_full)] <- 0
beta_full[1:10,1:10]
# 3. Perform Classical MDS/PCoA
pca_results <- cmdscale(beta_full, k = 2, eig = TRUE)
pca_points <- data.frame(Names = as.character(comps),Dim.1 = pca_results$points[,1], Dim.2 = pca_results$points[,2])
beta_sample <- cbind(sample[pca_points$Names,], pca_points)


ggscatter(beta_sample, x = "Dim.1", y = "Dim.2", 
          # label = rownames(swiss),
          # color = "Study_Short",
          # palette = "jco",
          color = "Status",
          size = 2, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE) + pub + ylab("PCoA Dim 1") + xlab("PCoA Dim 2")
# 4. View results (scores and variance explained)
head(pca_results$points) # Principal component scores

#4. Calculate a model for average abundance to percentage in case vs control based on Term/not term
# case_abundance <- rowSums(clonrep[,casept])/length(casept)
# ctrl_abundance <- rowSums(clonrep[,controlpt])/length(controlpt)
# case_lm <- lm(y~x, data = data.frame(y = case_abundance,x = ctrl_abundance))
# pred <- predict(case_lm,newdata =  data.frame(x = seq(from = 0,to = 1,by = 0.01)),interval = "confidence")
# rm(case_abundance,ctrl_abundance,case_lm)
# gc()

#Actually just do densities
# density(clonal[,1])

status_rep <- data.frame(Case = rowSums(clonrep[,casept])/length(casept),Control = rowSums(clonrep[,controlpt])/length(controlpt))

stripGlmLR <- function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  return(cm)
}
calcSlope <- function(rep,case,ctrl){
  case_rs <- rowSums(rep[,case])
  ctrl_rs <- rowSums(rep[,ctrl])
  df <- data.frame(Case = case_rs/length(case),Control = ctrl_rs/length(ctrl))
  lim <- lm(Case ~ Control,df,model = FALSE, 
            x = FALSE, 
            y = FALSE)
  return(stripGlmLR(lim))
}

lim <- calcSlope(clonrep,casept,controlpt)

split_names <- lapply(split(sample,paste(sample$Terminal,sample$risk)),function(x) lapply(split(x,x$Status),rownames))

split_lim <- lapply(split_names,function(x) calcSlope(clonal,x$case,x$control))
inters <- unlist(lapply(split_lim,function(x) coef(x)[1]))
slopes <- unlist(lapply(split_lim,function(x) coef(x)[2]))

ggplot(data.frame(x = c(1,0,1),y=c(0,1,1)), aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = inters, slope = slopes, 
              color = rep(c("red","green","brown","blue"),2),
              linetype = rep(c("solid","dashed"),each=4), linewidth = 1)

#Plot a sample of status
sam <- status_rep[sample.int(nrow(status_rep),100000),]
plot(Case ~ Control,sam,xlab = "% of Control",ylab = "% of Case")
abline(lim,lwd = 2)

calcEnrichment <- function(a,b,na,nb){
  df <- data.frame(NULL)#list(NULL)
  use <- (a > 5 & b > 5)
  a <- a[use]
  b <- b[use]
  print(paste("Doing ",sum(use)," clones",sep=""))
  for(i in 1:length(a)){
    if (i %% 10000 == 0 ){print(i)}
    mat <- matrix(c(a[i],na - a[i],b[i],nb - b[i]),2,2)
    test <- fisher.test(mat,alternative = "two.sided")
    # print(test$statistic)
    df <- rbind(df,data.frame(a = a[i], b = b[i], na = na, nb = nb, 
                          pval = test$p.value))#,x2 = test$statistic,ares = test$residuals[1,1],bres = test$residuals[1,2]))
  }
  return(df)
}


#Make sums of clonal usage (base of fishers exact test)
casesum <- rowSums(clonrep[,colnames(clonrep) %in% casept])
controlsum <- rowSums(clonrep[,colnames(clonrep) %in% controlpt])
# tst <- matrix(c(casesum[1],length(casept) - casesum[1],controlsum[1],length(controlpt) - controlsum[1]),2,2)

case_enrich <- calcEnrichment(casesum,controlsum,length(casept),length(controlpt))
case_enrich <- within(within(case_enrich, ratio <- a/b ),padj <- pval *450000)
sig <- case_enrich$padj < 0.05 & abs(log2(case_enrich$ratio)) > 1
case_enrich <- within(within(case_enrich,color <- "grey"),label <- NA)
case_enrich$color[sig] <- "red"
case_enrich$label[sig] <- rownames(case_enrich)[sig]


ggplot(case_enrich,aes(x = log2(ratio),y = -log(padj),label = label)) +
  geom_point(colour = case_enrich$color) +
  ylab("Log10 Fisher p-adjusted") +
  xlab("Log2 Fold change patient detection") +
  geom_hline(yintercept = 2.99,lty = "dashed") +
  geom_vline(xintercept = c(-1,1),lty = "dashed") +
  ggrepel::geom_text_repel() +
  pub
rwb <- colorRampPalette(colors = c("red", "white", "blue"))
ggplot(case_enrich[abs(log2(case_enrich$ratio)) > 1,],aes(x = a/na,y = -log(padj),label = label)) +
  geom_point(aes(color = log2(ratio))) +
  ylab("-Log10 Fisher p-adjusted") +
  xlab("% of cases") +
  scale_color_continuous(palette = c("blue","white","red"))+
  geom_hline(yintercept = 2.99,lty = "dashed") +
  # geom_vline(xintercept = c(-1,1),lty = "dashed") +
  ggrepel::geom_text_repel() +
  labs(color = "Log2FC") +
  ggtitle("CDR3 enrichment") +
  pub



#######################################
#######Read in  clone information#####
######################################


#read in the summarized clonal attributes (we need this to get gene usage)
cloneinfo <- read.table(".//Data/IReceptor/Sequences/All/cdr3.tsv",sep="\t",header = T)
cloneinfo <- cloneinfo[!duplicated(cloneinfo$junction_aa),]
rownames(cloneinfo) <- cloneinfo$junction_aa
cloneinfo <- cloneinfo[rownames(clonal),]
#remove ambiguous TCR gene calls by taking the first (highest scoring) alignment
cloneinfo$v_call <- gsub(",.*","",cloneinfo$v_call)
cloneinfo$j_call <- gsub(",.*","",cloneinfo$j_call)


################################
#########get Gene usage#############
################################

#So we can recalculate per TCR gene and we will define clonality as the total number of clones detected containing a TCR gene rather than unique clones (not proportional) but will have to normalize probably

TCR_genes <- c(unique(cloneinfo$v_call),unique(cloneinfo$j_call))
# length(unique(c(grep("or",cloneinfo$v_call),grep("or",cloneinfo$j_call)))) #Apparently there are 25,576,495 clones with ambiguous calls
TCR_genes <- TCR_genes[-grep("IGH",TCR_genes)]
TCR_genes <- TCR_genes[TCR_genes != ""]
TCR_genes <- TCR_genes[!is.na(TCR_genes)]

#Now we want to collect a list of
TCR_gene_rows <- lapply(as.list(TCR_genes),function(x) unique(c(which(cloneinfo$v_call == x),which(cloneinfo$j_call == x))))
TCR_gene_rows <- lapply(TCR_gene_rows,function(x) cloneinfo$junction_aa[x])
# TCR_gene_rows <- lapply(TCR_gene_rows,function(x) which(rownames(clonal) %in% ))

names(TCR_gene_rows) <- TCR_genes
#Now we can run a function on each to get the colsums of

#Which is the faster method
system.time(
  clonal[test,]
) #3.56
system.time(
  which(rownames(clonal) %in% test)
) #2.68
system.time(
  match(test,rownames(clonal))
) #2.27

sumUp <- function(x,clones,sums) {
  sub <- clones[match(x,rownames(clones)),]
  if (is.null(nrow(sub))){
    return(NA)
  }
  if (nrow(sub) == 0){
    return(NA)  
  } else if (nrow(sub) < 2){
    mat <- sub/sums
  } else {
    mat <- colSums(sub,na.rm = T)/sums
  }
  return(mat)
}

GUlist <- lapply(TCR_gene_rows,function(x) sumUp(x,clonal,clonsums))
gene_usage <- do.call(rbind,GUlist)
rownames(gene_usage) <- TCR_genes

#save the gene usage dataframe (since it is small)
write.table(gene_usage,".//Data/IReceptor/Sequences/All/geneusage.tsv",row.names = T,sep = "\t")
gene_usage <- read.table(".//Data/IReceptor/Sequences/All/geneusage.tsv",header = T,sep = "\t",check.names = F)
rm(TCR_gene_rows)
rm(GUlist)

TCRrep <- gene_usage > 0
colSums(TCRrep,na.rm = T)



####################################
#########Get cluster usage##############
####################################

clust.dat <- read.table(file = ".//Data/IReceptor/Sequences/All/GLIPHClusters.txt",header = F,sep = "\t")
clust.dat <- clust.dat[,1:2]
clust.dat <- distinct(clust.dat) #Got rid of about 20 million redundant rows
colnames(clust.dat) <- c("CDR3","Cluster")
length(unique(clust.dat$Cluster)) #2,018,873 unique clusters
#Remove the removed cdr3 quickly
clust.dat <- clust.dat[clust.dat$CDR3 %in% rownames(clonal),]
dim(clust.dat)
clustnum <- table(clust.dat$Cluster)
length(clustnum)
summary(as.numeric(clustnum))
sum(clustnum == 1) #141,492 or 7%
sum(clustnum == 2) #910,931 or 
sum(clustnum == 3) #295,867
hist(clustnum[clustnum > 3]) #670,133 clusters bigger than 3

#Sample how many patients we can detect in based on sample size to pick a minimum cluster size
# meanptnum_ext <- NULL
# sd_s <- NULL
# i_s <- NULL
# for (i in c(0.1,1:50)){
#   print(i)
#   min <- i*1000
#   max <- min + 999
#   use <- clustnum[clustnum > min & clustnum < max]
#   if (length(use) > 3 ){ } else { next }
#   clusters <- sample(as.numeric(names(use)),3)
#   #sample 10 of them each
#   rep <- NULL
#   for (q in 1:3) {
#     cdr3 <- clust.dat[clust.dat$Cluster == clusters[q],]$CDR3
#     if (length(cdr3) == 0) {
#       rep <- c(rep,0)
#     } else if (length(cdr3) == 1) {
#       rep <- c(rep,sum(clonrep[cdr3,]))
#     } else { 
#       rep <- c(rep,sum(colSums(clonrep[cdr3,]) > 1))
#     }
#   }
#   meanptnum_ext <- c(meanptnum_ext,mean(rep))
#   sd_s <- c(sd_s,sd(rep))
#   i_s <- c(i_s,i)
# }
# plot(1:100,meanptnum,main = "Average Patient coverage per cluster size",ylab = "Average # patients covered",xlab = "Number of CDR3 in cluster")
# plot(i_s,meanptnum_ext,main = "Average Patient coverage per cluster size",ylab = "Average # patients covered",xlab = "Number of CDR3 in cluster (1,000's)")


#Get rid of 1-2 sequence clusters (no better than just having the same sequence)
# keepclust <- as.numeric(names(clustnum[clustnum > 10])) #262,874
# keepclust <- as.numeric(names(clustnum[clustnum > 70])) #60,791
keepclust <- as.numeric(names(clustnum[clustnum > 3000])) #1,149 (coverage in at least 90% of patients)

clust.dat <- clust.dat[clust.dat$Cluster %in% keepclust,]
tmp <- data.frame(C = unique(clust.dat$Cluster), New = 1:length(unique(clust.dat$Cluster)))
clust.dat <- within(clust.dat, New_Cluster <- match(clust.dat$Cluster,tmp$C))

subclon <- clonal[clust.dat$CDR3,]

plot(as.numeric(clustnum[keepclust]))
mean(as.numeric(clustnum[keepclust]))
sd(as.numeric(clustnum[keepclust]))

clus_usage <- Matrix(data = 0, nrow = length(unique(clust.dat$Cluster)),ncol = ncol(subclon))
pb <- progress_bar$new(format = "  running [:bar] :percent eta: :eta",total = nrow(tmp))
for (i in 1:nrow(tmp)){
  rows <- which(clust.dat$New_Cluster == i)
  clus_usage[i,] <- colSums(subclon[rows,])
  # pb$tick()
  print(i)
}


#Get a list of which rows of clonal contain each cluster
# Clus_gene_rows <- lapply(split(clust.dat,clust.dat$Cluster),function(x) x[[1]])
# Clus_gene_rows <- lapply(Clus_gene_rows,function(x) x[!is.na(x)])
# CUlist <- lapply(Clus_gene_rows,function(x) sumUp(x,clonal,clonsums))
# clus_usage <- do.call(rbind,CUlist)

rownames(clus_usage) <- tmp$C
colnames(clus_usage) <- colnames(clonal)
# clus_usage <- Matrix(clus_usage,sparse = T)


saveRDS(clus_usage,".//Data/IReceptor/Sequences/All/clusterusage.rds")
clus_usage <- readRDS(".//Data/IReceptor/Sequences/All/clusterusage.rds")
# rm(Clus_gene_rows)
# rm(CUlist)

###############################################
########Get Antigen Specificity################
###############################################
db <- read.csv(".//Data/VDJdb/vdjdb.txt",header = T,sep = "\t")
head(db)
dim(db) #226,494 entries
#Find the T1D related antigen genes 
t1dantigens <- c("INS","GAD","G6P","PTPRN","IA","ICA","CHG","ZNT","IAPP","HSP","IGRP")
humanantigens <- unique(db[db$antigen.species == "HomoSapiens",]$antigen.gene)
humanantigens[grep(paste(t1dantigens,collapse="|"),humanantigens)]
ants <- c("INS","INS-DRiP","PTPRN","ZNT8","IGRP","IAPP","PREINS","GAD65","INSDRIP","GAD2","G6PC2") #based on https://pmc.ncbi.nlm.nih.gov/articles/PMC3312399/
#Match the rows of clonal to the database
matching <- match(rownames(clonal),db$cdr3)
matcher <- match(tcrcdr3$aaSeqCDR3,db$cdr3)
sum(!is.na(matcher))
sum(!is.na(matching)) #Only 37,009 matching records out of  
37009/47999684 * 100 #0.07% of clones have some match
matched <- db[matching[!is.na(matching)],]
clonmatched <- clonal[matched$cdr3,]
#Basic look at the matched set
length(unique(matched$antigen.epitope))
length(unique(matched$antigen.gene)) #Found 245 epitope gene specific
sort(table(matched$antigen.epitope),decreasing = T)
sort(table(matched$antigen.gene),decreasing = T)
sort(table(matched$antigen.species),decreasing = T)
sort(table(matched[matched$antigen.species == "HomoSapiens",]$antigen.gene),decreasing = T) #698 ins entries
#Make a small table based on gene identity
sumUp <- function(x,clones,sums) {
  sub <- clones[x,]
  if (is.null(nrow(sub))){
    return(NA)
  }
  if (nrow(sub) == 0){
    return(NA)  
  } else if (nrow(sub) < 2){
    mat <- sub/sums
  } else {
    mat <- colSums(sub,na.rm = T)/sums
  }
  return(mat)
}
epgene_rows <- lapply(split(matched,matched$antigen.gene),function(x) which(rownames(clonmatched) %in% x$cdr3))
epgenelist <- lapply(epgene_rows,function(x) sumUp(x,clonmatched,clonsums))
ep_usage <- do.call(rbind,epgenelist)
rownames(ep_usage) <- names(epgene_rows)
antgenes <- db[db$antigen.gene %in% rownames(ep_usage),]
antgenes <- antgenes[!duplicated(antgenes$antigen.gene),]
antspecies <- rep("Not Human",nrow(antgenes))
antspecies[antgenes$antigen.species == "HomoSapiens"] <- "Human"
antspecies[antgenes$antigen.species == "HomoSapiens" & antgenes$antigen.gene %in% ants] <- "Human T1D"
antrowdata <- data.frame(Epitope = antgenes$antigen.epitope, Gene = antgenes$antigen.gene, Species = antspecies)
rownames(antrowdata) <- antrowdata$Gene
#Save the dataframe for antigen specificity
write.table(ep_usage,file = ".//Results/IRepertoire/Tables/AntigenGeneUsage.tsv",sep = "\t")
write.table(antrowdata,file = ".//Results/IRepertoire/Tables/AntigenGeneRowdat.tsv",sep = "\t")



################################################
######### Sequence Abundance test ##################
################################################

calcSeqDiff <- function(dat,case,control,present = .5){
  case <- case[case %in% colnames(dat)]
  control <- control[control %in% colnames(dat)]
  in_case <- round(length(case) * present)
  in_control <- round(length(control) * present)
  
  case_arr <- dat[,case]
  control_arr <- dat[,control]
  
  #Calc representation by rowsum
  case_rep <- case_arr > 0
  control_rep <- control_arr > 0
  
  use_case <- rowSums(case_rep) > in_case
  use_control <- rowSums(control_rep) > in_control
  use <- use_case & use_control
  
  print(paste("Found ",sum(use_case)," in ", in_case," case samples",sep=""))
  print(paste("Found ",sum(use_control)," in ", in_control," case samples",sep=""))
  print(paste(" Found ",sum(use)," Total",sep=""))
  
  case_arr <- case_arr[use,]
  control_arr <- control_arr[use,]
  
  df <- data.frame(NULL)
  for (i in 1:nrow(case_arr)){
    x <- removeNA(case_arr[i,])
    y <- removeNA(control_arr[i,])
    gene <- rownames(case_arr)[i]
    fc <- mean(x)/mean(y)
    # sign <- (mean(x) - mean(y))/fc
    lfc <- log2(fc)
    if (var(x) == 0 | var(y) == 0){ pval <- 1 } else {
      test <- t.test(x,y)
      pval <- test$p.value
    }
    row <- data.frame(Gene = gene, 
                      mean_case = mean(x), mean_control = mean(y),
                      sd_case = sd(x),sd_control = sd(y),
                      logfc = lfc,#log2(abs(mean(x))/abs(mean(y))),
                      p.value = pval,
                      padj = pval * nrow(case_arr))
    df <- rbind(df,row)
  }
  return(df)
} 
volcanoPlot <- function(res){
  labels <- res$Gene
  labels[abs(res$logfc) < 1] <- NA
  labels[abs(res$padj) > 0.05] <- NA
  col <- rep("grey",nrow(res))
  col[res$logfc > 1 & res$padj < 0.05] <- "red"
  col[res$logfc < -1 & res$padj < 0.05] <- "blue"
  g <- ggplot(data = res, aes(x = logfc, y = -log10(padj),label = labels)) +
    geom_point(alpha = 0.8, size = 3,colour = col) +
    # scale_color_manual(values = c("brown", "yellow", "orange")) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(max.overlaps = Inf,size = 3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-Log10 P-value") +
    pub
  return(g)
}

#Lets just do differential on samples in greater than 50% of any patients
seqdiff <- calcSeqDiff(clonal,casept,controlpt,present = 0.5)
seqdiff[seqdiff$padj < 0.05,] #None
volcanoPlot(seqdiff) + ggtitle("Overrepresented CDR3")

allids <- lapply(split(sample,paste(sample$Terminal,sample$risk)),function(x) lapply(split(x,x$Status),function(y) rownames(y)))
allseqcomp <- lapply(allids,function(x) calcSeqDiff(clonal,x$case,x$control,present = 0.5))
lapply(allseqcomp,function(x) nrow(x[x$padj < 0.05,]))

ggarrange(plotlist = lapply(allseqcomp,volcanoPlot),nrow = 2, ncol = 4,labels = names(allseqcomp))

allseqcomp[["All"]] <- seqdiff
write.table(rbindlist(allseqcomp,idcol = "Subset"),file = "D:/Lab/DirtyRepertoire/Results/IRepertoire/Tables/Supptable_CDR3.tsv",sep = "\t",quote = F)
###################################
#######Antigen Specificity Differential############
###################################

#Read back in the data
ep_usage <- read.table(".//Results/IRepertoire/Tables/AntigenGeneUsage.tsv",sep = "\t",header = T,check.names = F)
antrowdata <- read.table(".//Results/IRepertoire/Tables/AntigenGeneRowdat.tsv",sep = "\t",header = T)

#Functions to calculate the differential expression based on pt ids of trt vs ctrl
removeNA <- function(x) {
  return(x[!is.na(x)])
}
calcSig <- function(x,y){
  pval <- NULL
  for (i in 1:nrow(x)){
    if (length(removeNA(x[i,])) < 1 | length(removeNA(y[i,])) < 1 ){ pval <- c(pval,NA); next } 
    pval <- c(pval,t.test(na.omit(x[i,]),na.omit(y[i,]))$p.value)
  }
  return(pval)
}
calcDiffAbundance <- function(df,trt,ctrl){
  dat <- data.frame(NULL)
  for(i in 1:nrow(df)){
    name <- rownames(df)[i]
    a <- removeNA(df[i,trt])
    b <- removeNA(df[i,ctrl])
    if(length(a) < 10 | length(b) <10) {
      next
    }
    amean <- mean(a,na.rm = T)
    asd <- sd(a)
    bmean <- mean(b,na.rm = T)
    bsd <- sd(b)
    LFC <- log2(amean/bmean)
    pval <- t.test(a,b)
    dat <- rbind(dat,data.frame(Gene = name,amean = amean,asd = asd,bmean = bmean,bsd = bsd,
                          LFC = LFC,pval = pval$p.value,padj = pval$p.value * nrow(df)))
  }
  return(dat)
}
######Calculate differential expression of the antigen specific using case_control term not term
AntSpecdiff <- calcDiffAbundance(ep_usage,casept,controlpt)
AntSpecdiff <- AntSpecdiff[AntSpecdiff$Gene != "",]
AntSpecdiff
table(AntSpecdiff[AntSpecdiff$pval < 0.05,]$LFC > 0) #20 up 2 down

termASdiff <- calcDiffAbundance(ep_usage,
                                 rownames(sample)[sample$Status == "case" & sample$Terminal],
                                 rownames(sample)[sample$Status == "control" & sample$Terminal])
ntermASdiff <- calcDiffAbundance(ep_usage,
                                 rownames(sample)[sample$Status == "case" & !sample$Terminal],
                                 rownames(sample)[sample$Status == "control" & !sample$Terminal])

table(termASdiff[termASdiff$pval < 0.05,]$LFC > 0) #24 up
table(ntermASdiff[ntermASdiff$pval < 0.05,]$LFC > 0) #6 up 4 down

#######Plots

antvolcanoPlot <- function(res){
  labels <- res$Gene
  labels[abs(res$LFC) < 1] <- NA
  g <- ggplot(data = res, aes(x = LFC, y = -log10(pval),color = Species,label = labels)) +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = c("brown", "yellow", "orange")) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(max.overlaps = Inf) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-Log10 P-value") +
    pub
  return(g)
}

tmp <- within(AntSpecdiff, Species <- antrowdata[AntSpecdiff$Gene,]$Species)
antvolcanoPlot(within(AntSpecdiff, Species <- antrowdata[AntSpecdiff$Gene,]$Species))

write.table(rbindlist(list(All = within(AntSpecdiff,Species <- antrowdata[AntSpecdiff$Gene,]$Species),
                           Terminal = within(termASdiff,Species <- antrowdata[termASdiff$Gene,]$Species),
                           Nonterminal = within(ntermASdiff,Species <- antrowdata[ntermASdiff$Gene,]$Species)),idcol = "Subset"),
            file = ".//Results/IRepertoire/Tables/Supptable_antspec.tsv",sep = "\t",quote = F)

########replot
# ASrelist <- list(Case = AntSpecdiff[,!colnames(AntSpecdiff) %in% c('amean','asd')],
#                  Control = AntSpecdiff[,!colnames(AntSpecdiff) %in% c('bmean','bsd')])
ASrelist <- list(Case = ntermASdiff[,!colnames(ntermASdiff) %in% c('amean','asd')],
                 Control = ntermASdiff[,!colnames(ntermASdiff) %in% c('bmean','bsd')])
ASre <- rbindlist(lapply(ASrelist,function(x) setNames(x,c("Gene","Mean","SD","LFC","pval","padj"))),idcol = "Status")
ggplot(ASre[ASre$Gene %in% ants,],aes(x = Gene,y = Mean,fill = Status)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),position = position_dodge(width = 1),width = 0.5) +
  coord_cartesian(ylim = c(0, 0.0015)) +
  xlab("Epitope Gene") +
  ylab("Epitope target abundance") +
  pubrot


#####################################
#######TCR Differential #########
#####################################

#Function to calculate differential abundance
calcDiff <- function(dat,case,control){
  case <- case[case %in% colnames(dat)]
  control <- control[control %in% colnames(dat)]
  case_arr <- dat[,case]
  control_arr <- dat[,control]
  df <- data.frame(NULL)
  for (i in 1:nrow(case_arr)){
    x <- removeNA(case_arr[i,])
    y <- removeNA(control_arr[i,])
    gene <- rownames(case_arr)[i]
    fc <- mean(x)/mean(y)
    # sign <- (mean(x) - mean(y))/fc
    lfc <- log2(fc)
    if (var(x) == 0 | var(y) == 0){ pval <- 1 } else {
      test <- t.test(x,y)
      pval <- test$p.value
    }
    row <- data.frame(Gene = gene, 
                      mean_case = mean(x), mean_control = mean(y),
                      sd_case = sd(x),sd_control = sd(y),
                      logfc = lfc,#log2(abs(mean(x))/abs(mean(y))),
                      p.value = pval,
                      padj = pval * nrow(case_arr))
    df <- rbind(df,row)
  }
  return(df)
} 
volcanoPlot <- function(res){
  labels <- res$Gene
  labels[abs(res$logfc) < 1] <- NA
  labels[abs(res$padj) > 0.05] <- NA
  col <- rep("grey",nrow(res))
  col[res$logfc > 1 & res$padj < 0.05] <- "red"
  col[res$logfc < -1 & res$padj < 0.05] <- "blue"
  g <- ggplot(data = res, aes(x = logfc, y = -log10(padj),label = labels)) +
    geom_point(alpha = 0.8, size = 3,colour = col) +
    # scale_color_manual(values = c("brown", "yellow", "orange")) +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(max.overlaps = Inf,size = 3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-Log10 P-value") +
    pub
  return(g)
}

#Normalize the TCR gene expression across patients
gene_usage <- read.table(".//Data/IReceptor/Sequences/All/geneusage.tsv",header = T,sep = "\t",check.names = F)
gene_usage <- gene_usage[,rownames(sample)]
boxplot(gene_usage) #Looks like there are some columns that contain no data?
plot(colSums(gene_usage,na.rm = T))
boxplot(t(gene_usage))

#Look at the TCR gene expression 
TCRnorm <- gene_usage[!apply(gene_usage,1,function(x) any(is.na(x))),]#got rid of one gene with NA

plot(colSums(TCRnorm)) #We should get rid of patients that are not normalized
TCRnorm <- TCRnorm[,colSums(TCRnorm) > 1] #770 -> 689
plot(colSums(TCRnorm)) #Removed studies patients that are outliers
boxplot(TCRnorm) #TCR expression relatively normal
boxplot(t(TCRnorm))

####################################Results
TCR_res <- calcDiff(TCRnorm,casept,controlpt)
table(abs(TCR_res[TCR_res$padj < 0.05,]$logfc) > 1)
volcanoPlot(TCR_res)


Term_TR <- calcDiff(TCRnorm,case_term,ctrl_term)
nTerm_TR <- calcDiff(TCRnorm,case_nonterm,ctrl_nonterm)
table(abs(Term_TR[Term_TR$padj < 0.05,]$logfc) > 1)
table(abs(nTerm_TR[nTerm_TR$padj < 0.05,]$logfc) > 1)

ggarrange(plotlist = list(Terminal = volcanoPlot(Term_TR),nonTerminal = volcanoPlot(nTerm_TR)),ncol = 1,labels = c("Term","nTerm"))


riskTR <- calcDiff(TCRnorm,case_risk,ctrl_risk)
protTR <- calcDiff(TCRnorm,case_prot,ctrl_prot)
table(abs(riskTR[riskTR$padj < 0.05,]$logfc) > 1)
table(abs(protTR[protTR$padj < 0.05,]$logfc) > 1)

ggarrange(plotlist = list(Terminal = volcanoPlot(riskTR),nonTerminal = volcanoPlot(protTR)),ncol = 1,labels = c("Risk","Prot"))

TCR_comps <- lapply(compsamps,function(x) calcDiff(TCRnorm,x[[1]],x[[2]]))
lapply(TCR_comps,function(x) table(abs(x[x$padj < 0.05,]$logfc) > 1))
ggarrange(plotlist = lapply(TCR_comps,volcanoPlot),ncol = 2,nrow = 2,labels = names(TCR_comps))

tcrdiff_list <- list(All = TCR_res, 
     Term = Term_TR, 
     nTerm = nTerm_TR,
     Risk = riskTR,
     Prot = protTR)
lapply(tcrdiff_list,function(x) dim(x[x$padj < 0.05,]))

TCRdiff_res <- rbindlist(tcrdiff_list,idcol = "Comparison")

write.table(TCRdiff_res,file = ".//Results/IRepertoire/Tables/SupptableDET.tsv",quote = F,sep = "\t")


genes <- lapply(tcrdiff_list,function(x) dim(x[x$padj < 0.05,]))
#Make the umaps for all of these
#######################Prep the TCR features
scaleTCR <- apply(TCRnorm,2,scale)
rownames(scaleTCR) <- rownames(TCRnorm)
boxplot(scaleTCR) #Check that scaleTCR is correctly normalized (remove outliers)
# plot(apply(scaleTCR,2,median),col = c("red","blue")[d.df[colnames(scaleTCR),]$response + 1])
# abline(h = c(-0.24,-0.36), col = "red", lty = 2)
outliers <- which(apply(scaleTCR,2,median) > -0.24 | apply(scaleTCR,2,median) < -0.37)
norm_scaleTCR <- scaleTCR[,-outliers]


allumap <- umap(t(norm_scaleTCR))
allumap <- allumap$layout
colnames(allumap) <- c("UMAP1","UMAP2")
allumap <- cbind(allumap,sample[rownames(allumap),])
a <- ggplot(allumap,aes(x = UMAP1,y = UMAP2,colour = Study_Short)) +
  geom_point() +
  ggtitle("All patients") +
  pub +
  guides(colour = F)
a
riskpt <- rownames(sample[sample$risk == "Risk",])
allumap <- umap(t(norm_scaleTCR[,colnames(norm_scaleTCR) %in% riskpt]))
allumap <- allumap$layout
colnames(allumap) <- c("UMAP1","UMAP2")
allumap <- cbind(allumap,sample[rownames(allumap),])
b <- ggplot(allumap,aes(x = UMAP1,y = UMAP2,colour = Status)) +
  geom_point() +
  ggtitle("Risk patients") +
  pub +
  guides(colour = F) +
  ylab('')
b
protpt <- rownames(sample[sample$risk == "Protective",])
allumap <- umap(t(norm_scaleTCR[,colnames(norm_scaleTCR) %in% protpt]))
allumap <- allumap$layout
colnames(allumap) <- c("UMAP1","UMAP2")
allumap <- cbind(allumap,sample[rownames(allumap),])
c <- ggplot(allumap,aes(x = UMAP1,y = UMAP2,colour = Status)) +
  geom_point() +
  ggtitle("Protective patients") +
  pub + 
  ylab('')
c
ggarrange(plotlist = list(a,b,c),ncol = 3,widths = c(.5,.5,.75))


#############Explore the scaled data
um <- umap(t(norm_scaleTCR))
pc <- prcomp(t(norm_scaleTCR))

# subsample_ex <- cbind(subsample[colnames(norm_scaleTCR),], setnames(as.data.frame(pc$x[,1:2]),new = c("UMAP1","UMAP2")))
um <- umap(pc$x[,1:12])
subsample_ex <- cbind(subsample[rownames(um$layout),], setnames(as.data.frame(um$layout),new = c("UMAP1","UMAP2")))
ggplot(subsample_ex,aes(x = UMAP1, y = UMAP2, colour = Status))+
  geom_point() + facet_wrap(~Study_Short)


##########################################
######  Cluster Usage Differential#########
##########################################

clus_usage <- readRDS(".//Data/IReceptor/Sequences/All/clusterusage.rds")
clus_usage <- clus_usage[,colnames(clus_usage) %in% rownames(sample)]

###############Normalize the cluster usage
normFunc <- function(x) {
  return(x / max(x, na.rm = TRUE)) # na.rm = TRUE handles potential missing values
} #Function to normalize cluster usage values across samples
plot(colSums(clus_usage)) #again we see a study batch effect so let's normalize by average
clusNorm <- apply(clus_usage,2,normFunc)
plot(colSums(clusNorm),xlab = "",ylab = "Sum Normalized Expression",main = "Normalized clonality of clusters per patient") #Normalized roughly
abline(h = 7,lty = 2,lwd = 2)
clusNorm <- clusNorm[,colSums(clusNorm) > 7]

################DIfferential cluster usage
clusres <- calcDiff(clusNorm,casept[casept %in% colnames(clusNorm)],controlpt[controlpt %in% colnames(clusNorm)])
# clussig <- clusres[abs(clusres$logfc) > 1 & clusres$padj < 0.05,]
clussig <- clusres[clusres$padj < 0.05,]
plot(clusres$mean_case,clusres$mean_control,xlab = "Mean Case",ylab = "Mean Control",main = "Diabetic vs Non-diabetic Cluster expression") #Not really anything

write.table(clusres,".//Results/IRepertoire/Tables/supptable5DEC.tsv",sep = "\t")

volcanoPlot(clusres)


termCR <- calcDiff(clusNorm,case_term,ctrl_term)
ntermCR <- calcDiff(clusNorm,case_nonterm,ctrl_nonterm)
riskCR <- calcDiff(clusNorm,case_risk,ctrl_risk)
protCR <- calcDiff(clusNorm,case_prot,ctrl_prot)

table(abs(termCR[termCR$padj < 0.05,]$logfc) > 1)
table(abs(ntermCR[termCR$padj < 0.05,]$logfc) > 1)

ggarrange(plotlist = list(Terminal = volcanoPlot(termCR),nonTerminal = volcanoPlot(ntermCR)),ncol = 1)

table(abs(riskCR[riskCR$padj < 0.05,]$logfc) > 1)
table(abs(protCR[protCR$padj < 0.05,]$logfc) > 1)

ggarrange(plotlist = list(Terminal = volcanoPlot(riskCR),nonTerminal = volcanoPlot(protCR)),ncol = 1)

clus_list <- list(All = clusres,Terminal = termCR, nonterminal = ntermCR,atRisk = riskCR, Protected = protCR)

write.table(rbindlist(clus_list,idcol = "Subset"),file = ".//Results/IRepertoire/Tables/supptable_cluster.tsv",sep="\t",quote = F)





##########################################
########        Functions      ############
##########################################

rotate <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

stat.test <- function(x) {
  x %>%
  group_by(gene) %>%
  t_test(clones ~ Status) %>%
  # adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  add_x_position(x = "gene", dodge = 0.8) %>%
  add_y_position(fun = "mean_ci")
}

repd <- function(x){
  return(names(table(x[,1]))[table(x[,1])>1])
}

plotGenes <- function(meta,gene,mhc,allele,agg = "clones",normalize = F,name=NA){
  rownames(meta) <- meta$repertoire_id
  if(is.na(name)) {
    name <- allele
  }
  if (any(allele %in% "all")) {
    ids <- rownames(mhc)
  } else {
    ids <- rownames(mhc)[apply(mhc[,allele],1,any)]
  }
  match <- meta[ids,]
  cols <- c("Status","gene","Subject","clones")
  case <- match[match$Status == "case",]#ids[ids %in% rownames(meta[meta$Status == "case",])]
  control <- match[match$Status == "control",]#ids[ids %in% rownames(meta[meta$Status == "control",])]
  locus <- unique(gene$locus)
  ########First do TRB##################################################
  if (any("TRB" %in% locus)) {
    trb <- gene[gene$locus == "TRB" & !is.na(gene$d_gene) & !gene$d_gene == "",]#gene[(gene$locus == "TRB"),]#
    trbcase <-  trb[trb$repertoire == rownames(case),]
    trbcontrol <- trb[trb$repertoire == rownames(control),]
    if(nrow(trbcase) > 0) { 
      if(agg == "clones"){
        casev <- aggregate(clones ~ v_gene + subject, data = trbcase, FUN = sum)
        casej <- aggregate(clones ~ j_gene + subject, data = trbcase, FUN = sum)
        cased <- aggregate(clones ~ d_gene + subject, data = trbcase, FUN = sum)
        casenorm <- aggregate(clones ~ subject, data = trbcase, FUN = sum)
      } else {
        casev <- aggregate(Count ~ v_gene + subject, data = trbcase, FUN = sum)
        casej <- aggregate(Count ~ j_gene + subject, data = trbcase, FUN = sum)
        cased <- aggregate(Count ~ d_gene + subject, data = trbcase, FUN = sum)
        casenorm <- aggregate(Count ~ subject, data = trbcase, FUN = sum)
      }
    } else {
      casev <- NULL
      casej <- NULL
      cased <- NULL
    }
    if(nrow(trbcontrol) > 0) { 
      if(agg == "clones"){
        conv <- aggregate(clones ~ v_gene + subject, data = trbcontrol, FUN = sum)
        conj <- aggregate(clones ~ j_gene + subject, data = trbcontrol, FUN = sum)
        cond <- aggregate(clones ~ d_gene + subject, data = trbcontrol, FUN = sum)
        connorm <- aggregate(clones ~ subject, data = trbcontrol, FUN = sum)
      } else {
        conv <- aggregate(Count ~ v_gene + subject, data = trbcontrol, FUN = sum)
        conj <- aggregate(Count ~ j_gene + subject, data = trbcontrol, FUN = sum)
        cond <- aggregate(Count ~ d_gene + subject, data = trbcontrol, FUN = sum)
        connorm <- aggregate(Count ~ subject, data = trbcontrol, FUN = sum)
      }
    } else {
      conv <- NULL
      conj <- NULL
      cond <- NULL
    }
    #Now combine the case vs control
    vint <- intersect(repd(casev),repd(conv))
    trbv <- rbindlist(list(Case = casev[casev$v_gene %in% vint,],Control = conv[conv$v_gene %in% vint,]),idcol = "Status")
    vint <- intersect(repd(casej),repd(conj))
    trbj <- rbindlist(list(Case = casej[casej$j_gene %in% vint,],Control = conj[conj$j_gene %in% vint,]),idcol = "Status")
    vint <- intersect(repd(cased),repd(cond))
    trbd <- rbindlist(list(Case = cased[cased$d_gene %in% vint,],Control = cond[cond$d_gene %in% vint,]),idcol = "Status")
    colnames(trbv) <- cols
    colnames(trbj) <- cols
    colnames(trbd) <- cols
    #normalize 
    label <- "Num Clones"
    if (normalize){
      norm <- rbind(casenorm,connorm)
      # print(norm)
      rownames(norm) <- norm$subject
      trbv$clones <- trbv$clones/(norm[trbv$Subject,]$clones)
      trbd$clones <- trbd$clones/(norm[trbd$Subject,]$clones)
      trbj$clones <- trbj$clones/(norm[trbj$Subject,]$clones)
      label <- "% Clones"
    }
    #Make plots
    Stat <- data.frame(gene = NA,y.position = NA, p.signif = NA)
    try(Stat <- stat.test(trbv))
    a <- ggplot(trbv,aes(x = gene,y = clones,fill = Status)) +
      stat_summary(fun = mean ,geom = "bar",position = "dodge") +
      stat_summary(fun.data = "mean_se",position = position_dodge(width = 0.9), fun.args = list(mult = 1), geom = "errorbar", width = 0.5) +
      annotate(geom = "text",x = Stat$gene, label = Stat$p.signif,y = Stat$y.position*0.75,) +
      # stat_compare_means(comparisons = c("Case","Control")) +
      pubrot +
      ggtitle("TRBV")+
      ylab(label) +
      guides(fill = F)
    try(Stat <- stat.test(trbj))
    c <- ggplot(trbj,aes(x = gene,y = clones,fill = Status)) +
      stat_summary(fun = mean ,geom = "bar",position = "dodge") +
      stat_summary(fun.data = "mean_se",position = position_dodge(width = 0.9), fun.args = list(mult = 1), geom = "errorbar", width = 0.5) +
      annotate(geom = "text",x = Stat$x,y = Stat$y.position*0.75, label = Stat$p.signif) +
      # stat_compare_means(comparisons = c("Case","Control")) +
      pubrot +
      ggtitle("TRBJ")+
      ylab(label) +
      guides(fill = F)
    Stat <- stat.test(trbd)
    e <- ggplot(trbd,aes(x = gene,y = clones,fill = Status)) +
      stat_summary(fun = mean ,geom = "bar",position = "dodge") +
      stat_summary(fun.data = "mean_se",position = position_dodge(width = 0.9), fun.args = list(mult = 1), geom = "errorbar", width = 0.5) +
      # geom_jitter(colour = "black",position = position_jitterdodge(jitter.width = 0.1)) +
      annotate(geom = "text",x = Stat$x,y = Stat$y.position*0.75, label = Stat$p.signif) +
      # stat_compare_means(comparisons = c("Case","Control"), lable = "p.signif") +
      pubrot +
      ggtitle("TRBD")  +
      ylab(label) 
  } else {
    a <- NULL
    c <- NULL
    e <- NULL
  }

  ##########Now create a principle component analysis to see how each allele differs
  tr <- rbind(trbv,trbj,trbd)#,trav,traj)
  tr <- tr[!tr$gene == "",]
  tr <- tr %>% pivot_wider(names_from = gene,values_from = clones,values_fill = 0)
  pc <- prcomp(tr[,c(3:ncol(tr))])#,scale. = T)
  f <- autoplot(pc,tr,colour = "Status") + guides(colour = F) + ggtitle("PCA")
  um <- umap(tr[,c(3:ncol(tr))])
  um <- um$layout
  colnames(um) <- c("x1","x2")
  g <- ggplot(cbind(tr[,c(1:2)],um),aes(x = x1,y = x2,color = Status)) +
          geom_point() +
          guides(colour = F) +
          ggtitle("UMAP")
  ####Plot it out
  # arr <- ggarrange(plotlist = list(a,b,c,d,e,f),ncol = 2, nrow = 3,widths = c(2,1))
  arr <- ggarrange(plotlist = list(a,c,ggarrange(plotlist = list(e,f,g),nrow = 1, ncol = 3)),ncol = 1, nrow = 3)#,widths = c(2,1))
  ann <- annotate_figure(arr, top = text_grob(paste("TCR usage for ",agg, " in ",name,sep = ""), color = "red", face = "bold", size = 14))
  print(ann)
  ls <- list(a,c,e,f,g,ann,tr)
  return(ls)
}

#Make a function that will take a model and plot the coefficients of the model/extract the coefficient and CI table and signficance

getCoefs <- function(l){
  best_coefs <- coef(l)#, s = lim$lambda.min) #min lambda
  coef_df <- data.frame(Term = rownames(best_coefs),Coef = as.numeric(best_coefs))
  coef_df <- coef_df[rev(order(coef_df$Coef)),]
  coef_df$Term <- factor(coef_df$Term,levels = coef_df$Term)
  return(coef_df)
}

get_ROC <- function(x,y){
  proc_obj <- roc(predictor = x,response = y,ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE,quiet = TRUE)
  net_roc <- data.frame(y = rev(proc_obj$sensitivities),x = rev(proc_obj$specificities))
  return(net_roc)
}
get_AUC <- function(y,x){
  if (length(levels(factor(y))) < 2){
    return(NA)
  }
  proc_obj <- roc(predictor = x,response = y,ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE,quiet= TRUE)
  net_auc <- auc(proc_obj)
  return(net_auc)
}
tjurs <- function(y,pred){ #Uses the Tjurs psuedo R2 to calculate an error
  probs_failure <- pred[y == 0]
  probs_failure[is.na(probs_failure)] <- 1
  probs_success <- pred[y == 1]
  probs_success[is.na(probs_success)] <- 0
  mean_failure <- mean(probs_failure,na.rm = T)
  mean_success <- mean(probs_success,na.rm = T)
  if (is.nan(mean_failure)){ mean_failure <- 0 }
  if (is.nan(mean_success)){ mean_success <- 0 }
  return(mean_success - mean_failure)
}
zap <- function(x, y){
  stopifnot(length(x) == length(y))
  z <- list()
  for (i in seq_along(x)){
    z[[i]] <- list(x[[i]], y[[i]])
  }
  return(z)
}
errorRate <- function(y,pred){
  if(length(unique(y)) < 2){ return(list(sensitivity = NA, specificity = NA,FPR = NA,FNR = NA))}
  if(length(unique(pred > 0.5)) < 2){ return(list(sensitivity = NA, specificity = NA,FPR = NA,FNR = NA)) }
  cm_caret <- confusionMatrix(as.factor(as.numeric(pred>0.5)), as.factor(y), positive = "1")
  sensitivity <- cm_caret$byClass["Sensitivity"]
  specificity <- cm_caret$byClass["Specificity"]
  FPR <- 1 - specificity
  FNR <- 1 - sensitivity
  return(list(sensitivity = sensitivity, specificity = specificity,FPR = FPR,FNR = FNR))
}
calcHLATCR <- function(prob,test){
  hlas <- grep("HLA",colnames(test),value = T)
  tcr <- grep("TR",colnames(test),value = T)
  mat <- matrix(data = NA,nrow = length(hlas),ncol = length(tcr))
  for(i in 1:length(hlas)){
    idx <- which(test[,hlas[i]] > 0)
    for(j in 1:length(tcr)){
      cor <- tryCatch(cor(prob[idx],test[idx,tcr[j]]),error = function(e) {return(NA)})
      mat[i,j] <- cor
    }
  }
  rownames(mat) <- hlas
  colnames(mat) <- tcr
  return(mat)
}
getPredictors <- function(y,pred,hla){
  # print("Yes")
  acc <- sum(y == (pred > 0.5),na.rm = T)/length(y)
  auc_all <- get_AUC(pred,y)
  r2 <- tjurs(y,pred)
  error <- errorRate(y,pred)#tryCatch(errorRate(y,pred))
  accuracies <- NULL
  r2s <- NULL
  rocs <- NULL
  aucs <- NULL
  fprs <- NULL
  fnrs <- NULL
  n <- NULL
  for (i in 1:ncol(hla)) {
    rows <- which(hla[,i] == 1)
    # print(rows)
    # print(data.frame(y = y[rows], pred = pred[rows]))
    n <- c(n,length(rows))
    trues <- y[rows] == (pred[rows] > 0.5)
    trues[is.na(trues)] <- FALSE
    #Calc accuracy for hla only
    accs <- sum(trues)/length(rows)
    if(is.nan(accs)) { accs <- 0 }
    accuracies <- c(accuracies,accs)
    #Calc ROC and AUC for hla only
    auc <- get_AUC(y[rows],pred[rows])
    aucs <- c(aucs,auc)
    #Calc R2 for hla only
    r2s <- c(r2s,tjurs(y[rows],pred[rows]))
    #Calc error rates
    errors <- errorRate(y[rows],pred[rows])#tryCatch(errorRate(y[rows],pred[rows]),
                                 #error = {return(list(sensitivity = NA,specificity = NA,FPR = NA,FNR = NA))})
    fprs <- c(fprs,errors$FPR)
    fnrs <- c(fnrs,errors$FNR)
  }
  df <- data.frame(Allele = c('All',colnames(hla)),n = c(length(y),n),
                   Acc = c(acc,accuracies),R2 = c(r2,r2s),AUC = c(auc_all,aucs),
                   FPR = c(error$FPR,fprs),FNR = c(error$FNR,fnrs))
  # print(df)
  return(df)
}

log.regress <- function(train,test,hla){
  #Train on train
  hla <- hla[rownames(test),]
  lim <- cv.glmnet(as.matrix(train[,2:ncol(train)]),train[,1],family = "binomial",alpha = 1)
  tpred <- predict(lim,as.matrix(train[,2:ncol(train)]),s = "lambda.1se",type = "response")
  pred <- predict(lim,as.matrix(test[,2:ncol(test)]),s = "lambda.1se",type = "response")
  if(length(grep("HLA",colnames(test))) == 0 | length(grep("TR",colnames(test))) == 0) {
   cor = NA 
  } else {
    cor <- calcHLATCR(pred,test)
  }
  lim_res <- list(Coefs = getCoefs(lim),Results = getPredictors(test[,1],pred,hla),Corr = cor,trainpred = tpred[,1], trainresponse = train$response) 
  return(lim_res)
}
mv.regress <- function(train,test,hla){
  train_hla <- hla[rownames(train),]
  hla <- hla[rownames(test),]
  k <- t(t(train_hla) * train[,1])
  hla_index <- grep("HLA",colnames(train))
  
  mv <- multivariateGlm.fit(k,train[,-c(1,hla_index)],family = rep("binomial",ncol(k)),size = NULL)
  coefs <- sapply(mv,coef)
  colnames(coefs) <- colnames(k)
  mvpred <- predict(mv,test[,-c(1,hla_index)],type = "response")
  mvlist <- lapply(mvpred,function(x) getPredictors(test[,1],x,hla))
  allacc <- do.call(cbind,lapply(mvlist,function(x) x$Acc))
  allr2 <- do.call(cbind,lapply(mvlist,function(x) x$R2))
  allauc <- do.call(cbind,lapply(mvlist,function(x) x$AUC))
  allfpr <- do.call(cbind,lapply(mvlist,function(x) x$FPR))
  allfnr <- do.call(cbind,lapply(mvlist,function(x) x$FNR))
  
  mv_res <- list(Coefs = coefs,
                 Inter = data.frame(Term = colnames(coefs),Intercept = coefs[1,]),
                 #as.data.frame(t(data.frame(apply(coefs,1,function(x) t.test(x)$conf.int),row.names = c("CImin","CImax"))))),
                 Results = cbind(mvlist[[1]][,1:2],
                                 data.frame(Acc = apply(allacc,1,mean),
                                            R2 = apply(allr2,1,mean), 
                                            AUC = apply(allauc,1,function(x) mean(x,na.rm=T)),
                                            FPR = apply(allfpr,1,function(x) mean(x,na.rm=T)),
                                            FNR = apply(allfnr,1,function(x) mean(x,na.rm=T)))))
  return(mv_res)
}
make_long <- function(df){
  df <- within(df,sample <- rownames(df))  %>%
    pivot_longer(cols = colnames(df)[grep("HLA",colnames(df))],names_to = "HLA",values_to = "Contains")
  df <- df[df$Contains >= 1,]
  df <- df[,!colnames(df) == "Contains"]
  tmp <- as.data.frame(df)
  colnames(tmp) <- gsub("\\.","_",colnames(tmp))
  tmp[,1] <- as.factor(tmp[,1])
  tcrindex <- grep("TR",colnames(df))
  tmp[,tcrindex] <- apply(tmp[,tcrindex],1,function(x) as.numeric(x)*1)
  return(tmp)
}
lmm.regress <- function(train,test,hla){
  train.df <- make_long(train)
  test.df <- make_long(test)
  matchhla <- intersect(unique(train.df$HLA),unique(test.df$HLA))
  train.df <- train.df[train.df$HLA %in% matchhla,]
  test.df <- test.df[test.df$HLA %in% matchhla,]
  
  form <- as.formula(paste("response ~", paste(colnames(train.df)[4:ncol(train.df)-2], collapse = " + ")," + (1 | HLA) "))#+ (1 | Age)" )) #Make this on the test set?
  lmm <- glmer(form,data = train.df,
               family = binomial,control = glmerControl(optimizer = "nloptwrap",calc.derivs = FALSE,optCtrl=list(maxfun=2e5)),nAGQ = 0)
  coefs <- coef(lmm)
  coef_table <- summary(lmm)$coefficients
  #Predict
  train_pred <- lmm_pred <- predict(lmm,train.df,type = "response")
  train_pred <- split(train_pred,train.df$sample)#Get the median prediction 
  train_pred <- unlist(lapply(train_pred,function(x) mean(as.numeric(x))))
  
  lmm_pred <- predict(lmm,test.df,type = "response")
  lmm_pred <- split(lmm_pred,test.df$sample)#Get the median prediction 
  lmm_pred <- unlist(lapply(lmm_pred,function(x) mean(as.numeric(x))))
  cor <- calcHLATCR(lmm_pred,test)
  lmm_res <- list(Coef = data.frame(Term = rownames(coef_table),Coef = coef_table[,1], Error = coef_table[,2],pvalue = coef_table[,4]),
                  Inter = data.frame(Term = rownames(coefs$HLA),Intercept = coefs$HLA[,1]),
                  Results = getPredictors(test[names(lmm_pred),]$response,lmm_pred,hla[names(lmm_pred),]),
                  Corr = cor,
                  trainpred = train_pred,
                  trainresponse = train[names(train_pred),]$response)
  return(lmm_res)
}

maketraintest <- function(df){
  features <- df
  set.seed(12345)
  train_index <- createDataPartition(1:nrow(features), p=0.70, list=FALSE)
  train <- features[train_index,]
  validation <- features[-train_index,]
  # age_index <- 2
  hla_index <- grep("HLA",colnames(train))
  tcr_index <- grep("TR",colnames(train))
  return(list(train = train, test = validation, hla = hla_index,tcr = tcr_index))
}

prepFeatures <- function(status,pheno,tcr,alleles,ids){
  ids <- ids[ids %in% colnames(tcr)]
  sub <- pheno[ids,]
  status <- status[ids]
  response <- ifelse(factor(status) == "case", 1, 0)
  names(response) <- names(status)
  df <- cbind(as.data.frame(response),
              sub,
              data.frame(alleles[ids,],
                         t(tcr)[ids,]))
  rownames(df) <- ids
  print(dim(df))
  keep <- apply(tcr[,ids],2,function(x) any(is.na(x))) | apply(alleles[ids,],1,function(x) any(is.na(x))) | apply(pheno[ids,],1,function(x) any(is.na(x)))
  df <- df[!keep,]
  print(dim(df))
  return(df)
} #Function to prep get the data into a format suitable for statistical modeling

plotModels <- function(ls){
  #Get separate dataframes for Acc, R2, AUC
  res <- rbindlist(lapply(ls,function(x) x$Results[x$Results$n > 0,]),idcol = "Model")
  print("yes")
  reslong <- res %>% pivot_longer(cols = c("Acc","R2","AUC","FPR","FNR"),names_to = "Measurement",values_to = "Value")
  
  # g0 <- ggplot(reslong[reslong$Allele == "All",],aes(y = Value,x = Model)) +
  g0 <- ggplot(res[res$Allele == "All",],aes(y = Acc,x = Model)) +
    geom_bar(stat = "identity") +
    ylab("Accuracy %") +
    xlab("Model") +
    scale_y_continuous(limits = c(0,1)) + 
    pubrot
  
  cols <- rep("black",nrow(res))
  cols[res$Allele == "All"] <- "green"
  cols[res$Allele %in% prott] <- "blue"
  cols[res$Allele %in% riskt] <- "red"
  res <- within(res,col <- cols)
  
  
  g1 <- ggplot(res,aes(y = Acc, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("Accuracy (%)") +
    xlab ("Model") +
    pubrot
  g2 <- ggplot(res,aes(y = R2, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Correlation") +
    ylab("Tjur's R^2") +
    xlab ("Model") +
    pubrot
  g3 <- ggplot(res,aes(y = AUC, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("AUC") +
    xlab ("Model") +
    pubrot
  g4 <- ggplot(res,aes(y = FPR, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("FPR (%)") +
    xlab ("Model") +
    pubrot
  g5 <- ggplot(res,aes(y = FNR, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("FNR (%)") +
    xlab ("Model") +
    pubrot
  #Also do for coefs
  coefs <- rbindlist(coefs <- list(GLM = ls$GLM$Coefs,GLM_TCR = ls$GLM_TCR$Coefs,
                                   # MV = (within(as.data.frame(ls$MV$Coefs),Term <- rownames(ls$MV$Coefs)) %>% pivot_longer(cols = colnames(ls$MV$Coefs),
                                   #                                                                                                        names_to = "HLA",values_to = "Coef"))[,c(1,3)],
                                   LMM = ls$LMM$Coef[,1:2]),
                     idcol = "Model",fill = T)
  coefs$Term <- gsub("_","\\.",coefs$Term)
  coefs <- coefs[grep("TR",coefs$Term),]   #!grep("HLA",coefs$Term),]
  # coefs <- coefs[!coefs$Term %in% c("(Intercept)","Age"),]
  terms <- unique(unlist(lapply(split(coefs,coefs$Model),function(x){
    term <- x[rev(order(x$Coef)),]$Term
    return(c(term[1:10],tail(term,n = 10)))
  })))
  
  terms <- terms[order(coefs[coefs$Term %in% terms & coefs$Model == "LMM",]$Coef)]
  coefs$Term <- factor(coefs$Term,levels = terms)
  g6 <- ggplot(coefs[coefs$Term %in% terms,],aes(x = Coef,y = Term)) +
    # geom_point(size = 2.5) +
    geom_bar(stat="identity") +
    ggtitle("Coefficient") +
    ylab("TCR gene") +
    # geom_vline(xintercept = 0,linetype = "dashed",color = "black")+
    facet_wrap(~Model,nrow = 1,scales = "free_x") +
    ylab ("Importance") +
    pub
  #Also for intercepts
  inters <- rbindlist(list(GLM = setNames(ls$GLM$Coefs[grep("HLA",ls$GLM$Coefs$Term),],c("Term","Intercept")),
                           GLM_HLA = setNames(ls$GLM_HLA$Coefs[grep("HLA",ls$GLM_HLA$Coefs$Term),],c("Term","Intercept")),
                           # MV = Full_list$MV$Inter,
                           LMM = ls$LM$Inter),idcol = "Model")
  cols <- rep("black",nrow(inters))
  # cols[inters$Term == "All"] <- "green"
  cols[inters$Term %in% prott] <- "blue1"
  cols[inters$Term %in% riskt] <- "red1"
  inters <- within(inters,col <- cols)
  terms <- unique(unlist(lapply(split(inters,inters$Model),function(x){
    term <- x[rev(order(x$Intercept)),]$Term
    return(c(term[1:10],tail(term,n = 10)))
  })))
  g7 <- ggplot(inters[inters$Term %in% terms,],aes(x = Intercept,y = Term)) +
    # geom_point(colour = inters[inters$Term %in% terms,]$col,size = 2.5) +
    geom_bar(stat="identity",fill = inters[inters$Term %in% terms,]$col)+
    ggtitle("Intercept") +
    ylab("HLA allele") +
    # geom_vline(xintercept = 0,linetype = "dashed",color = "black")+
    facet_wrap(~Model,nrow = 1,scales = "free_x") +
    pub
  return(list(g0,g1,g2,g3,g4,g5,g6,g7))
}

runStatTests <- function(train,test,ls,h){
  hla_index <- ls$hla
  tcr_index <- ls$tcr
  model_list = list(NULL)
  model_list[["GLM"]] <- log.regress(train,test,h)
  print(1)
  model_list[['GLM_HLA']] <- log.regress(train[,-tcr_index],test[,-tcr_index],h)
  print(2)
  model_list[['GLM_TCR']] <- log.regress(train[,-hla_index],test[,-hla_index],h)
  print(3)
  # model_list[['MV']] <- tryCatch(mv.regress(train,test,h),error = function(e) {print("MV Fail"),return(NA)})
  print(4)
  model_list[['LMM']] <- lmm.regress(train,test,h)
  return(model_list)
  # return(list(
  #   GLM = tryCatch(log.regress(train,test,h),error = { print("GLM Fail")}),#; return(NA)}),
  #   GLM_HLA = tryCatch(log.regress(train[,-tcr_index],test[,-tcr_index],h),error = { print("GLM HLA Fail")}),#; return(NA)}),
  #   GLM_TCR = tryCatch(log.regress(train[,-hla_index],test[,-hla_index],h),error = { print("GLM TCR Fail")}),#; return(NA)}),
  #   MV = tryCatch(mv.regress(train,test,h),error = { print("MV Fail")}),#; return(NA)}),
  #   LMM = tryCatch(lmm.regress(train,test,h),error = { print("LMM Fail")}),#}; return(NA)})
  # ))
}
##############################################
######### Run with train/test split###########
##############################################

########################Get the samples ready
subsample <- sample[colnames(TCRnorm),]
term_ids <- rownames(subsample)[subsample$Terminal] #ID's of samples that are the oldest for a patient (terminal)
non_term_ids <- colnames(TCRnorm)[!colnames(TCRnorm) %in% term_ids] #ID's of samples that are not the oldest (presumably before diabetic onset)
subterm <- subsample[rownames(subsample) %in% rownames(term_meta),]
cont_ids <- rownames(subterm[subterm$Status == "control",])
ca_ids <- rownames(subterm[subterm$Status == "case",])
equal_ids <- c(cont_ids,
               sample(ca_ids,length(cont_ids)))
#######################Prep the TCR features
scaleTCR <- apply(TCRnorm,2,scale)
rownames(scaleTCR) <- rownames(TCRnorm)
boxplot(TCRnorm) #Check that scaleTCR is correctly normalized (remove outliers)
plot(apply(scaleTCR,2,median))
abline(h = c(-0.24,-0.36), col = "red", lty = 2)
outliers <- which(apply(scaleTCR,2,median) > -0.24 | apply(scaleTCR,2,median) < -0.37)
norm_scaleTCR <- scaleTCR[,-outliers]
plot(apply(norm_scaleTCR,2,median))



########################Get the features into a format for modeling
response <- subsample$Status
names(response) <- rownames(subsample)
sampallele <- alleles[subsample$subject_id,]
rownames(sampallele) <- rownames(subsample)

d.df <- prepFeatures(response,subsample[,c(3,4)],scaleTCR,sampallele,rownames(subsample))
hla <- d.df[,grep("HLA",colnames(d.df))]

term.df <- prepFeatures(response,subsample[,c(3,4)],scaleTCR,sampallele,term_ids)
nterm.df <- prepFeatures(response,subsample[,c(3,4)],scaleTCR,sampallele,non_term_ids)

equal.df <- prepFeatures(response,subsample[,c(3,4)],scaleTCR,sampallele,equal_ids)

#Check for correlation between the 
# corrplot(cor(d.df/colSums(d.df)))

#there is obvious bias in the allele model (probably because there are )
#########################Separate training and testing datasets (70 train/30 test)

Full <- maketraintest(d.df)
Term <- maketraintest(term.df)
nTerm <- maketraintest(nterm.df)
Equal <- maketraintest(equal.df)

Fullnoage <- Full
Fullnoage[["train"]] <- Fullnoage[["train"]][,-3]
Fullnoage[["test"]] <- Fullnoage[["test"]][,-3]

Equalnoage <- Equal
Equalnoage[["train"]] <- Equalnoage[["train"]][,-3]
Equalnoage[["test"]] <- Equalnoage[["test"]][,-3]
#Make a function for each model type that will take in the train and predict on the test


# log.regress(rnaFull$train,rnaFull$test,rnahla)
# mv.regress(Equal[["train"]],Equal[["test"]],hla)
# lmm.regress(Full[["train"]],Full[["test"]],hla)

##########################################################################
#####make lists containing the different models we want to develop########
##########################################################################

Full_list <- runStatTests(Full$train,Full$test,Full,hla) 
Full_list_noage <- runStatTests(Fullnoage$train,Fullnoage$test,Fullnoage,hla)
Term_list <- runStatTests(Term$train,Term$test,Term,hla) 
nTerm_list <- runStatTests(nTerm$train,nTerm$test,nTerm,hla) 

Term_nTerm_list <- runStatTests(rbind(Term$train,Term$test),rbind(nTerm$train,nTerm$test),Term,hla) 
nTerm_Term_list <- runStatTests(rbind(nTerm$train,nTerm$test),rbind(Term$train,Term$test),nTerm,hla) 

Equal_list <- runStatTests(Equal$train,Equal$test,Equal,hla)
Equal_noage <- runStatTests(Equalnoage$train,Equalnoage$test,Equal,hla)

##############################################
#### Compare models ##########################
##############################################

#Lets write a function that will create all the plots within each of the lists. 
#We want to create a plot with each models Accuracy, R2, AUC and Coef (ranked top 20)

# tmp <- rbindlist(lapply(Term_nTerm_list,function(x) x$Results[x$Results$n > 0,]),idcol="Model")
# ggplot(tmp,aes(y = AUC, x = Model)) +
#   geom_boxplot()+
#   geom_jitter()+
#   # geom_point() +
#   # facet_wrap(~Model)+
#   ggtitle("Accuracy")

Full_plots <- plotModels(Full_list)
Full_noage_plots <- plotModels(Full_list_noage)
Term_plots <- plotModels(Term_list)
nTerm_plots <- plotModels(nTerm_list)
T_nt_plots <- plotModels(Term_nTerm_list)
T_t_plots <- plotModels(nTerm_Term_list)
equal_plots <- plotModels(Equal_list)
equal_noage_plots <- plotModels(Equal_noage)
#Now make the plots out
#First all the accuracies side by side
Datasets <- c("Full","Full no age","Terminal","nonTerminal","TermxnTerm","nTermxTerm","Equal")
x <- 6
ggarrange(plotlist = list(Full_plots[[x]],
                          Full_noage_plots[[x]],
                          Term_plots[[x]],
                          nTerm_plots[[x]],
                          T_nt_plots[[x]],
                          T_t_plots[[x]],
                          equal_plots[[x]]),
          labels = Datasets,nrow = 1,label.x = 0.6,label.y = 0.2)


# equal_noage_plots[[1]]
#Plot the accuracies comparison
ggarrange(plotlist = list(Full_plots[[1]] + ggtitle("train 70/test 30"),
                          T_t_plots[[1]] + ggtitle("train non-terminal/test terminal") + ylab("")))
ggsave(".//Results/IRepertoire/TCRfigs/Fig3Acc_comp.svg",width = 10,height = 5,units = "in",dpi = 600 )
#Plot the coefficients for the full model

Full_plots[[7]]
ggsave(".//Results/IRepertoire/TCRfigs/Fig3Stat_coef_all.svg",width = 7,height = 10,units = "in",dpi = 600 )

#Plot inter
Full_plots[[8]]
ggsave(".//Results/IRepertoire/TCRfigs/Fig3Stat_hla_inter.svg",width = 7,height = 10,units = "in",dpi = 600 )


#Plot the typeIandtype2 error
tmp.df <- rbindlist(lapply(Term_nTerm_list,function(x) x$Results[x$Results$Allele == "All",]),idcol="Model")
g1 <- ggplot(tmp.df[tmp.df$Model != "MV",],aes(x = Model,y = FPR)) +
  geom_bar(stat = "identity") +
  ggtitle("Type 1 Error") +
  ylab("Rate (%)") +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  theme(axis.text.x = element_blank() )+
  pubrot
g2 <- ggplot(tmp.df[tmp.df$Model != "MV",],aes(x = Model,y = FNR)) +
  geom_bar(stat = "identity") +
  ggtitle("Type II Error") +
  ylab("Rate (%)") +
  xlab("Model") +
  pubrot
ggarrange(plotlist = list(g1 = g1,g2 = g2),nrow = 2,heights = c(0.7,0.9))

#############Overfit
trainacc <- rbindlist(lapply(Full_list,function(x) data.frame(Accuracy = sum((x$trainpred > 0.5) == x$trainresponse)/length(x$trainpred))),idcol = "Model")
a <- ggplot(trainacc[!is.nan(trainacc$Accuracy),],aes(x = Model, y = Accuracy)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1))+
  pubrot
b <- Full_plots[[1]]
ggarrange(plotlist = list(a = a,b = b),nrow = 1,labels = c("Training Accuracy","Testing Accuracy"))

#Compile the top features from each model
get_top <- function(x){
  t <- c(as.character(x$GLM$Coefs[c(1:10,(nrow(x$GLM$Coefs)-10):nrow(x$GLM$Coefs)),]$Term))
  return(t)
}



##############################################
#### Machine Learning ########################
##############################################

#The problem is that some patients have 0 recorded HLA so we remove those patient files first
features <- data.frame(response = subsample$Status,
                       Age = subsample$age_max,
                       alleles[subsample$subject_id,],
                       t(scaleTCR)[rownames(subsample),],row.names = rownames(subsample))
# features <- d.df
# colnames(features) <- gsub("\\.","_",colnames(features))
features <- features[,!apply(features,2,function(x) any(is.na(x)))] #get rid of any NA in features
# features <- features[,!duplicated(t(features))]
features <- features[,-which(apply(features,2,var) == 0)]
#Now use machine learning code to develop multiple models predicting the response based on TCR and allele usage
#Make a train test set
set.seed(12345)
train_index <- createDataPartition(1:nrow(features), p=0.70, list=FALSE)
train <- features[train_index,]
validation <- features[-train_index,]

trainModels <- function(train,validation){
  train_age <- train$Age
  val_age <- validation$Age
  train <- train[,!colnames(train) == "Age"]
  validation <- validation[,!colnames(validation) == "Age"]
  
  #Iterate through models
  
  control <- trainControl(method="repeatedcv", number=10, repeats=3,classProbs = T)
  
  m_bag <- train(as.factor(response) ~ ., data=train, method="cforest", metric="ROC", 
                 trControl=control, preProcess = c("center", "scale") )
  m_sda <- train(as.factor(response)~., data=train, method="sda", metric="ROC", 
                 trControl=control, preProcess = c("center", "scale") )
  m_mlp <- train(as.factor(response)~., data=train, method="mlp", metric="ROC", 
                 trControl=control, preProcess = c("center", "scale") )
  m_glm <- train(as.factor(response)~., data=train, method="glm", metric="ROC", 
                 trControl=control, preProcess = c("center", "scale") )
  m_boost <- train(as.factor(response)~., data=train, method="gbm", metric="ROC", 
                   trControl=control, preProcess = c("center", "scale") )
  m_svm <- train(as.factor(response)~., data=train, method="svmLinear2", metric="ROC", 
                 trControl=control, preProcess = c("center", "scale") )
  
  # resample_results <- resamples(list(Bag=m_bag,SDA=m_sda, 
                                     # MLP=m_mlp, GLM=m_glm, Boost=m_boost, SVM=m_svm))
  model_list <- list(Bag=m_bag,SDA=m_sda, 
                 MLP=m_mlp, GLM=m_glm, Boost=m_boost, SVM=m_svm)
  # summary(resample_results,metric = c("Kappa","Accuracy"))
  # bwplot(resample_results , metric = c("Accuracy","ROC"))
  
  #Then predict on the validation set 
  y <- rep(0,nrow(validation))#as.character(validation$encoding)
  y[validation$response == "case"] <- 1
  calcAcc <- function(x,y){ sum((x > 0.5) == (y == 1))/length(y)}
  
  p1 <- predict(m_bag, validation, preProcess = c("center", "scale"),type = "prob")
  p2 <- predict(m_sda, validation, preProcess = c("center", "scale"),type = "prob")
  p3 <- predict(m_mlp, validation, preProcess = c("center", "scale"),type = "prob")
  p4 <- predict(m_glm, validation, preProcess = c("center", "scale"),type = "prob")
  p5 <- predict(m_boost, validation, preProcess = c("center", "scale"),type = "prob")
  p6 <- predict(m_svm, validation, preProcess = c("center", "scale"),type = "prob")
  models <- c("Bagging","SDA","MLP","GLM","Boost","SVM")
  accuracies <- c(calcAcc(p1$case,y),calcAcc(p2$case,y),calcAcc(p3$case,y),calcAcc(p4$case,y),calcAcc(p5$case,y),calcAcc(p6$case,y))
  pos_accuracies <- unlist(lapply(list(p1,p2,p3,p4,p5,p6),function(x) calcAcc(x = x$case[y == 1],y = y[y==1])))
  neg_accuracies <- unlist(lapply(list(p1,p2,p3,p4,p5,p6),function(x) calcAcc(x = x$case[y == 0],y = y[y==0])))
  
  plist <- list(bag = p1,sda = p2,mlp = p3,glm = p4,boost = p5,svm = p6)
  plist <- lapply(plist,function(x) within(cbind(validation[,1],val_age,x),acc <- (x$case > 0.05) == (validation$response == "case")))
  
  #ROC
  get_ROC <- function(x,y){
    proc_obj <- roc(predictor = x,response = y,ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                    print.auc=TRUE, show.thres=TRUE)
    net_roc <- data.frame(y = rev(proc_obj$sensitivities),x = rev(proc_obj$specificities))
    return(net_roc)
  }
  get_AUC <- function(x,y){
    proc_obj <- roc(predictor = x,response = y,ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                    print.auc=TRUE, show.thres=TRUE)
    net_auc <- auc(proc_obj)
    return(net_auc)
  }
  
  predlist <- list(p1,p2,p3,p4,p5,p6)
  names(predlist) <- models
  roclist <- lapply(predlist,function(x) list(roc = get_ROC(x$case,y), auc = get_AUC(x$case,y)))
  names(roclist) <- models
  roc.df <- rbindlist(lapply(roclist,function(x) x[[1]]),idcol = "Model")
  
  a <- ggplot(roc.df,aes(x = x,y = y, colour = Model)) +
    geom_line(linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 1,linetype = "dashed") +
    scale_x_reverse() +
    # scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    pub
  
  df <- data.frame(Model = models,Accuracy = accuracies,AUC = unlist(lapply(roclist,function(x) x[[2]])),Pos_acc = pos_accuracies, Neg_acc = neg_accuracies)
  return(list(models = model_list,AUC = a,data = df,plist = plist,ROC = roc.df))
}

#Now train for different feature levels
#May need to hide age
hlaonly <- trainModels(train[,-grep("TR",colnames(features))],validation[,-grep("TR",colnames(features))])
tcronly <- trainModels(train[,-grep("HLA",colnames(features))],validation[,-grep("HLA",colnames(features))])
both <- trainModels(train,validation)

#Plot
resample_results <- resamples(hlaonly[[1]])
hla_plot <- bwplot(resample_results , metric = c("Accuracy","ROC"),main = "HLA only",xlim = c(0.5,1))

resample_results <- resamples(tcronly[[1]])
tcr_plot <- bwplot(resample_results , metric = c("Accuracy","ROC"),main = "TCR only",xlim = c(0.5,1))


resample_results <- resamples(both[[1]])
both_plot <- bwplot(resample_results , metric = c("Accuracy","ROC"),main = "HLA + TCR",xlim = c(0.5,1))

par(mfrow = c(1, 3))
hla_plot
tcr_plot
both_plot
par(mfrow = c(1, 1))

#AUC
hlaonly$AUC + ggtitle("HLA only")https://ansc-ilamb.ansc.illinois.edu/graphics/plot_zoom_png?width=774&height=900
tcronly$AUC + ggtitle("TCR only")
both$AUC + ggtitle("TCR+HLA")

#Feature importance 

hlafeatures <- varImp(hlaonly$models$Bag,scale = F)
tcrfeatures <- varImp(tcronly$models$Bag,scale = F)
bothfeatures <- varImp(both$models$GLM,scale = F)

plot(hlafeatures,top=20,main = "HLA features")
plot(tcrfeatures,top=20,main = "TCR features")
plot(bothfeatures,top=20,main = "HLA + TCR features")


training <- both$models$SVM$trainingData
filterVarImp(x = both$models$SVM, y = training$.outcome)
varImp(both$models$SVM,scale = F,useModel = T)


hlaonly$models$Bag
hlaonly$models$SDA$finalModel$beta

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(features[,2:ncol(features)], as.factor(features[,1]), sizes=round(seq(from = 1, to = ncol(features)-1, length.out = 20)), rfeControl=control)
a <- ggplot(results) + theme_bw() + pub
best <- results$variables[results$variables$Variables == 101,]
best <- best[!duplicated(best$var),]
best <- best[order(best$Overall,decreasing = T),]
b <- ggplot(best[1:20,],aes(x = Overall,y = factor(var,levels = var))) +
  geom_bar(stat = "identity") +
  ylab("") +
  pub
ggarrange(plotlist=list(a,b),ncol = 1,heights = c(0.4,0.6))
#############################################################
########bulk RNA statistics ###################################
#############################################################

rnameta <- read.csv(".//Data/RNA/metadata.tsv",header = T,sep = "\t")
rownames(rnameta) <- rnameta$Run
head(rnameta)
rnameta$Risk[rnameta$Risk == "Protective"] <- "Protected"
rnameta$Risk[rnameta$Risk == "Risk"] <- "At-risk"


rnameta <- within(rnameta,Terminal <- rnameta$Age_at_collection > 12)
table(rnameta$Status,rnameta$Terminal)


rnahla <- read.table(".//Data/RNA/HLA.tsv",header = T,sep = "\t")

# rnatcr <- read.table(".//Data/RNA/TCRexp.tsv",header = T,sep = "\t")
rnatcr <- read.table(".//Data/RNA/TCRnorm.tsv",header = T,sep = "\t")

boxplot(rnatcr)
plot(colSums(rnatcr)) #fairly even (I believe we used the sct transformation which develops a negative binomial)

#Diff exp 
calcDiff <- function(dat,case,control){
  case <- case[case %in% colnames(dat)]
  control <- control[control %in% colnames(dat)]
  case_arr <- dat[,case]
  control_arr <- dat[,control]
  df <- data.frame(NULL)
  for (i in 1:nrow(case_arr)){
    x <- removeNA(case_arr[i,])
    y <- removeNA(control_arr[i,])
    gene <- rownames(case_arr)[i]
    fc <- mean(abs(x))/mean(abs(y))
    lfc <- log2(fc)
    if (var(x) == 0 | var(y) == 0){ pval <- 1 } else {
      test <- t.test(x,y)
      pval <- test$p.value
    }
    row <- data.frame(Gene = gene, 
                      mean_case = mean(x), mean_control = mean(y),
                      sd_case = sd(x),sd_control = sd(y),
                      logfc = lfc,#log2(abs(mean(x))/abs(mean(y))),
                      p.value = pval,
                      padj = pval * nrow(case_arr))
    df <- rbind(df,row)
  }
  return(df)
} 
diff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case",]),rownames(rnameta[rnameta$Status == "control",]))
sigdiff <- diff[diff$padj < 0.05,]
termdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Age_at_collection > 12,]),rownames(rnameta[rnameta$Status == "control" & rnameta$Age_at_collection > 12,]))
ntermdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Age_at_collection < 12,]),rownames(rnameta[rnameta$Status == "control" & rnameta$Age_at_collection < 12,]))
bothdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Risk == "Both",]),rownames(rnameta[rnameta$Status == "control" & rnameta$Risk == "Both",]))
neitherdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Risk == "Neither",]),rownames(rnameta[rnameta$Status == "control" & rnameta$Risk == "Neither",]))
riskdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Risk == "At-risk",]),rownames(rnameta[rnameta$Status == "control" & rnameta$Risk == "At-risk",]))
protdiff <- calcDiff(rnatcr,rownames(rnameta[rnameta$Status == "case" & rnameta$Risk == "Protected",]),rownames(rnameta[rnameta$Status == "control" & rnameta$Risk == "Protected",]))
Bulk_res <- rbindlist(list(All = diff[diff$padj < 0.05,],
                           Terminal = termdiff[termdiff$padj < 0.05,],
                           Nonterminal = ntermdiff[ntermdiff$padj < 0.05,],
                           Both = bothdiff[bothdiff$padj < 0.05,],
                           Neither = neitherdiff[neitherdiff$padj < 0.05,],
                           Risk = riskdiff[riskdiff$padj < 0.05,],
                           Protective = protdiff[protdiff$padj < 0.05,]),idcol = "Comparison")

write.table(Bulk_res,file = ".//Results/RNAseq/Table/SuppTableDET.tsv",sep = "\t")

vl <- lapply(list(Both = bothdiff,Neither = neitherdiff,Protected = protdiff, At_risk = riskdiff),volcanoPlot)
ggarrange(plotlist = vl,ncol =1,labels = names(vl))

########set ids
subrna_meta <- rnameta[colnames(rnatcr),]
allrna_ids <-  rownames(subrna_meta)
termrna_ids <- rownames(subrna_meta[subrna_meta$Terminal == T & !is.na(subrna_meta$Terminal),])
ntermrna_ids <- rownames(subrna_meta[subrna_meta$Terminal == F & !is.na(subrna_meta$Terminal),])
equaltermids <- sample(termrna_ids,length(ntermrna_ids),replace = F)
#Now put together the features
response <- subrna_meta$Status
names(response) <- rownames(subrna_meta)

rna.df <- prepFeatures(response,subrna_meta[,c(6,7)],rnatcr,rnahla,allrna_ids)


rnaFull <- maketraintest(rna.df)
rnaTerm <- maketraintest(prepFeatures(response,subrna_meta[,c(6,7)],rnatcr,rnahla,termrna_ids))
rnanTerm <- maketraintest(prepFeatures(response,subrna_meta[,c(6,7)],rnatcr,rnahla,ntermrna_ids))
equalrnaTerm <- maketraintest(prepFeatures(response,subrna_meta[,c(6,7)],rnatcr,rnahla,equaltermids))

rnaFull_list <- runStatTests(rnaFull$train,rnaFull$test,rnaFull,rnahla) 
rnaTerm_list <- runStatTests(rnaTerm$train,rnaTerm$test,rnaTerm,rnahla) 
# rnanTerm_list <- runStatTests(rnanTerm$train,rnanTerm$test,rnanTerm,rnahla) 

rnaT_nT_list <- runStatTests(rbind(rnaTerm$train,rnaTerm$test),rbind(rnanTerm$train,rnanTerm$test),rnaTerm,rnahla)
# rnanT_T_list <- runStatTests(rbind(rnanTerm$train,rnanTerm$test),rbind(rnaTerm$train,rnaTerm$test),rnaTerm,rnahla)
equal_list <- runStatTests(rbind(equalrnaTerm$train,equalrnaTerm$test),rbind(rnanTerm$train,rnanTerm$test),equalrnaTerm,rnahla)

plotrnaModels <- function(ls){
  #Get separate dataframes for Acc, R2, AUC
  res <- rbindlist(lapply(ls,function(x) x$Results[x$Results$n > 0,]),idcol = "Model")
  print("yes")
  reslong <- res %>% pivot_longer(cols = c("Acc","R2","AUC","FPR","FNR"),names_to = "Measurement",values_to = "Value")
  
  # g0 <- ggplot(reslong[reslong$Allele == "All",],aes(y = Value,x = Model)) +
  g0 <- ggplot(res[res$Allele == "All",],aes(y = Acc,x = Model)) +
    geom_bar(stat = "identity") +
    ylab("Accuracy %") +
    xlab("Model") +
    scale_y_continuous(limits = c(0,1)) + 
    pubrot
  
  cols <- rep("black",nrow(res))
  matchal <- unlist(lapply(strsplit(x = res$Allele,"\\."),function(x) paste(x[-length(x)],collapse=".")))
  cols[res$Allele == "All"] <- "green"
  cols[matchal %in% prott] <- "blue"
  cols[matchal %in% riskt] <- "red"
  res <- within(res,col <- cols)
  
  
  g1 <- ggplot(res,aes(y = Acc, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("Accuracy (%)") +
    xlab ("Model") +
    pubrot
  g2 <- ggplot(res,aes(y = R2, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Correlation") +
    ylab("Tjur's R^2") +
    xlab ("Model") +
    pubrot
  g3 <- ggplot(res,aes(y = AUC, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("AUC") +
    xlab ("Model") +
    pubrot
  g4 <- ggplot(res,aes(y = FPR, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("FPR (%)") +
    xlab ("Model") +
    pubrot
  g5 <- ggplot(res,aes(y = FNR, x = Model)) +
    geom_boxplot() +
    geom_jitter(colour = res$col,alpha = 0.6) +
    # ggtitle("Accuracy") +
    ylab("FNR (%)") +
    xlab ("Model") +
    pubrot
  #Also do for coefs
  coefs <- rbindlist(coefs <- list(GLM = ls$GLM$Coefs,GLM_TCR = ls$GLM_TCR$Coefs,
                                   # MV = (within(as.data.frame(ls$MV$Coefs),Term <- rownames(ls$MV$Coefs)) %>% pivot_longer(cols = colnames(ls$MV$Coefs),
                                   #                                                                                                        names_to = "HLA",values_to = "Coef"))[,c(1,3)],
                                   LMM = ls$LMM$Coef[,1:2]),
                     idcol = "Model",fill = T)
  coefs$Term <- gsub("_","\\.",coefs$Term)
  coefs <- coefs[grep("TR",coefs$Term),]   #!grep("HLA",coefs$Term),]
  # coefs <- coefs[!coefs$Term %in% c("(Intercept)","Age"),]
  terms <- unique(unlist(lapply(split(coefs,coefs$Model),function(x){
    term <- x[rev(order(x$Coef)),]$Term
    return(c(term[1:10],tail(term,n = 10)))
  })))
  
  terms <- terms[order(coefs[coefs$Term %in% terms & coefs$Model == "LMM",]$Coef)]
  coefs$Term <- factor(coefs$Term,levels = terms)
  g6 <- ggplot(coefs[coefs$Term %in% terms,],aes(x = Coef,y = Term)) +
    # geom_point(size = 2.5) +
    geom_bar(stat="identity") +
    ggtitle("Coefficient") +
    ylab("TCR gene") +
    # geom_vline(xintercept = 0,linetype = "dashed",color = "black")+
    facet_wrap(~Model,nrow = 1,scales = "free_x") +
    ylab ("Importance") +
    pub
  #Also for intercepts
  inters <- rbindlist(list(GLM = setNames(ls$GLM$Coefs[grep("HLA",ls$GLM$Coefs$Term),],c("Term","Intercept")),
                           GLM_HLA = setNames(ls$GLM_HLA$Coefs[grep("HLA",ls$GLM_HLA$Coefs$Term),],c("Term","Intercept")),
                           # MV = Full_list$MV$Inter,
                           LMM = ls$LM$Inter),idcol = "Model")
  cols <- rep("black",nrow(inters))
  matchal <- unlist(lapply(strsplit(x = inters$Term,"\\."),function(x) paste(x[-length(x)],collapse=".")))
  cols[matchal %in% prott] <- "blue1"
  cols[matchal %in% riskt] <- "red1"
  inters <- within(inters,col <- cols)
  terms <- unique(unlist(lapply(split(inters,inters$Model),function(x){
    term <- x[rev(order(x$Intercept)),]$Term
    return(c(term[1:10],tail(term,n = 10)))
  })))
  g7 <- ggplot(inters[inters$Term %in% terms,],aes(x = Intercept,y = Term)) +
    # geom_point(colour = inters[inters$Term %in% terms,]$col,size = 2.5) +
    geom_bar(stat="identity",fill = inters[inters$Term %in% terms,]$col)+
    ggtitle("Intercept") +
    ylab("HLA allele") +
    # geom_vline(xintercept = 0,linetype = "dashed",color = "black")+
    facet_wrap(~Model,nrow = 1,scales = "free_x") +
    pub
  return(list(g0,g1,g2,g3,g4,g5,g6,g7))
}


Full_plots <- plotModels(rnaFull_list)
Term_plots <- plotModels(rnaTerm_list)
rnaT_nT_plots <- plotModels(rnaT_nT_list)\
equal_plots <- plotModels(equal_list)

Full_plots[[8]]

Datasets <- c("Full","Terminal","TermxnTerm","Equal")
x <- 1
ggarrange(plotlist = list(Full_plots[[x]],
                          Term_plots[[x]],
#                           nTerm_plots[[x]],
                          rnaT_nT_plots[[x]],
#                           T_t_plots[[x]],
                          equal_plots[[x]]
                      ),
          labels = Datasets,nrow = 1,label.x = 0.6,label.y = 0.2)


tmp.df <- rbindlist(lapply(rnaFull_list,function(x) x$Results[x$Results$Allele == "All",]),idcol="Model")
g1 <- ggplot(tmp.df,aes(x = Model,y = FPR)) +
  geom_bar(stat = "identity") +
  ggtitle("Type 1 Error") +
  ylab("Rate (%)") +
  xlab("") +
  scale_x_discrete(labels = NULL) +
  theme(axis.text.x = element_blank() )+
  pubrot
g2 <- ggplot(tmp.df,aes(x = Model,y = FNR)) +
  geom_bar(stat = "identity") +
  ggtitle("Type II Error") +
  ylab("Rate (%)") +
  xlab("Model") +
  pubrot
ggarrange(plotlist = list(g1 = g1,g2 = g2),nrow = 2,heights = c(0.7,0.9))


#####################We want a function that will calc the model prob correlation to the TCR number for 

mat <- rnaFull_list$GLM$Corr
dim(mat)
mat <- mat[apply(mat,1,function(x) any(!is.na(x))),]
mat <- mat[,apply(mat,2,function(x) any(!is.na(x)))]
dim(mat)

mat <- mat[order(rownames(mat)),]
#Cluster columns (TCR genes) and also make annotation
rowclust <- cutree(fastcluster::hclust(dist(mat)),k = 12)
colclust <- cutree(fastcluster::hclust(dist(t(mat))),k = 6)
rowanno <- data.frame(Group = str_replace_all(substr(rownames(mat), 1, 9), "[^[:alpha:]]", ""),Cluster = as.character(rowclust[rownames(mat)]))
rownames(rowanno) <- rownames(mat)
colanno <- data.frame(Cluster = as.character(colclust))
rownames(colanno) <- colnames(mat)
pheatmap(mat,annotation_row = rowanno,annotation_col = colanno,cluster_rows = T,fontsize_row = 5)

names(colclust[colclust == 3])
names(rowclust[!rowclust %in% c(1,4)])









############################################################
######Comparison to known TCR ##############################
############################################################

knowndf <- read.csv(".//Data/KnownTCR.txt",header = T,sep = "\t")
head(knowndf)

knowndf <- within(knowndf,gene <- unlist(lapply(strsplit(knowndf$Target,":"),function(x) x[[1]])))
genespec <- rbindlist(lapply(split(knowndf,knowndf$gene),function(x) data.frame(TCRgene = unique(c(x[,5],x[,6],x[,8],x[,9])))),idcol = "Gene")
genespec <- genespec[genespec$TCRgene != "",]
pheatmap(table(genespec$Gene,genespec$TCRgene))


knowndf <- knowndf[knowndf$TRBV != "TRBV7\x969",]

#analyze average expression of 
TCRsig <- unique(TCRdiff_res$Gene)
TCRstatfeatures <- unique(c(get_top(Full_list),get_top(Term_list),get_top(nTerm_list),get_top(Term_nTerm_list),get_top(nTerm_Term_list)))
TCRMLfeat <- rownames(bothfeatures$importance)[rev(order(abs(bothfeatures$importance$Overall)))][1:20]
RNAseg <- unique(Bulk_res$Gene)
RNAstatfeatures <- unique(c(get_top(rnaFull_list),get_top(rnaTerm_list),get_top(rnaT_nT_list)))
RNAcorclust <- names(colclust[colclust == 3])
scTCR <- c("TRBV21-1","TRBV28","TRAV13-1","TRAV12-1","TRAV26-1","TRBV5-1","TRBV12-4",
           "TRAV19","TRBV7-3","TRAV12-1","TRAV17","TRBV28","TRAV8-2","TRAV4","TRBV3-1",
           "TRAV13-1","TRAV36DV7","TRAV13-2","TRBV18","TRAV14DV4","TRAV2","TRAV26-1","TRAV25","TRAV5")

matchtcr <- function(x,y){
  nx <- unlist(lapply(strsplit(x, split = "[[:punct:]]"),function(z) z[1]))#gsub("[^[:alnum:] ]", "", x)
  ny <- unlist(lapply(strsplit(y, split = "[[:punct:]]"),function(z) z[1]))#gsub("[^[:alnum:] ]", "", y)
  px <- unlist(lapply(strsplit(x, split = "[[:punct:]]"),function(z) z[2]))#gsub("[^[:alnum:] ]", "", x)
  py <- unlist(lapply(strsplit(y, split = "[[:punct:]]"),function(z) z[2]))#gsub("[^[:alnum:] ]", "", y)
  # print(data.frame(nx,as.numeric(px)))
  # print(data.frame(ny,as.numeric(py)))
  # nx <- nx[nx != ""]
  # ny <- ny[ny != ""]
  # print(paste(nx,px,sep="_"))
  # print(paste(ny,py,sep="_"))
  matches <- lapply(nx,function(z) grep(z,ny,value=F,fixed = T))
  # print(matches)
  m <- NULL
  for(t in 1:length(matches)){
    sec <- px[t]
    check <- py[ matches[[t]] ]
    if(any(sec %in% check)){
      m <- c(m,x[t])
    }
  }
                    # matches <- (paste(nx,px,sep="_") %in% paste(ny,py)) #unlist(lapply(nx,function(z) agrepl(z, ny,costs = list(substitutions = 5, deletions = 2,insertions = 5), max.distance = 0.1)))
  # print(matches)
  # m <- x[matches]
  return(m[!is.na(m)])
}
matchtcr(tcr,TCRsig)


impdf <- data.frame(NULL)
percdf <- data.frame(NULL)
for(g in unique(knowndf$gene)){
  sub <- knowndf[knowndf$gene == g,]
  names <- paste(sub$Name,collapse=" & ")
  tcr <- unique(c(sub$TRAV,sub$TRBJ,sub$TRBV,sub$TRAJ))
  impdf <- rbind(impdf,data.frame(Gene = g,
             Clones = names,
             TCRDiff = paste(matchtcr(tcr,TCRsig),collapse=" & "),
             TCRStat = paste(matchtcr(tcr,TCRstatfeatures),collapse=" & "),
             TCRML = paste(matchtcr(tcr,TCRMLfeat),collapse=" & "),
             RNADiff = paste(matchtcr(tcr,RNAseg),collapse=" & "),
             RNAstat = paste(matchtcr(tcr,RNAstatfeatures),collapse=" & "),
             RNACor = paste(matchtcr(tcr,RNAcorclust),collapse=" & "),
             scTCR = paste(matchtcr(tcr,scTCR),collapse=" & ")))
  percdf <- rbind(percdf,data.frame(Gene = g,
                                   TCRDiff = length(matchtcr(tcr,TCRsig))/length(tcr),
                                   TCRStat = length(matchtcr(tcr,TCRstatfeatures))/length(tcr),
                                   TCRML = length(matchtcr(tcr,TCRMLfeat))/length(tcr),
                                   RNADiff = length(matchtcr(tcr,RNAseg))/length(tcr),
                                   RNAstat = length(matchtcr(tcr,RNAstatfeatures))/length(tcr),
                                   RNACor = length(matchtcr(tcr,RNAcorclust))/length(tcr),
                                   scTCR = length(matchtcr(tcr,scTCR))/length(tcr)))
}

write.table(impdf,file = ".//Results/importance.tsv",sep = "\t")



alltcr <- grep("TR",colnames(d.df),value = T)
alldf <- data.frame(TCR = alltcr,
          TCRDiff = alltcr %in% matchtcr(alltcr,TCRsig),
          TCRStat = alltcr %in% matchtcr(alltcr,TCRstatfeatures),
          TCRML = alltcr %in% matchtcr(alltcr,TCRMLfeat),
          RNADiff = alltcr %in% matchtcr(alltcr,RNAseg),
          RNAstat = alltcr %in% matchtcr(alltcr,RNAstatfeatures),
          RNACor = alltcr %in% matchtcr(alltcr,RNAcorclust),
          scTCR = alltcr %in% matchtcr(alltcr,scTCR))
rownames(alldf) <- alldf$TCR
allmat <- apply(as.matrix(alldf[,-1]),2,as.numeric)
rownames(allmat) <- rownames(alldf)
pheatmap(t(allmat[rowSums(allmat) > 0,]),cluster_rows = F)
######################################################
#################dPCR  results######################
######################################################


dpcr <- read.table(file = ".//Results/dPCR/dPCR_results2.csv",header = T,sep = ",")
dpcr


dpcr_long <- dpcr %>% pivot_longer(cols = c("TRBV14","TRBV29","TRAV16","TRAV7"))

ggplot(dpcr_long,aes(x = Cells * Mult,y = value,colour = name)) +
  geom_point()

ggplot(dpcr_long,aes(x = RNA,y = value,colour = name)) +
  geom_point()


ggplot(dpcr_long,aes(x = Time, y = value,color = Group, group = Group)) +
  stat_summary(fun = "mean",geom = "line") +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width = 0.5) +
  facet_wrap(~name)


ggplot(dpcr_long,aes(x = Time, y = value/B2M,color = Group, group = Group)) +
  stat_summary(fun = "mean",geom = "point",size = 3) +
  stat_summary(fun = "mean",geom = "line",size = 1.5) +
  stat_summary(fun.data = "mean_se",geom = "errorbar",width = 1, size = 1.2) +
  facet_wrap(~name,nrow = 1,scales = "free_y") +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test",size = 5) +
  ylab("Normalized concentration") +
  xlab("Age (weeks)") +
  geom_vline(xintercept = 12,linetype = "dashed")+
  guides(labels = NULL) +
  pub



##########################################################
######## Test GLM on each risk/prot HLA ##################
##########################################################

rnahla2 <- rnahla
newrisk <- unlist(lapply(c(riskt,prott),function(x) grep(x,colnames(rnahla),value = T)[1]))
newrisk <- newrisk[!is.na(newrisk)]

# newrisk <- c(riskt,prott)

calcLMM <- function(x,k){
  hlaindex <- grep("HLA",colnames(x),value = T)
  tcrindex <- grep("TR",colnames(x),value = T)
  x <- x %>% pivot_longer(cols = hlaindex,values_to = "Identity",names_to = "HLA")
  splitx <- split(x,x$HLA)
  df <- NULL
  ls <- list(NULL)
  for(i in 1:length(splitx)){
    emat <- splitx[[i]]
    if(length(unique(emat$response)) < 2){
      next
    }
    mat <- as.matrix(emat[,c(5:ncol(emat)-2,ncol(emat))])
    lim <- cv.glmnet(x = mat,y = as.numeric(as.matrix(emat[,1])),family = "binomial",alpha = 1)
    pred <- predict(lim,mat,s = "lambda.1se",type = "response")[,1]
    coefs <- getCoefs(lim)
    result <- within(emat[,!colnames(emat) %in% tcrindex], Pred <- pred)
    res <- rbindlist(lapply(split(result,result$Identity),function(y) data.frame(Acc = sum((y$Pred > 5) == y$response)/nrow(y))),idcol = "present")
    comp <- result$HLA[1]
    df <- rbind(df,within(res,Compare <- comp))
    ls[[comp]] <- coefs
  }
  df <- within(df,Allele <- k)
  return(list(Acc = df, Coef = ls))
}
# result <- calcLMM(x,k)

fullresult <- list(NULL)#data.frame(NULL)
for(k in newrisk){
  print(k)
  sub <- hla[which(hla[,colnames(hla) == k] > 0),]
  tcrindex <- grep("TR",colnames(d.df),value = T)
  x <- cbind(rna.df[,1:3],d.df[,colnames(d.df) %in% c(tcrindex,newrisk)])
  x <- x[rownames(sub),!colnames(x) == k]
  result <- tryCatch(calcLMM(x,k),error = function(e){return(NA)})
  fullresult[[k]] <- result
  
  # fullresult <- rbind(fullresult,result)
}

fullresult <- fullresult[-1]
combdf <- rbindlist(lapply(fullresult,function(x) if("Acc" %in% names(x)){ x$Acc }))

ggplot(combdf[combdf$present != 2,],aes(x = factor(Compare,levels = newrisk), y = factor(Allele,levels = newrisk),fill = Acc)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +  # Add white borders between tiles
  scale_fill_gradient(low = "white", high = "red") +     # Set custom color gradient
  coord_fixed() +
  facet_wrap(~present,nrow = 2) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1),axis.text.y = element_text(size = 12))

ggsave(".//Results/IRepertoire/TCRfigs/Fig6HLA_HLAplot.svg",width = 7,height = 12)

calcTop <- function(x,tcrgenes){
  dat <- data.frame(Term = tcrgenes,Num = 0)
  rownames(dat) <- tcrgenes
  top <- lapply(x[-1],function(x) as.character(x[order(abs(x$Coef),decreasing = T),]$Term[1:10] ))
  for(i in top){
    dat[i,]$Num <- dat[i,]$Num + 1
  }
  return(dat[!is.na(dat$Num),])
}
tcrindex <- grep("TR",colnames(d.df),value = T)
coefdf <- rbindlist(lapply(fullresult,function(x) if("Coef" %in% names(x)){ calcTop(x$Coef,tcrindex)}),idcol = "HLA") 

ggplot(coefdf,aes(x = Term, y = factor(HLA,levels = c(riskt,prott)),fill = Num)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +  # Add white borders between tiles
  scale_fill_gradient(low = "white", high = "red") +     # Set custom color gradient
  coord_fixed() +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1),axis.text.y = element_text(size = 12))

coefmat <- (coefdf %>% pivot_wider(names_from = "Term",values_from = "Num"))
rown <- as.data.frame(coefmat)[,1]
coefmat <- as.matrix(coefmat[,-1])
rownames(coefmat) <- as.character(rown)

coefmat <- coefmat[match(c(riskt,prott),rownames(coefmat)),]
coefmat <- coefmat[!is.na(coefmat[,1]),]
ann <- data.frame(Risk = c(rep("risk",13),rep("protective",17)))
rownames(ann) <- rownames(coefmat)

pheatmap(t(coefmat[,colSums(coefmat)>1]),cluster_cols = F,annotation_col = ann)

