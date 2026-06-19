#Code to process the files and extract the stats and 
# install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))
require(ggplot2)
require(data.table)
require(ggpubr)
require(rstatix)
require(tidyr)
require(ggfortify)
require(umap)
require(rlist)
require(dplyr)
require(immunarch)
require(Matrix)
require(stringdist)
require(e1071)

folders <- c("anderson_01","gomez_01","mitchell_01","seay_01","IR-T1D-01","IR-T1D-02","IR-T1D-03")
files <- paste("~/Projects/T1D_RNA/repertoire/data/Sequences/",folders,"/t1d.tsv",sep="")
allseq <- data.frame(NULL)
stats <- data.frame(NULL)

for(x in files){
  print(paste("grabbing",x))
  tmp <- fread(file = x,header = T,sep = "\t")[,c(10,11,14,22,146,151)]
  print(dim(tmp))
  allseq <- rbind(allseq,tmp)
  tmp.ls <- split(tmp,tmp$repertoire_id)
  names <- names(tmp.ls)
  unique <- unlist(lapply(tmp.ls,nrow))
  total <- unlist(lapply(tmp.ls,function(y) sum(y$duplicate_count)))
  statmp <- data.frame(Subject = names,Unique = unique,Total = total,Study = rep(x,length(names)))
  stats <- rbind(stats,statmp)
  rm(statmp)
  rm(tmp.ls)
  rm(tmp)
  gc()
}

print("Found")
print(length(unique(allseq$repertoire_id)))
print("samples")

# allseq <- allseq[complete.cases(allseq),]
allseq <- allseq[!grepl("[^A-Za-z]",allseq$junction_aa),]
allseq <- allseq[!nchar(allseq$junction_aa) < 3,]

print("Remaining")
print(length(unique(allseq$repertoire_id)))
print("samples")

#Get a list of all unique CDR3
cdr3 <- distinct(allseq[,c(2,3,4)]) #unique(allseq$juction_aa)
cdr3 <- cdr3[complete.cases(cdr3),]
print(paste(nrow(cdr3),"unique cdr3"))
write.table(cdr3[,c(3,1,2)],file = "~/Projects/T1D_RNA/repertoire/data/Sequences/All/cdr3.tsv",sep = "\t",quote = F,row.names = F)

#make a sparse dataset containing which patients have it
sparse <- sparseMatrix(i = as.numeric(as.factor(allseq$junction_aa)),j = as.numeric(as.factor(allseq$repertoire_id)),x = allseq$duplicate_count)
dimnames(sparse) <- list(levels(as.factor(allseq$junction_aa)),levels(as.factor(allseq$repertoire_id)))
saveRDS(sparse, file = "~/Projects/T1D_RNA/repertoire/sparse_repertoire.rds")
#Now we can delete the original because we got all the info we need
rm(allseq)
gc()


#Make a cluster sparse dataset

clonal <- readRDS(file = "~/Projects/T1D_RNA/repertoire/sparse_repertoire.rds")
clusters <- read.table(file = "~/Projects/T1D_RNA/repertoire/data/Sequences/All/cdr3--RotationEncodingBL62.txt",header = F,sep = "\t")
clusters <- clusters[,1:2]
colnames(clusters) <- c("CDR3","Cluster")
clusters <- clusters[!duplicated(clusters$CDR3),]
clusters <- clusters[order(clusters$Cluster),]
numclust <- table(clusters)

subclonal <- clonal[clusters$CDR3,]

#Now let's go ahead and cluster/get targets for the sequences

#Calc the distance between every CDR3 efficiently


# cdr3dat <- list(data = list(dat = cdr3),meta = data.table(Sample = "Dat",Lane = "A"))
# 
# distTCR <- seqDist( cdr3dat$data, .col = 'junction_aa',.methods = "lv", .group_by = c('v_call'),.group_by_seqLength = T)
# clustTCR <- seqCluster(cdr3dat$data, distTCR, .perc_similarity = 0.9)
# print("found")
# print(clustTCR$"A2-i129" %>% .$Cluster %>% unique() %>% length())
# print("unique clusters")
# 
# vdjdb <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = "HomoSapiens")
# ann <- vdjdb[vdjdb$cdr3 %in% cdr3$junction_aa,]

#What is the average levenstein distance between same target
# vdjls <- split(vdjdb, vdjdb$antigen.epitope)
# length(vdjls)
# vdjls <- vdjls[unlist(lapply(vdjls,nrow)) > 1]
# length(vdjls)

#Make a function to compare and get levenstein for each pair also calc length distance and diff in chain, v, and j seg
# calc_diffs <- function(df){
#   print(df$antigen.epitope[1])
#   df <- df[complete.cases(df),]
#   if (nrow(df)<2 | nrow(df)>200) {
#     return(NULL)
#   }
#   combs <- combn(c(1:nrow(df)), 2)
#   #randomly sample up to 1000 diffs
#   if (ncol(combs) > 1000){
#     combs <- combs[,sample(1:ncol(combs),1000)]
#   }
#   print(paste("Processing",ncol(combs),"pairs"))
#   dat <- data.frame(NULL)
#   for (i in 1:ncol(combs)){
#     iter <- combs[,i]
#     x <- df[iter[1],]
#     y <- df[iter[2],]
#     lev <- stringdist(x$cdr3,y$cdr3,method = "lv")
#     # lev <- hamming.distance(x$cdr3,y$cdr3)
#     ldiff <- abs(nchar(x$cdr3) - nchar(y$cdr3))
#     all_a <- x$gene == "TRA" & y$gene == "TRA"
#     all_b <- x$gene == "TRB" & y$gene == "TRB"
#     chain_match <- x$gene == y$gene
#     v_match <- x$v.segm == y$v.segm
#     j_match <- x$j.segm == y$j.segm
#     mhc <- x$mhc.class
#     hla <- x$mhc.a
#     epitope <- x$antigen.epitope
#     dat <- rbind(dat,data.table(x = x$cdr3,y = y$cdr3,
#                                 dist = lev, len_diff = ldiff, allA = all_a, allB = all_b, chaim_match = chain_match, 
#                                 v_match = v_match, j_match = j_match,mhc = mhc, hla = hla,epitope = epitope))
#   }
#   return(dat)
# }
# vdj_dist <- lapply(vdjls,calc_diffs)
# vdj_dist <- rbindlist(vdj_dist)
# random <- calc_diffs(vdjdb[sample(1:nrow(vdjdb),200),])
# chain <- rep("Cross",nrow(vdj_dist))
# chain[vdj_dist$allA == T] <- "TCRA"
# chain[vdj_dist$allB == T] <- "TCRB"
# vdj_dist <- within(vdj_dist,Chain <- chain)
# random <- within(random,Chain <- rep("Random",nrow(random)))
# 
# vdj <- rbind(vdj_dist,random)
# 
# ggplot(vdj[vdj$Chain != "Cross",],aes(x = dist,fill = Chain)) +
# #  geom_boxplot()
#   geom_density(alpha = 0.5) +
#   facet_grid(v_match~j_match)
# 
# ggplot(vdj_dist[vdj_dist$chain_match == T,],aes(x = epitope,y = dist)) +
#   geom_boxplot()
# 
# 
# 
# ggplot(vdj[vdj$Chain != "Cross",],aes(x = dist,y = len_diff,colour = Chain)) +
#   geom_jitter(alpha = 0.5) +
#   facet_grid(v_match~j_match)
# 
# ggplot(vdj[vdj$Chain != "Cross" & len_diff == 0,],aes(x = dist,fill = Chain)) +
#   #  geom_boxplot()
#   geom_density(alpha = 0.5) +
#   facet_grid(v_match~j_match) +
#   ggtitle("Lv dist for equal length")

#Giana code
#

  
