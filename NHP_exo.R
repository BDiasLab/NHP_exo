 
### LOADING PACKAGES
library(reshape2)
library(gridExtra)
library(Rsamtools)
library(Rsubread)
library(ggplot2)
library(matlab)
library(matrixStats)
library(dendextend)
library(tidyverse) 
library(edgeR) 
library(matrixStats) 
library(gt)
library(limma)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(plyr)
`%!in%` <- Negate(`%in%`)

### LOADING DATA 
setwd("NHP_data")
#only sense reads were analyzed
star_sense = read.csv("STAR_sense.csv")
#removing first rows of descriptive statistics
star_sense = star_sense[5:35436,]

#loading metadata
meta_monkey =  read.csv("meta_monkey.csv")
meta_monkey$meta = paste(meta_monkey$sample_short, meta_monkey$group, sep = "_") 
sampleLabels = meta_monkey$meta

### reads with coverage higher than 80%
survived_transcripts = read.csv("survived_transcripts.csv")
survived_transcripts = survived_transcripts$x
length(survived_transcripts)

stat_survived = read.csv("filtered_genes_count.csv")
ggplot(data=stat_survived, aes(x=as.factor(case), y=fraction), color = "#949494")  +geom_bar(stat = "identity")+ ggtitle("Reads with coverage higer than 80%")+ theme(text = element_text(size = 20), axis.text.x = element_text(size=20, angle=90)) + ylab("% of reades")+geom_text(aes(label = paste0(round(fraction), "%")), position = position_stack(vjust = 0.5)) + xlab("case ID")

###mapping statistics obtained from samtools
mapping_stat = read.csv("mapping.csv")
ggplot(data=mapping_stat, aes(x=as.factor(case), y=percentage, fill = reads))  +geom_bar(stat = "identity") +geom_text(aes(label = paste0( round(percentage), "%")), 
                                                                                                                       position = position_stack(vjust = 0.5))+scale_fill_manual(values=c( "light grey","grey", "dark grey"))+ theme(text = element_text(size = 15), axis.text.x = element_text(angle=90))  +ylab("% of reads") +xlab("case ID") + ggtitle("Mapping statistics")

star_sense = star_sense[star_sense$X %in% survived_transcripts,]
star_sense_tr = star_sense %>% dplyr::group_by(X, Dias17,Dias18,Dias19,Dias20,Dias21,Dias22,Dias23,Dias24,Dias25,Dias26,   
                                               Dias27,Dias28,Dias29)

row.names(star_sense) = star_sense$X
star_sense$X = NULL
names(star_sense) = c("Dias17","Dias18","Dias19","Dias20","Dias21","Dias22","Dias23","Dias24","Dias25","Dias26",   
                      "Dias27","Dias28","Dias29")

### PERFORMING DGE (Differential gene expression) ANALYSIS
myDGEList_star_strand <- DGEList(star_sense)

# counts per million
cpm_star_strand <- cpm(myDGEList_star_strand) 
colSums(cpm_star_strand)
log2.cpm_star_strand <- cpm(myDGEList_star_strand, log=TRUE)
log2.cpm.df_star_strand <- as_tibble(log2.cpm_star_strand, rownames = "geneID")
sampleLabels_star_strand = names(log2.cpm.df_star_strand)[-c(1)]

log2.cpm.df_star_strand = as.data.frame(log2.cpm.df_star_strand)
names(log2.cpm.df_star_strand) = c("geneID", sampleLabels)
log2.cpm.df.pivot_star_strand <- pivot_longer(log2.cpm.df_star_strand, # dataframe to be pivoted
                                              cols = as.vector(unique(sampleLabels)),
                                              names_to = "samples", 
                                              values_to = "expression") 

## filtering 
table(rowSums(myDGEList_star_strand$counts==0)==10)
keepers_star_strand <- rowSums(cpm_star_strand>1)>=6
keepers_to_filter_malt <- rowSums(cpm_star_strand>1)>=3
table(rowSums(cpm_star_strand>1)>=3)

myDGEList.filtered_star_strand <- myDGEList_star_strand[keepers_star_strand,]
dim(myDGEList.filtered_star_strand)
dim(myDGEList_star_strand) 

to_filter_genes = myDGEList_star_strand[keepers_to_filter_malt,]
to_filter_genes = to_filter_genes$counts

row.names(myDGEList.filtered_star_strand$samples) = sampleLabels
row.names(myDGEList_star_strand$samples) = sampleLabels

## TMM normalization
myDGEList.filtered.norm_star_strand <- calcNormFactors(myDGEList.filtered_star_strand, method = "TMM")


### PERFORMING VOOM
# Transform RNA-Seq Data Ready for Linear Modelling
# Filtering most differentially expressed genes
###### VOOM 
group  <- factor(meta_monkey$group)
sex<- factor(meta_monkey$Sex)

design <- model.matrix(~group*sex)

v.DEGList.filtered_star_strand <- voom(myDGEList.filtered.norm_star_strand, design, plot = FALSE, normalize.method ="none")
fit_star_strand <- lmFit(v.DEGList.filtered_star_strand, design)
fit2_star_strand <- eBayes(fit_star_strand)

TopHits_star_strand_group <- topTable(fit2_star_strand, adjust ="fdr", coef=2, number=1000, sort.by="p" )
head(TopHits_star_strand_group)
dim(TopHits_star_strand_group)

myTopHits.df_star_strand_group_filt =  TopHits_star_strand_group[TopHits_star_strand_group$P.Value < 0.05,]
myTopHits.df_star_strand_group_filt$logabs= abs(myTopHits.df_star_strand_group_filt$logFC)
myTopHits.df_star_strand_group_filt = myTopHits.df_star_strand_group_filt[myTopHits.df_star_strand_group_filt$logabs > 2,]
myTopHits.df_star_strand_group_filt$geneID = row.names(myTopHits.df_star_strand_group_filt)
dim(myTopHits.df_star_strand_group_filt)
myTopHits.df_star_strand_group_filt[1:20,]

sign_genes = myTopHits.df_star_strand_group_filt[myTopHits.df_star_strand_group_filt$adj.P.Val < 0.05,]$geneID

v.DEGList.filtered_star_strand_sign = v.DEGList.filtered_star_strand$E
v.DEGList.filtered_star_strand_sign = as.data.frame(v.DEGList.filtered_star_strand_sign)
v.DEGList.filtered_star_strand_sign$geneID = row.names(v.DEGList.filtered_star_strand_sign)
v.DEGList.filtered_star_strand_sign = v.DEGList.filtered_star_strand_sign[v.DEGList.filtered_star_strand_sign$geneID %in% sign_genes,]

### visualizing survived transcripts 
genes_sign = read.csv("genes_sign.csv")
line_size = 1
ggplot(data=genes_sign, aes(x=as.factor(geneID), y=logFC, fill=group)) +
  stat_summary(
    fun.data = mean_se,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.3,
    size = line_size, position=position_dodge(width=0.5)) +stat_summary(
      fun.data = mean_se,
      fun.args = list(mult = 1),
      geom = "bar",position=position_dodge(width=0.5),
      width = 0.5)+ labs(title = "logfc sign genes", x ="", y = "logfc")  + theme(text = element_text(size = 15))  +scale_fill_manual(values = c("#0080FF", "#FF8886"))

### PERFORMING GO ANALYSIS 
# Functional profiles for genes and gene clusters
library(clusterProfiler)
library(enrichplot)
library(biomaRt)

df_qbio = myTopHits.df_star_strand_group_filt
original_gene_list <- df_qbio$logFC

# name the vector
names(original_gene_list) <- df_qbio$geneID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


names = names(gene_list)

library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
org = get("org.Hs.eg.db")

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             #keyType = "ENSEMBL", 
             #keyType = "GENENAME",
             keyType = "SYMBOL",
             minGSSize = 5, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org, 
             pAdjustMethod = "fdr")

#read dataset with filtered unique GO terms 
gse_qbio = read.csv("gse_qbio_names_oct22_BP.csv")

row.names(gse_qbio) = gse_qbio$ID
gse@result =gse_qbio 

require(DOSE)
enrichplot::dotplot(gse, x = "GeneRatio",
                    color = "p.adjust", font.size = 15, label_format = 10, 
                    title="Gene Ontology Enrichment (GO) Analyses", showCategory = 12) + facet_grid(.~.sign)

## changing the order of presented groups 
gse@result$num = c(1,5,2,8,10,3,7,12,6,11,4,9)
enrichplot::dotplot(gse, x = "GeneRatio", orderBy = "num",
                    color = "p.adjust", font.size = 15, label_format = 20, 
                    title="Gene Ontology Enrichment (GO) Analyses", showCategory = 12) + facet_grid(.~.sign)


heatplot(gse, foldChange=gene_list)+ggplot2::theme(axis.text.x = element_text(angle = 90, size = 6), axis.text.y = element_text(size = 20))+scale_fill_gradient(low = "orange",high = "#570861",guide = "colorbar")  

### BRACKEN METAGENOMEANALYSIS 

#loading data 
#plotting a rarefaction curve of one of the samples
brack = read.csv("21.bracken_f.csv")

cat ("Total number of reads is ",sum(brack$new_est_reads), "\n") 
cat ("The number of unique species is ",length(unique(brack$taxonomy_id)),"\n")

#this function plots the curve of reads sampling over # of detected species 
sampling = function (data, fraction_of_reads){
  species_num = rep(0,length(fraction_of_reads))
  reads = sum(data$new_est_reads)
  prob_sp = data$new_est_reads/reads
  num = 0
  for (i in fraction_of_reads){
    num = num +1
    
    pois = rpois(nrow(data),lambda= prob_sp*i*reads)
    positive = length(pois[pois > 0])
    #print(positive)
    species_num[num] = positive 
    
    cat("For fraction", i, " estimated number of different species is ", positive, "\n")
  }
  
  plot(fraction_of_reads, species_num)
  
}

fraction_of_reads = c(0.0001, 0.001, 0.01, 0.1, 0.2, 0.5)
data = brack

#applying the function with given fraction parameters
sampling(brack, fraction_of_reads)

# Maltreated vs Control METAGENOME ANALYSIS using Bracken output

#loading bracken full output dataframe with metadata
brack.df = read.csv("monkey_bracken.csv")

#### SAMPLING UNIQUE

unique_control_all = c()
unique_malt_all = c()

control_bac_all = brack.df[brack.df$group=="control" & brack.df$fraction_total_reads >0,]
malt_bac_all = brack.df[brack.df$group=="malt" & brack.df$fraction_total_reads >0,]

for (i in unique(control_bac_all$name)) 
{ 
  sub.df = control_bac_all[control_bac_all$name==i & control_bac_all$fraction_total_reads>0 ,]
  
  if (nrow(sub.df) >= 7) {
    unique_control_all <- c(unique_control_all, i)
  }
}


for (i in unique(malt_bac_all$name)) 
{ 
  sub.df = malt_bac_all[malt_bac_all$name==i & malt_bac_all$fraction_total_reads>0,]
  
  if (nrow(sub.df) >= 6) {
    unique_malt_all <- c(unique_malt_all, i)
  }
}

malt_bac_all_name=unique(malt_bac_all[malt_bac_all$fraction_total_reads>0,]$name) 

cont_bac_all_name=unique(control_bac_all[control_bac_all$fraction_total_reads>0,]$name) 



##### venn.diagram of overlapped speccies

venn.diagram(
  x = list(unique_control_all, unique_malt_all),
  category.names = c("CONT", "MALT"),
  filename = 'venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 680 , 
  width = 720 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("red", "blue"),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans", # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

common_all = unique_control_all[unique_control_all %in% unique_malt_all]


#Defining distinct classes of microbiome profile in control and maltreated groups

# this step takes long time. you can access the data from a file 
classification_common_all_total = read.csv("classification_common_all_total.csv")


#or run the following code: 
#classification_common_all = unique(brack.df[brack.df$fraction_total_reads >0,]$name)
# classification_common_all_total = NULL
# 
# for (i in classification_common_all) {
#   data_spec = myTAI::taxonomy(organism = i, db= "ncbi",output   = "classification")
#   data_spec$species = i
#   classification_common_all = rbind(classification_common_all, data_spec)
# }


classification_common_all_total$count = 1
classification_common_all_total_class = classification_common_all_total[classification_common_all_total$rank=="class",]
total_classification_full_common = classification_common_all_total_class %>% group_by(name) %>% summarise(count = sum(count))
total_classification_full_common_filt = total_classification_full_common[total_classification_full_common$count>10,]

classification_common_all_total_class_filt = classification_common_all_total_class[classification_common_all_total_class$name %in% total_classification_full_common_filt$name ,]


#plotting common for CONT and MALT microbiome classes 
df.pie = total_classification_full_common_filt %>% 
  arrange(desc(name)) %>%
  mutate(prop = count / sum(total_classification_full_common_filt$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


ggplot(df.pie, aes(x="", y=count, fill=name)) +
  geom_bar(stat="identity", width=2,  color="white") + coord_polar("y", start=0) + geom_text(aes(x = 1.6, label = count), position = position_stack(vjust = .5), color = "white", size=4) +theme_bw()+theme(text = element_text(size = 15))


# plotting unique for CONT animals microbiom classes 

classification_cont_perm_class = read.csv("classification_cont_perm_class.csv")

ggplot(classification_cont_perm_class[classification_cont_perm_class$count>2,], aes(x="", y=count, fill=name)) +
  geom_bar(stat="identity", width=2,  color="white") + coord_polar("y", start=0) + geom_text(aes(x = 1.6, label = count, color = name), position = position_stack(vjust = .5), color = "white", size=4) +theme_bw()+theme(text = element_text(size = 15))


classification_common_all_total_phylum = classification_common_all_total[classification_common_all_total$rank=="phylum",]
total_classification_full_common_phylum = classification_common_all_total_phylum %>% group_by(name) %>% summarise(count = sum(count))

unique(classification_common_all_total[classification_common_all_total$rank=="phylum",]$name)


### PLOTTING bracken output stats for aligned not aligned to microbiome reads
##### mapped vs unmapped 
homo_stat = read.csv("homo_stat.csv")
ggplot(data=homo_stat, aes(x=as.factor(case), fill=filtering, y=mapped)) +geom_bar(stat="identity", position = "dodge") + labs(x = "animal ID", y = "reads mapped to bacteria, %") + scale_fill_manual(values = c("lightgrey", "darkgrey")) + theme(text = element_text(size = 20))+geom_text(aes(label = paste(round(mapped,digits=0), "%")),size=4, color = "black")+ggtitle("Reads aligned to microbiome")


### richness comparison, median
bac_richness = read.csv("bac_richness.csv")
bac_richness = bac_richness[bac_richness$name != "Homo sapiens",]
bac_richness = bac_richness[bac_richness$fraction_total_reads > 0.000001,]
meta_micro = read.csv("meta_monkey_micro.csv")
bac_richness = merge(bac_richness, meta_micro)

mu_total = ddply(bac_richness, "group", summarise, fraction.median=median(log))

ggplot(bac_richness, aes(x=log, color=group, fill=group)) +geom_density(alpha=0.6)+
  geom_vline(data=mu_total, aes(xintercept=fraction.median, color=group),
             linetype="dashed")+
  scale_color_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  labs(title="Species density",x="Fraction, log", y = "Density")+
  theme_classic() + theme(text = element_text(size = 15))


#### enthropy median mean 

library(DescTools)
library(modeest)
library(lsr)

sum.df = NULL
cases = unique(bac_richness$case)
win = 0.25
for (i in 1:13 ){
  c = cases[i]
  print(c)
  vals = bac_richness[bac_richness$case==c,]$log
  out = data.frame(case=c)
  out$entropy = Entropy(vals)
  out$mean = mean(vals)
  out$sd = sd(vals)
  out$median = median(vals)
  out$iqr = IQR(vals)
  out$mode = mlv(vals, method = "meanshift")
  out$density = length(vals[vals > out$mode - win & vals < out$mode + win])
  sum.df = rbind(sum.df, out)
}

sum.df = merge(sum.df, meta_micro)
sum.df$group_shuffle=NULL

wilcox.test(sum.df$mode ~ sum.df$group)
t.test(sum.df$mode ~ sum.df$group)
cohensD(sum.df$mode[sum.df$group=="control"], sum.df$mode[sum.df$group=="malt"])


ggplot(sum.df, aes(x=group, y=mode)) + geom_boxplot(aes(fill=group), alpha = 0.5)+ggtitle("Microbiome RNA abundancy") + stat_compare_means(size=5.5,label.y = -4.1) +theme(axis.text.x=element_text(angle=90, hjust=1))  + theme(text = element_text(size = 20)) +ylab("mode of reads fraction, log") + scale_fill_manual(values = c("blue", "red" ),name = "Group", labels = c("CAU", "MALT")) + theme(text = element_text(size = 20))

ggplot(bac_richness,aes(x=log, fill=group, color = group, group=case))+ geom_density(alpha = 0.1, size=1)+ggtitle("Species density") +theme(axis.text.x=element_text(angle=90, hjust=1))  + theme(text = element_text(size = 20)) +ylab("Density") + xlab("Fraction, log")+scale_fill_manual(values = c("blue", "red" )) + scale_color_manual(values = c("blue", "red"), name = "Group", labels = c("CAU", "MALT"))+theme(text = element_text(size = 20))



#COMPARING GUT AND EXO (phyli)

### list of genus  
classification_genus_clean = read.csv("classification_genus_upd.csv")
names(classification_genus_clean) = c("name_rank","rank", "name")
unique(classification_genus_clean[classification_genus_clean$rank=="phylum",]$name_rank)

#filtering median reads 
brack.df2 = brack.df[brack.df$fraction_total_reads > median(brack.df$fraction_total_reads),]

brack.df2_classification = merge(x=brack.df2, y=classification_genus_clean, by=c("name"), all.x=TRUE)

brack.df2_classification_up = read.csv("brack.df2_classification_phylum_class.csv")
brack.df2_phylum = brack.df2_classification_up[brack.df2_classification_up$rank=="phylum",]

exo_phylum=read.csv("exo_phylum.csv")
exo_phylum$sample = "exo"
exo_class = read.csv("exo_class.csv")

classification_gut = read.csv("gut_micro_counts_class_phylum.csv")

data_long <- tidyr::gather(classification_gut, case, fraction_total_reads, X22.unmapped:X20.unmapped, factor_key=TRUE)

classification_gut_long = read.csv("classification_gut_long.csv")

classification_gut_long = merge(classification_gut_long,meta_micro)

gut_phylum=read.csv("gut_phylum.csv")

gut_phylum_sum = gut_phylum %>% dplyr::group_by(Phylum, group, sample) %>% dplyr::summarise(fraction_total_reads = sum(fraction_total_reads))
gut_phylum_sum_upd_mean = gut_phylum %>% dplyr::group_by(Phylum, case, group, sample) %>% dplyr::summarise(fraction_total_reads = sum(fraction_total_reads))

exo_class_sum = exo_class %>% dplyr::group_by(Class, group, sample) %>% dplyr::summarise(fraction_total_reads = sum(fraction_total_reads))

exo_phylum_sum = exo_phylum %>% dplyr::group_by(Phylum, group, sample) %>% dplyr::summarise(fraction_total_reads = sum(fraction_total_reads))
phylum_sum = read.csv("phylum_sum_upd.csv")

#### plotting phylum

library(RColorBrewer)
library(ggsci)
ggplot(phylum_sum, aes(x = group, y = fraction_norm, fill = Phylum)) + geom_bar(stat = "identity")+facet_grid(.~sample) + theme(text = element_text(size = 15)) +ylab("Normalized abundance") +xlab("") +scale_fill_manual(values = pal_jco()(8))


