Skip to content
This repository
Search
Pull requests
Issues
Marketplace
Explore
 @sargdavid
 Sign out
 Watch 0
  Star 0  Fork 0 renyiwu/bioseq
 Code  Issues 0  Pull requests 0  Projects 0  Wiki  Insights
Branch: master Find file Copy pathbioseq/RNA-seq_general_hisat2.sh
26f687d  on Sep 27
@renyiwu renyiwu new
1 contributor
RawBlameHistory    
116 lines (95 sloc)  4.61 KB
#RNAs-seq workflow. R Wu Sept 2017
#Hisat2 -> picard -> feartureCounts
#########################
#This version uses Hisat2
#########################
#
#1, align reads to genome with hisat2 and save aligns to bam file with samtools. 
#2, optional, concatnate multiple files if they belong to one sample.
#3, sort .bam files
#4, deduplicate with picard.
#5, get genes counts file with featureCounts. use annotation file here.


#0.2 Rename (shorten) fastq file names if needed.
rename 's/_S.*gz/.fastq.gz/' *.fastq.gz
#eg, C75_S1_R1_001_.fastq.gz to C75.fastq.gz

#1 hisat2
for i in *.fastq.gz; do hisat2 -p 6 -x /path/to/hisat2/index/ -U $i | samtools view -bh -o ${i%.fastq.gz}.bam -; done #note: piping to samtolls reqires version 1.15 or above. The last "-" (dash) may be ommitted.

#2 concatnate
samtools cat -o RW04.bam RW04?.bam

 ~/bin/samtools1 sort -@ 6 -o WJ1.sorted.bam WJ1.bam

#3 sort
for i in *.bam; do ~/bin/samtools1 sort -@ 6 -o ${i%.bam}.sorted.bam $i; done

#4 deduplicate
#picard
for i in *.sorted.bam; do java -jar ~/tools/picard-2.10.5/picard.jar MarkDuplicates I=$i O=${i%.sorted.bam}.dedup.bam M=${i%.sorted.bam}.dedup.txt REMOVE_DUPLICATES=true; done
#or samtools
for i in *.sorted.bam; do samtools rmdup -s $i ${i%.sorted.bam}.rmdup.bam; done

#5 FeatureCounts
featureCounts --primary -T 8 -a /path/to/genes.gtf -o featurecounts.results.csv *dup.bam


###### one line command for all:

rename 's/_S.*gz/.fastq.gz/' *.fastq.gz && for i in *.fastq.gz; do hisat2 -p 6 -x ~/genome/mm10/genome -U $i | samtools view -bh -o ${i%.fastq.gz}.bam -; samtools sort -@ 6 -o ${i%.fastq.gz}.sorted.bam ${i%.fastq.gz}.bam; rm ${i%.fastq.gz}.bam; samtools rmdup -s ${i%.fastq.gz}.sorted.bam ${i%.fastq.gz}.rmdup.bam; rm ${i%.fastq.gz}.sorted.bam; done && featureCounts --primary -T 8 -a ~/genome/mm10/genes.gtf -o featurecounts.results.csv *dup.bam

####################
#   Version 2      #
####################
#
#Rename (shorten) fastq file names if needed.
rename 's/_S.*gz/.fastq.gz/' *.fastq.gz
#eg, C75_S1_R1_001_.fastq.gz to C75.fastq.gz
#
#concatnate. eg, RW041 and RW042.
zcat RW04?.fastq.gz | gzip >RW04.fastq.gz

#Count reads
for i in RW??.fastq.gz; do zcat $i | echo $((`wc -l`/4)); done

#run
for i in *.fastq.gz; do hisat2 -p 6 -x ~/genome/mm10/genome -U $i | samtools sort -@ 6 --output-fmt SAM | samtools rmdup -s --output-fmt BAM - ${i%.fastq.gz}.rmdup.bam; done && featureCounts --primary -T 8 -a ~/genome/mm10/genes.gtf -o featurecounts.results.csv *dup.bam

#########################################################
#Further analyses with counting results by featureCounts:
#Run in R
R
library(DESeq2)
library(data.table)
library("gplots")
# load matrix of read counts
df <- read.table("/home/administrator/Documents/bgi_mes13/mes13_featurecounts_Dec2017_david.csv", header = T, row.names = "Geneid")
df <- df[,-(1:4)]
colnames(df) <- c("length", paste("WJ", 1:7, sep = ""))
df_7 <- df[-1]
# create sample matrix
mat <- matrix(paste("WJ",1:7, sep = ""), nrow = 7)
colnames(mat) <- "condition"
mat
rownames(mat) <- colnames(df[-1])
mat_7 <- mat
dds <- DESeqDataSetFromMatrix(countData=df_7, colData=mat_7, design= ~ condition)
dds <- dds[ rowSums(counts(dds) >= 50) >= 2, ]  
mcols(dds)$basepairs <- df$length
dds_7 <- DESeq(dds)

write.table(fpkm(dds_7), file = "/home/administrator/Documents/bgi_mes13/mes13_fpkm_Dec2017_david.csv", sep = "\t", quote = T, col.names = NA)

#
fpkm_7 <- fpkm(dds_7)
fpkm_71 <- fpkm_7+1 #Add pseudo 1 to all values to avoid log10(0) error.
fpkm_72 <- log10(fpkm_7+1) #log10 transformation
mat_f <- fpkm_72[!rowSums(fpkm_72 <1),] #Only keep rows which have values greater than or equal to 2 (in all groups/columns)
write.table(log10(fpkm(dds_7)+1), file = "data/shan_rna/RW_all_log10_fpkm.csv", sep = "\t", quote = T, col.names = NA)
write.table(mat_f, file = "data/shan_rna/RW_all_log10_fpkm_gt1.csv", sep = "\t", quote = T, col.names = NA)

#tmp <- read.table(file = "data/shan_rna/RW_all_log10_fpkm_gt11.csv", sep = "\t", header = T, row.names = "X")
#heatmap
lmat <- rbind(c(4,1),c(3,1),c(2,1))
lmat
dev.off()
heatmap.2(mat_f,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat,col=bluered(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none", key = FALSE, labCol = NA,
          trace = "none",
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f)


dds$condition <- relevel(dds$condition, ref="WJ1")
# filter for genes with at least one count in at least two samples:
dds <- dds[ rowSums(counts(dds) >= 50) >= 2, ]  
#Run deseq
dds_7 <- DESeq(dds)
# pairwise comparisons
res1 <- results(dds_7,contrast=c("condition", "WJ2", "WJ1"))
#write out
write.table(res1, file="/home/administrator/Documents/bgi_mes13/WJ2-WJ1-g50.csv", sep="\t", quote=T, col.names=NA)

© 2017 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
API
Training
Shop
Blog
About