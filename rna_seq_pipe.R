#library calls
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(Biostrings)
library(GOstats)
#library(BiocParallel)
library(ape)
library(edgeR)
library(AnnotationDbi)
library(DESeq2)
source("~/rnaSeq_scripts/overlapper.R")
source("~/rnaSeq_scripts/utilities.R")
source("~/rnaSeq_scripts/geneSetAnalysis.R")
#fastxtoolkit
fastx = function(name, fq){
	stats = paste("fastx_quality_stats -i ",fq," -o ./quality_analysis/",name,".txt", sep="")
	system(stats)
	box = paste("fastq_quality_boxplot_graph.sh -i ./quality_analysis/",name,".txt -o ./quality_analysis/",name,"_boxplot.png -t ",name, sep="")
	system(box)
	nuc = paste("fastx_nucleotide_distribution_graph.sh -i ./quality_analysis/",name,".txt -o ./quality_analysis/",name,"_nuc.png -t ",name, sep="")
	system(nuc)  

}


#call to tophat
tophat = function(name, fq1, fq2, ref, gff){
	tp2 = "~/Downloads/tophat-2.0.14.Linux_x86_64/tophat"
	sys_call = paste(tp2," -g 1 --segment-length 25 -i 30 -I 3000 -p 10 -o ",name,"_tophat ",ref," ", fq1, " ",fq2," -G ", gff, sep="")
	#print(sys_call)
	#?readline
	system(sys_call)
}
###Taken from SystemPiper by T. Girke see GitHub
returnRPKM = function(counts, ranges){
	geneLengthsInKB <- sum(width(reduce(ranges)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)

}
 
################run section for samples####################
arg = commandArgs(trailingOnly = TRUE)
#target file
targets = read.delim(arg[1], comment.char = "#")
#refrence genome
refG = arg[2]
#gff annotation file
antG = arg[3]
#sqlite file db 
## later insert check if file does not exitst create it from GFF file
dbPath = arg[4]
###check if all the varibles are all set except for 4 which can be generated from the GFF###


#if(dbPath == NULL){
#	txdb = makeTxDbFromGFF(file=antG, format="gff3", dataSource="TAIR", organism="A. thaliana")
#	need to rename the seqlevles of the txdb to match bam files
#	saveDb(txdb, file="./TAIR_trunc.sqlite") 
#}
#else{
	txdb = loadDb(dbPath)
#	txdb
#1

#need if statement to check if it already exists

#dir.create("./quality_analysis")


bamfl = c()
#while targets get line
for(i in 1:nrow(targets)){
	sample=targets[i,]$SampleLong
	fastq1=targets[i,]$FileName1
	fastq2=targets[i,]$FileName2
	if(!file.exists(paste("./quality_analysis/",sample,"_1.txt",sep=""))){
	#fastx on first fastq file
	#	fastx(paste(sample,"_1",sep=""), fastq1)
	}
	if(!file.exists(paste("./quality_analysis/",sample,"_2.txt",sep=""))){
	#fastx on second fastq file
	#	fastx(paste(sample,"_2",sep=""), fastq2)
	}
	#top hat on pair of the files
	if(!file.exists(paste("./",sample,"_tophat/accepted_hits.bam",sep=""))){
		tophat(sample,fastq1,fastq2, refG, antG)
	}
	#use the name to get bam file
	bamfl = c(bamfl, BamFile(paste("./",sample,"_tophat/accepted_hits.bam",sep=""), yieldSize=50000))
	#print(i)
}
#return RPKM

names(bamfl) = targets$SampleLong
#bamfl[1]
core_db =  renameSeqlevels(txdb, seqlevels(bamfl[[1]]))
#core_db
exonXgene = exonsBy(core_db, by="gene")

#bamfiles = BamFileList(bamfl, yieldSize=50000, index=character())
countEbyGene = lapply(bamfl, function(x) summarizeOverlaps(exonXgene, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=FALSE))

#change singele end true when doing SE vs PE sequencing
countDFeByG = sapply(seq(along=countEbyGene),function(x) assays(countEbyGene[[x]])$counts)

rownames(countDFeByG) = names(rowRanges(countEbyGene[[1]])); colnames(countDFeByG) = names(bamfl)
rpkmDFeByG = apply(countDFeByG, 2, function(x) returnRPKM(counts=x, ranges=exonXgene))
write.table(countDFeByG, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByG, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
#RPKM mean and clustering
#library(ape)
library(DESeq2, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
colData <- data.frame(row.names=targets$SampleLong, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("./results/spearman_sample_cor.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()
#use edgeR determine differential expression
#library(edgeR)
countDF <- read.delim("./results/countDFeByg.xls", row.names=1, check.names=FALSE)
targets <- read.delim(arg[1], comment="#")
#cmp <- readComp(file="targets.txt", format="matrix", delim="-")
cmp = readComp(arg[1], format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
edgeDF2 <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[2]], independent=FALSE, mdsplot="")
write.table(edgeDF2, "./results/edgeRglm_allcomp2.xls", quote=FALSE, sep="\t", col.names = NA)

library(ggplot2)
edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)

#venn diagram comparison
vennsetup <- overLapper(DEG_list$Up, type="vennsets")
vennsetdown <- overLapper(DEG_list$Down, type="vennsets")
#write.table(DEG_list$Up, "results/vensetup.txt", sep="\t")
#write.table(DEG_list$Down, "results/vensetdown.txt", sep="\t")
pdf("results/vennplot.pdf")
vennPlot(list(vennsetdown, vennsetup), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
#dev.off()

#clustering and heat maps
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
Ngeneids = geneids[-grep("ATC.*",geneids)]
Cgeneids = geneids[grep("ATC.*",geneids)]
y <- assay(rlog(dds))[geneids, ]
n = assay(rlog(dds))[Ngeneids, ]
c = assay(rlog(dds))[Cgeneids, ]
pdf("./results/heatmap_all.pdf")
pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()


#gene to go
library("biomaRt")
source("~/rnaSeq_scripts/geneSetAnalysis.R")
#listMarts() # To choose BioMart database
#m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m)
if(!file.exists("data/GO/catdb.RData")){
m <- useMart(host="plants.ensembl.org", biomart="plants_mart", dataset="athaliana_eg_gene")
listAttributes(m) # Choose data types you want to download
go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component", 3] <- "C"
go[1:4,]
dir.create("./data/GO")
write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
save(catdb, file="data/GO/catdb.RData")
}
#go enrichment
load("data/GO/catdb.RData")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1), plot=FALSE)
up_down <- DEG_list$UporDown; names(up_down) <- paste(names(up_down), "_up_down", sep="")
up <- DEG_list$Up; names(up) <- paste(names(up), "_up", sep="")
down <- DEG_list$Down; names(down) <- paste(names(down), "_down", sep="")
#DEGlist <- c(up_down, up, down)
DEGlist = c(up, down)
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type="gene", CLSZ=2, cutoff=0.9, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
library("biomaRt"); m <- useMart(host="plants.ensembl.org",biomart="plants_mart", dataset="athaliana_eg_gene")
goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type="gene", myslimv=goslimvec, CLSZ=10, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)


#plot go results ##### needs refinement
#gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
gos <- BatchResultslim
pdf("./results/GOslimbarplotMF.pdf", height=8, width=10);
 goBarplot(gos, gocat="MF")
 goBarplot(gos, gocat="BP")
 goBarplot(gos, gocat="CC")
dev.off()

#clustering and heat maps
##### Uncomment the following lines to break up large heat maps in to those that that are not in the chloroplast or mitochondria
pdf("./results/heatmap_nuc.pdf")
pheatmap(n, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()
pdf("./results/heatmap_chl.pdf")
pheatmap(c, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()
