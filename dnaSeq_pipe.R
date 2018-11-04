#Control pipeline for variant analysis of WGS data 
#
#Stephen Bolaris - UCR - GGB Program
#

#####Library imports
source("~/dnaSeq_scripts/dna_subroutines.R")
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(Biostrings)
library(AnnotationDbi)
library(ShortRead)
library(VariantAnnotation)
library(rtracklayer)
library(GenomicRanges)
library(VennDiagram)
library(ggplot2)
library(reshape2)

#use system args to get command line iput of target file, ref genome, and gff file / sqlite
arg = commandArgs(trailingOnly = TRUE)
#targets
targets = read.delim(arg[1],sep="\t", comment.char = "#")
ref = arg[2]
gff = arg[3]


if(!file.exists(paste(ref,".bwt", sep=""))){
	cmd = paste("bwa index -a bwstw ", ref, sep="")
	system(cmd)
}

#check for ref.fa -> ref.dic picard tools
#if(!file.exists(paste( 
#java -jar /opt/picard/1.112/CreateSequenceDictionary.jar REFERENCE=reference.fa OUTPUT=reference.dict

#check for samtools faidx
#samtools faidx reference.fa
if(!file.exists(paste(ref,".fai",sep=""))){
	cmd = paste("samtools faidx ", ref, sep="")
	system(cmd)
}

#makedir("./results/")

#loop through the tagets file
for(i in 1:nrow(targets)){
	sample_name = as.character(targets[i,]$SampleName)
	fastq1 = as.character(targets[i,]$FileName1)
	fastq2 = as.character(targets[i,]$FileName2)
	long_name = as.character(targets[i,]$SampleLong)
	#check if VCF files exist already if so skip
	if(!file.exists(paste("./results/",long_name,"_aln_pe.sam.gz",sep=""))){
		#call to bwa function(fastq1, fastq1, ref, sample_name, factor)
		bwa_alignment(fastq1, fastq2, ref, sample_name, as.character(targets[i,]$Factor))
	}
	if(!file.exists(paste("./results/",long_name,"_aln_pe_mk.bai",sep=""))){
		#call to picardtools(sample_name_long)
		picard_tools(long_name)
	}
	if(!file.exists(paste("./results/",long_name,"_filt_gindels.vcf",sep=""))){
		#call to GATK(sample_name_long)
		gatk(long_name,ref)
	}
		#call to Varscan(sample_name_long)
	if(!file.exists(paste("./results/",long_name,"_pileup_indels_filtered.vcf",sep=""))){
	#	varscan(long_name,ref)
	}
	if(!file.exists(paste("./results/",long_name,"_raw_gindels_filteRed.vcf",sep=""))){	
		#call to filteR(sample_name_long)
		#filteR(long_name,ref)
	}
	
}		
#end loop

#txdb = makeTxDbFromGFF(file=gff, format="gff3", dataSource="TAIR", organism="A. thaliana")
#core_db = renameSeqlevels(txdb, seqlevels(BamFile(paste("./results/",targets[1,]$SampleLong,"_re_alignment.bam", sep="")))

#make dir for report pdfs
dir.create("./report_pdfs/")
#vcf diagnostics
#reload targets file


#make a VCF list of each type (GATK)
snps_list = lapply(targets$SampleLong, function(x) readVcf(paste("./results/",x,"_filt_gsnps.vcf",sep=""),genome=ref))
names(snps_list) = targets$SampleLong
indel_list = lapply(targets$SampleLong, function(x) readVcf(paste("./results/",x,"_filt_gindels.vcf",sep=""), genome=ref))
names(indel_list) = targets$SampleLong
#pull out uniques of each for snpeff analysis
	#compare all vcf files for each line to pull out line snps
#get the filtered list snp_list$long_name[filt(snp_list$long_name) == PASS] can this be done with apply
filt_snp_list = lapply(snps_list, function(x) x[filt(x) == "PASS"])
filt_indel_list = lapply(indel_list, function(x) x[filt(x) == "PASS"])
##test for unfiltered Rate
#filt_snp_list = snps_list
#filt_indel_list = indel_list

names(filt_snp_list) = targets$SampleLong
names(filt_indel_list) = targets$SampleLong
rn_filt_snps = lapply(filt_snp_list, function(x) names(rowRanges(x)))
t_snp_freq = table(unlist(rn_filt_snps))
zero_snps = lapply(filt_snp_list[names(filt_snp_list) %in% targets[targets$Factor== "A0.0",]$SampleLong], function(x) names(rowRanges(x))) 
#check to see the that the line snps are at mendelian rates of at least 50% of the sampled population
line_snps = t_snp_freq[t_snp_freq >= (length(targets$SampleLong)*.4) & names(t_snp_freq) %in% unlist(zero_snps)]
#line_snps = unique(unlist(rn_filt_snps[targets[targets$Factor == "A0.0",]$SampleLong]))
rn_filt_indels = lapply(filt_indel_list, function(x) names(rowRanges(x)))
zero_indels = lapply(filt_indel_list[names(filt_indel_list) %in% targets[targets$Factor== "A0.0",]$SampleLong], function(x) names(rowRanges(x)))

t_indel_freq = table(unlist(rn_filt_indels))
#This could also be changed to 50% 
line_indels = t_indel_freq[t_indel_freq >= (length(targets$SampleLong)*.4) & names(t_indel_freq) %in% unlist(zero_indels)]
#line_indels = unique(unlist(rn_filt_indels[targets[targets$Factor == "A0.0",]$SampleLong]))
#line_snps = line_compare(line_snps, vcf3, "i")
u_filt_snp_list = lapply(1:length(filt_snp_list), function(x) setdiff(names(rowRanges(filt_snp_list[[x]])), names(line_snps)))
u_filt_indel_list = lapply(1:length(filt_indel_list), function(x) setdiff(names(rowRanges(filt_indel_list[[x]])), names(line_indels)))

names(u_filt_indel_list) = targets$SampleLong
names(u_filt_snp_list) = targets$SampleLong

#determine possible hotspots (variants that while not unique do not follow mendelian genetics)
hs = table(unlist(u_filt_indel_list[targets$SampleLong]))
hs0 = table(unlist(u_filt_indel_list[names(u_filt_indel_list) %in% targets[targets$Factor== "A0.0",]$SampleLong]))
hsA = table(unlist(u_filt_indel_list[names(filt_indel_list) %in% targets[targets$Factor!= "A0.0",]$SampleLong]))
hs = hs[hs>1]
hs0 = hs0[hs0 >1]
hsA = hsA[hsA >1]

names(hs) = targets$sampleLong
names(hs0) = targets[targets$Factor== "A0.0",]$SampleLong
names(hsA) = targets[targets$Factor!= "A0.0",]$SampleLong

#pull out unique variants and hotspots for further analysis
#unique / hotspot indels 
hs_indel_list = lapply(1:length(u_filt_indel_list), function(x) intersect(u_filt_indel_list[[x]], names(hs)))
unique_indel_list = lapply(1:length(u_filt_indel_list), function(x) setdiff(u_filt_indel_list[[x]], names(hs)))

#unique / hotspot snps
hss = table(unlist(u_filt_snp_list[targets$SampleLong]))
hss = hss[hss >1]

#sperate out hotspots (hs) from unique variants for further analysis
unique_snp_list = lapply(1:length(u_filt_snp_list), function(x) setdiff(u_filt_snp_list[[x]], names(hss)))
hs_snp_list = lapply(1:length(u_filt_snp_list), function(x) intersect(u_filt_snp_list[[x]], names(hss)))
names(unique_snp_list) = targets$SampleLong
names(unique_indel_list) = targets$SampleLong
names(hs_snp_list) = targets$SampleLong
names(hs_indel_list) = targets$SampleLong



#put variants in vcf form for further processing
uf_snp_vcfs = lapply(1:length(filt_snp_list), function(x) filt_snp_list[[x]][unique_snp_list[[x]]])
uf_indel_vcfs = lapply(1:length(filt_indel_list), function(x) filt_indel_list[[x]][unique_indel_list[[x]]])
hs_snp_vcfs = lapply(1:length(filt_snp_list), function(x) filt_snp_list[[x]][hs_snp_list[[x]]])
hs_indel_vcfs = lapply(1:length(filt_indel_list), function(x) filt_indel_list[[x]][hs_indel_list[[x]]])

names(uf_snp_vcfs) = targets$SampleLong
names(uf_indel_vcfs) = targets$SampleLong
names(hs_snp_vcfs) = targets$SampleLong
names(hs_indel_vcfs) = targets$SampleLong

#to check 2bp AT deletions,AT_del = lapply(uf_indel_vcfs, function(x) names(rowRanges(x))[grep('*._[ACTG](AT|TA)/[ACTG]', names(rowRanges(x)))])
#AT_vcfs = lapply(1:length(AT_del), function(x) uf_indel_vcfs[x][AT_del[x]])
#names(AT_vcfs) = targets$SampleLong
#lapply(AT_vcfs, function(x) homology_check(x, ref))

#break down variants for pattern processing.
unique_indel_breakdown = lapply(uf_indel_vcfs, function(x) indel_breakdown(x))
indel_stats = indel_metrics(unique_indel_breakdown)
write.table(indel_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_boxplot.pdf",sep=""),width=15);box_graph(indel_stats, s=1);dev.off()
### do indel stats with F-test and Poisson

dfm = t(indel_stats)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"s")))
sink()

indel_stats_small = indel_stats[,1:4]
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_1to4_boxplot.pdf",sep=""),width=15);box_graph(indel_stats_small, s=1);dev.off()

dfm = t(indel_stats_small)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_numeric_1to4_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"s")))
sink()



unique_snp_breakdown = lapply(uf_snp_vcfs, function(x) indel_breakdown(x))
snp_stats = breakdown_table(unique_snp_breakdown)
write.table(snp_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_boxplot.pdf",sep=""));box_graph(snp_stats);dev.off()
p_snp = pois_stat(snp_stats)
p_snp_i = pois_stat(snp_stats, "i")
f_snp = f_test(snp_stats, "i")
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_snp))
print("Poisson, no interaction\n")
print(summary(p_snp))
print("Poisson, with interaction \n")
print(summary(p_snp_i))
sink()

#repeat the process for the hotspots
hotspot_indel_breakdown = lapply(hs_indel_vcfs, function(x) indel_breakdown(x))
hs_indel_stats = indel_metrics(hotspot_indel_breakdown)
write.table(hs_indel_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_hotspot_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_hs_boxplot.pdf",sep=""),width=15);box_graph(hs_indel_stats, s=1);dev.off()

dfm = t(hs_indel_stats)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_hs_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"s")))
sink()

hs_indel_stats_small = hs_indel_stats[,1:4]
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_hs_1to4_boxplot.pdf",sep=""),width=15);box_graph(hs_indel_stats_small, s=1);dev.off()

dfm = t(hs_indel_stats_small)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_hs_numeric_1to4_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"s")))
sink()




hs_snp_breakdown = lapply(hs_snp_vcfs, function(x) indel_breakdown(x))
hs_snp_stats = breakdown_table(hs_snp_breakdown)
write.table(hs_snp_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_hotspot_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_hs_boxplot.pdf",sep=""));box_graph(hs_snp_stats);dev.off()
p_snp = pois_stat(hs_snp_stats)
f_snp = f_test(hs_snp_stats,"i")
p_snp_i = pois_stat(snp_stats, "i")
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_snp_hs_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_snp))
print("Poisson, no interaction\n")
print(summary(p_snp))
print("Poisson, with interaction \n")
print(summary(p_snp_i))
sink()


#seperate out insertions and deletions
unique_indel_breakdown_i = lapply(unique_indel_breakdown, function(x) x[x$type == "insert",])
unique_indel_breakdown_d = lapply(unique_indel_breakdown, function(x) x[x$type == "deletion",])
insertion_stats = indel_metrics(unique_indel_breakdown_i)
deletion_stats = indel_metrics(unique_indel_breakdown_d)
write.table(deletion_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_deletion_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(insertion_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_insertion_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_insert_boxplot.pdf",sep=""),width=15);box_graph(insertion_stats, s=1);dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_indel_deletion_boxplot.pdf",sep=""), width=15);box_graph(deletion_stats, s=1);dev.off()


dfm = t(insertion_stats)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_insertion_un_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"i")))
sink()
dfm = t(deletion_stats)
sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_deletion_un_numeric_stats.txt",sep=""))
print("ANOVA\n")
print(summary(f_test(dfm,"i")))
print("Poisson, no interaction\n")
print(summary(pois_stat(dfm)))
print("Poisson, with interaction \n")
print(summary(pois_stat(dfm,"i")))
sink()


#focus the breakdown on specific basepair sizes
bp1 = breakdown_table(unique_indel_breakdown, 1)
bp2 = breakdown_table(unique_indel_breakdown, 2)
bp3 = breakdown_table(unique_indel_breakdown, 3)
bp4 = breakdown_table(unique_indel_breakdown, 4)
box1 = box_graph(bp1) 
box2 = box_graph(bp2)
box3 = box_graph(bp3)
box4 = box_graph(bp4)
#table output of the counts
write.table(bp1, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_1bpVar_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp2, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_2bpVar_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp3, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_3bpVar_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp4, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_4bpVar_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
#pdfs of box plots
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_1bpVar_numeric_boxplot.pdf",sep=""),width=15)
box1
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_2bpVar_numeric_boxplot.pdf",sep=""),width=15)
box2
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_3bpVar_numeric_boxplot.pdf",sep=""),width=15)
box3
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_4bpVar_numeric_boxplot.pdf",sep=""),width=15)
box4
dev.off()
#Stat out put in text format
p1 = pois_stat(bp1)
p2 = pois_stat(bp2)
p3 = pois_stat(bp3)
p4 = pois_stat(bp4)
p1_i = pois_stat(bp1,"s")
p2_i = pois_stat(bp2,"s")
p3_i = pois_stat(bp3,"s")
p4_i = pois_stat(bp4,"s")
f1 = f_test(bp1, "i")
f2 = f_test(bp2, "i")
f3 = f_test(bp3, "i")
f4 = f_test(bp4, "i")


sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_poisson_stats.txt",sep=""))
print("ANOVA Testing\n")
print(summary(f1))
print(summary(f2))
print(summary(f3))
print(summary(f4))
print("Poisson No interaction\n")
print(summary(p1))
print(summary(p2))
print(summary(p3))
print(summary(p4))
print("Poisson with interaction\n")
print(summary(p1_i))
print(summary(p2_i))
print(summary(p3_i))
print(summary(p4_i))
sink()

#hotspot analytics based on bp size
bp1 = breakdown_table(hotspot_indel_breakdown, 1)
bp2 = breakdown_table(hotspot_indel_breakdown, 2)
bp3 = breakdown_table(hotspot_indel_breakdown, 3)
bp4 = breakdown_table(hotspot_indel_breakdown, 4)
box1 = box_graph(bp1)
box2 = box_graph(bp2)
box3 = box_graph(bp3)
box4 = box_graph(bp4)
write.table(bp1, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_1bpVar_hs_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp2, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_2bpVar_hs_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp3, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_3bpVar_hs_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
write.table(bp4, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_4bpVar_hs_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)

pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_1bpVar_hs_numeric_boxplot.pdf",sep=""),width=15)
box1
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_2bpVar_hs_numeric_boxplot.pdf",sep=""),width=15)
box2
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_3bpVar_hs_numeric_boxplot.pdf",sep=""),width=15)
box3
dev.off()
pdf(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_4bpVar_hs_numeric_boxplot.pdf",sep=""),width=15)
box4
dev.off()

#Stat out put in text format
p1 = pois_stat(bp1)
p2 = pois_stat(bp2)
p3 = pois_stat(bp3)
p4 = pois_stat(bp4)
p1_i = pois_stat(bp1,"s")
p2_i = pois_stat(bp2,"s")
p3_i = pois_stat(bp3,"s")
p4_i = pois_stat(bp4,"s")
f1 = f_test(bp1, "i")
f2 = f_test(bp2, "i")
f3 = f_test(bp3, "i")
f4 = f_test(bp4, "i")

sink(paste("./results/",gsub("R.*","",targets$SampleName[1]),"_hs_poisson_stats.txt",sep=""))
print("ANOVA Testing\n")
print(summary(f1))
print(summary(f2))
print(summary(f3))
print(summary(f4))
print("Poisson No interaction\n")
print(summary(p1))
print(summary(p2))
print(summary(p3))
print(summary(p4))
print("Poisson with interaction\n")
print(summary(p1_i))
print(summary(p2_i))
print(summary(p3_i))
print(summary(p4_i))
sink()

### Plotting boxplot with scatter plot #### ---> rolled in to box_graph subroutine
#dfm = melt(bp1)
#dfm$variable = gsub(".*_A", "", dfm$variable)
#dfm$variable = as.factor(dfm$variable)  
#names(dfm) = c("Change", "Al_mM", "Counts")
#pdf(sample_boxplot_1bp.pdf)
#ggplot(dfm aes(Change,Counts) ) + geom_boxplot(color=Al_mM)+theme(axis.text.x=element_text(angle=90,hjust=1))
#dev.off()

#### computing statistical significance ####
# f-test
#res = aov(Counts~Change*Al_mM,dfm)
#summary(res)
#TukeyHSD(res, which="Al_mM")
#
# Paired t-test
#t_test_pairs = lapply(unique(dfm$Change), function(x) pairwise.t.test(dfm[dfm$Change == x,]$Counts,dfm[dfm$Change ==x,]$Al_mM,p.adjust.method="bonf"))
#names(t_test_pairs) = unique(dfm$Change)


uf_snp_vcfs = lapply(targets$SampleLong, function(x) renameSeqlevels(uf_snp_vcfs[[x]],c("1","2","3","4","5","M","C")))
uf_indel_vcfs = lapply(targets$SampleLong, function(x) renameSeqlevels(uf_indel_vcfs[[x]],c("1","2","3","4","5","M","C")))
hs_snp_vcfs = lapply(targets$SampleLong, function(x) renameSeqlevels(hs_snp_vcfs[[x]],c("1","2","3","4","5","M","C")))
hs_indel_vcfs = lapply(targets$SampleLong, function(x) renameSeqlevels(hs_indel_vcfs[[x]],c("1","2","3","4","5","M","C")))

names(uf_snp_vcfs) = targets$SampleLong
names(uf_indel_vcfs) = targets$SampleLong
#write out vcf of unique and filtered snps
lapply(targets$SampleLong, function(x) writeVcf(uf_snp_vcfs[[x]], paste("./results/",x,"_unique_snps.vcf",sep="")))
#write vcf of unique and filtered indels
lapply(targets$SampleLong, function(x) writeVcf(uf_indel_vcfs[[x]], paste("./results/",x,"_unique_indels.vcf",sep="")))
#write out vcf of unique and filtered snps
lapply(targets$SampleLong, function(x) writeVcf(hs_snp_vcfs[[x]], paste("./results/",x,"_hs_snps.vcf",sep="")))
#write vcf of unique and filtered indels
lapply(targets$SampleLong, function(x) writeVcf(hs_indel_vcfs[[x]], paste("./results/",x,"_hs_indels.vcf",sep="")))


#use breakdown_table(breakdown, size<-optional)
#unique_breakdown = lapply(unique_indel_vcfs[targets[targets$Factor != "A0.0",]$SampleLong], function(x) indel_breakdown(x))
#(also called uf_breakdown further down....
#spec_indel_vcf = lapply(unique_indel_vcfs, function(x) x[rowRanges(x)$REF %in% c("AAT","CAT", "GAT", "CTA", "GTA", "TTA")])
#h_check = lapply(spec_indel_vcf[targets[targets$Factor != "A0.0",]$SampleLong], function(x) homology_check(x, ref))

tot_snp = lapply(unique_snp_list, function(x) length(x))
tot_indel = lapply(unique_indel_list, function(x) length(x)) 
names(tot_snp) = names(filt_snp_list)
names(tot_indel) = names(filt_indel_list)
m_snp = melt(tot_snp)
m_indel = melt(tot_indel)

#barplot comparing filtered snps of each line
pdf(paste("./report_pdfs/",targets$SampleName,"_Total_snp_compare.pdf",sep=""))
ggplot(m_snp,aes(x=L1,y=value,fill=L1))+geom_bar(stat="identity",position="dodge")+xlab("Sample")+ylab("SNP Frequency")+ scale_fill_brewer()+theme_dark()+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))+labs(title="Total Filt. SNPs Per Sample")
dev.off()
pdf(paste("./report_pdfs/",targets$SampleName,"_Total_indel_compare.pdf",sep=""))
ggplot(m_indel,aes(x=L1,y=value,fill=L1))+geom_bar(stat="identity",position="dodge")+xlab("Sample")+ylab("INDEL Frequency")+scale_fill_brewer()+theme_dark()+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))+labs(title = "Total Filt. INDELs per Sample")
dev.off()

#print("SNP EFF")

for(i in targets$SampleLong){
	if(!file.exists(paste("./results/",i,"_Uindels_ann.vcf",sep=""))){
		snpeff(i)
	}
}

#print("SNP EFF complete")
#this could be replaced by the Db.sqlite
#core_db =  renameSeqlevels(txdb, seqlevels(bamfl[[1]]))
#exonXgene = exonsBy(core_db, by="gene")
#bamfiles = BamFileList(bamfl, yeildSize=50000, index=character())
#countEbyGene = lapply(bamfl, function(x) summarizeOverlaps(exonXgene, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=FALSE))
gff = import.gff3(gff)
gff = renameSeqlevels(gff,c("1","2","3","4","5","M","C"))
gff_genes = gff[which(elementMetadata(gff)[,"type"] == "gene")]
gff_exons = gff[which(elementMetadata(gff)[,"type"] == "exon")]

for(i in targets$SampleLong){ 
	
	curr_vcf = uf_snp_vcfs[[i]]
	curr_vcf
	
	by_chrom = num_chrmChng(curr_vcf)
	pdf(paste("report_pdfs/",i,"_A_Chgs_by_chrom.pdf",sep=""))
	barplot(by_chrom, xlab="Chromosome", ylab="Frequency", main="SNP Freq. by Chromosome*", col=c("black", "yellow", "purple","blue", "cyan", "orange","gold"))
	dev.off()
	#nucleotide changes
	nc = nuc_changes(curr_vcf)
	pdf(paste("report_pdfs/",i,"_B_Ref_base_changes.pdf",sep=""))
	barplot(nc$A, xlab= "Alt Allele", ylab = "frequency", main="Ref A Allele", col=c("yellow", "blue", "purple"))
	barplot(nc$C, xlab= "Alt Allele", ylab = "frequency", main="Ref C Allele", col=c("black", "blue", "purple"))
	barplot(nc$G, xlab= "Alt Allele", ylab = "frequency", main="Ref G Allele", col=c("black","yellow", "purple"))
	barplot(nc$T, xlab= "Alt Allele", ylab = "frequency", main="Ref T Allele", col=c("black", "yellow", "blue"))
	dev.off()
	#genomic region dispersal
	total = dim(curr_vcf)[1]
	gol = subsetByOverlaps(rowRanges(curr_vcf), gff_genes)
	eol = subsetByOverlaps(gol, gff_exons)
	intergenic = total - length(gol)
	genes = length(gol)
	exons = length(eol)
	intron = genes - exons
	type = c("Exons","Introns/UTR","Intergenic")
	location = c("Gene","Gene","Intergenic")
	breakdownS = c(exons, intron, intergenic)
	df = data.frame(location,type, breakdownS)
	snp_rd = ggplot(df,aes(x=location, y=breakdownS, fill=type))+geom_bar(stat="identity")+labs(title="GATK Fltr. SNP Genomic Region Dispersal", x="Genomic Region", y="Frequency")+scale_fill_brewer(palette="YlOrBr")+theme_dark()	

	ggsave(paste("report_pdfs/",i,"_C_SNP_region_dispersal.pdf",sep=""), snp_rd)	
	#predict how many would lead to mutations
	#read in vcf file
	sE = readVcf(paste("./results/",i,"_Usnps_ann.vcf",sep=""),genome=ref)
	#compare the number that hit potein coding regions and the type of change

	####Write Mertric files for each PDF
	
}

insert_content = data.frame()
#row(insert_content) = c("AT","GC")
#colnames(insert_content) = c("sample", "AT","GC")   
deletion_content = data.frame()
#colnames(deletion_content) = c("sample", "AT","GC")   

for(i in targets$SampleLong){ 
	curr_vcf = uf_indel_vcfs[[i]]
	#size distrubtion / con
	size_dist = indel_breakdown(curr_vcf)
	size_dist_plot = ggplot(size_dist, aes(size,fill=type))+geom_bar(position="dodge")+scale_x_discrete(limits=unique(size_dist$size))+labs(title ="GATK filt. INDEL Size Distribution")+scale_fill_brewer(palette="PuOr")+theme_dark()

	ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D_indel_size_dist.pdf",sep=""))
	#add bottom title to specify 2-3 bp	
	size_dist_plot = ggplot(subset(size_dist,size <100), aes(size,fill=type))+theme_dark()+geom_bar(position="dodge")+scale_x_discrete(limits=unique(size_dist$size))+labs(title =paste(i," INDEL Size Distribution",sep=""), x="INDEL size (bp)", y="INDEL Frequency")+scale_fill_manual(values=c("light green","tan"))+ scale_y_continuous(limits = c(0, 275))+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))

        ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D2_indel_size_dist_100.pdf",sep=""))

 ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D_indel_size_dist.pdf",sep=""))
        #add bottom title to specify 2-3 bp     
        size_dist_plot = ggplot(subset(size_dist,size <50), aes(size,fill=type))+theme_dark()+geom_bar(position="dodge")+scale_x_discrete(limits=unique(size_dist$size))+labs(title =paste(i," INDEL Size Distribution",sep=""), x="INDEL size (bp)", y="INDEL Frequency")+scale_fill_manual(values=c("light green","tan"))+ scale_y_continuous(limits = c(0, 275))+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))

        ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D2_indel_size_dist_50.pdf",sep=""))
 ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D_indel_size_dist.pdf",sep=""))
        #add bottom title to specify 2-3 bp     
        size_dist_plot = ggplot(subset(size_dist,size <25), aes(size,fill=type))+theme_dark()+geom_bar(position="dodge")+scale_x_discrete(limits=unique(size_dist$size))+labs(title =paste(i," INDEL Size Distribution",sep=""), x="INDEL size (bp)", y="INDEL Frequency")+scale_fill_manual(values=c("light green","tan"))+ scale_y_continuous(limits = c(0, 275))+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))

        ggsave(size_dist_plot, file=paste("report_pdfs/",i,"_D2_indel_size_dist_25.pdf",sep=""))


	
	#frequency by chromosome
	by_chrom = num_chrmChng(curr_vcf)
	pdf(paste("report_pdfs/",i,"_E_Chgs_by_chrom.pdf",sep=""))
        barplot(by_chrom, xlab="Chromosome", ylab="Frequency", main="INDEL Freq. by Chromosome*", col=c("black", "yellow", "purple","blue", "cyan", "orange","gold"))
        dev.off()
	#base content
	bc = base_content(curr_vcf)
	
	data = bc[,c("type","AT","GC")]
	data.m = melt(data)
	g =ggplot(data.m[data.m$value>0,], aes(x=type, y=value,fill=variable))+geom_bar(stat="identity")+labs(title="Indel Base Content", x="INDEL type", y="Frequency")+scale_fill_brewer(palette="BrBG")+theme_dark()
	ggsave(paste("report_pdfs/",i,"_F_Indel_base_content.pdf",sep=""),g)
	
	insert_sub = subset(bc, type=="insert")
	itot = sum(insert_sub$AT)+sum(insert_sub$GC)
	insert_content = rbind(insert_content, name =  c((sum(insert_sub$AT)/itot*100),(sum(insert_sub$GC)/itot*100)))
	
	del_sub = subset(bc, type=="deletion")
	dtot = sum(del_sub$AT)+sum(del_sub$GC)
	deletion_content = rbind(deletion_content,name = c((sum(del_sub$AT)/dtot*100),(sum(del_sub$GC)/dtot*100)))
	
	#print(insert_content)
	#genomeic region dispersal
	total = dim(curr_vcf)[1]
        gol = subsetByOverlaps(rowRanges(curr_vcf), gff_genes)
        eol = subsetByOverlaps(gol, gff_exons)
        intergenic = total - length(gol)
        genes = length(gol)
        exons = length(eol)
        intron = genes - exons
        type = c("Exons","Introns/UTR","Intergenic")
        location = c("Gene","Gene","Intergenic")
        breakdownS = c(exons, intron, intergenic)
        df = data.frame(location,type, breakdownS)
        indel_rd = ggplot(df,aes(x=location, y=breakdownS, fill=type))+geom_bar(stat="identity")+labs(title="GATK Fltr. INDEL Genomic Region Dispersal", x="Genomic Region", y="Frequency")+scale_fill_brewer(palette="YlOrBr")+theme_dark()
        ggsave(paste("report_pdfs/",i,"_G_indel_region_dispersal.pdf",sep=""), indel_rd) 
	#predict how many would lead to mutations
	readVcf(paste("./results/",i,"_Uindels_ann.vcf",sep=""), genome=ref)

	####Write metric files for each PDF to have text versions
}
#Indel metrics for totals -> export this to the subroutines.... 
print("indel metrics")
#uf_breakdown = lapply(uf_indel_vcfs, function(x) indel_breakdown(x))
uf_breakdown = lapply(uf_indel_vcfs, function(x) indel_breakdown(x))
names(uf_breakdown) = names(targets$SampleLong)
#write out table showing distribution of indels in table format
indel_stats = indel_metrics(uf_breakdown)
#write out table
indel_stats
fn = paste("./results/",gsub("_.*","",targets$SampleName[1]),"_total_indel_numeric_metrics.xls",sep="")
print(fn)
write.table(indel_stats,fn, sep="\t", col.names=NA)

uf_breakdown_d = lapply(uf_breakdown, function(x) x[x$type == "deletion",])
names(uf_breakdown_d) = names(uf_breakdown)
uf_breakdown_i = lapply(uf_breakdown, function(x) x[x$type == "insert",])
names(uf_breakdown_i) = names(uf_breakdown)

deletion_stats = indel_metrics(uf_breakdown_d)

write.table(deletion_stats, paste("./results/",gsub("R.*","",targets$SampleName[1]),"_total_deletion_numeric_metrics.xls",sep=""), sep="\t", col.names=NA,append=FALSE)
print("table 1 ready")
#write out deletion table
insertion_stats = indel_metrics(uf_breakdown_i)

write.table(insertion_stats, paste("./results/",gsub("_.*","",targets$SampleName[1]),"_total_insertion_numeric_metrics.xls",sep=""), sep="\t", col.names=NA)
#write out insertion table


#PDF for insertion and deletion base content amoung the three samples


colnames(insert_content) = c("AT","GC")
colnames(deletion_content) = c("AT","GC")
inc = melt(insert_content)
inc["sample"] <- rep(targets$SampleLong,2)
dec = melt(deletion_content)
dec["sample"] <- rep(targets$SampleLong,2)
pdf("report_pdfs/sample_insertion_content_comp.pdf")
ggplot(inc, aes(x=sample, y=value,fill=variable))+geom_bar(stat="identity",position="dodge")+labs(title="Insertion Base Content", x="Sample", y="Content %")+scale_fill_brewer(palette="BrBG")+theme_dark()+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
pdf("report_pdfs/sample_deletion_content_comp.pdf")
ggplot(dec, aes(x=sample, y=value,fill=variable))+geom_bar(stat="identity",position="dodge")+labs(title="Deletion Base Content", x="Sample", y="Content %")+scale_fill_brewer(palette="BrBG")+theme_dark()+ theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
#for indels confirm that the T-DNA is present as an indel using the primer sequence
#LBa1  5’ tggttcacgtagtgggccatcg 3’


#gnerate PDF report
##pdfunite report_pdfs/*.pdf $dir-results.pdf
for(s in targets$SampleLong){
	cmd = paste("pdfunite ./report_pdfs/",s,"_*.pdf ./results/",s,"_results.pdf",sep="")
	system(cmd)
}

	

sessionInfo()

