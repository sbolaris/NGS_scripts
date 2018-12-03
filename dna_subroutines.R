#function to creat alignments with bwa
bwa_alignment = function(fq1, fq2, ref, sample_name, exfactor){
	#changed \" to ' in paste to see if compatability issue with @PG is fixed
	header = paste("'@RG\tID:",sample_name,"\tSM:",exfactor,"\tPL:illumina\tLB:lib1\tPU:unit1'",sep="")
	cmd = paste("bwa mem -aM -t 8 -R ", header," ",ref," ",fq1," ",fq2," ",">./results/",sample_name,"_",exfactor,"_aln_pe.sam",sep="")
	system(cmd)
	#picard has issues with the formating of the PG line that BWA puts in when doing the alignments
	cmd = paste("sed -i '/^@PG/d' ./results/",sample_name,"_",exfactor,"_aln_pe.sam",sep="")
	system(cmd)        
}

#function for QC with picard tools
picard_tools = function (sample_name){
	cmd = paste("picard-tools SortSam INPUT=./results/",sample_name,"_aln_pe.sam OUTPUT=./results/",sample_name,"_aln_pe.bam SO=coordinate",sep="")
	system(cmd)
	cmd = paste("picard-tools MarkDuplicates INPUT=./results/",sample_name,"_aln_pe.bam OUTPUT=./results/",sample_name,"_aln_pe_mk.bam METRICS_FILE=./results/",sample_name,"_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true",sep="")
	system(cmd)
	cmd = paste("picard-tools BuildBamIndex INPUT=./results/",sample_name,"_aln_pe_mk.bam",sep="")
	system(cmd)
	file = paste("./results/",sample_name,"_aln_pe.sam",sep="")
	gzip(file)
	file = paste("./results/",sample_name,"_aln_pe.bam",sep="")
	gzip(file)
}

gzip = function(file_name){
	cmd = paste("gzip ",file_name, sep="")
	system(cmd)
}

#function to call variants with GATK
gatk = function(sample_name, ref){
	gatk = "java -jar ~/GenomeAnalysisTK.jar"
	cmd = paste(gatk," -T RealignerTargetCreator -R ",ref," -I ./results/",sample_name,"_aln_pe_mk.bam -o ./results/",sample_name,"_indels.intervals",sep="")
	system(cmd)
	cmd = paste(gatk," -T IndelRealigner -R ",ref," -I ./results/",sample_name,"_aln_pe_mk.bam -targetIntervals ./results/",sample_name,"_indels.intervals -o ./results/",sample_name,"_re_alignment.bam --filter_bases_not_stored",sep="")
	system(cmd)
	cmd = paste(gatk," -T HaplotypeCaller -R ",ref," -I ./results/",sample_name,"_re_alignment.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ./results/",sample_name,"_raw_gvar.vcf", sep="")
	system(cmd)
	cmd = paste(gatk," -T SelectVariants -R ",ref," -V ./results/",sample_name,"_raw_gvar.vcf -selectType SNP -o ./results/",sample_name,"_raw_gsnps.vcf", sep="")
	system(cmd)
	fexp = "\"QD < 2.0 || FS > 60.0 || MQ < 40.0 ||HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 ||ReadPosRankSum < -8.0\""
	fname = "\"snp_filter\""
	cmd = paste(gatk," -T VariantFiltration -R ",ref," -V ./results/",sample_name,"_raw_gsnps.vcf --filterExpression ",fexp," --filterName ",fname," -o ./results/",sample_name,"_filt_gsnps.vcf", sep="")
	system(cmd)
	cmd = paste(gatk," -T SelectVariants -R ",ref," -V ./results/",sample_name,"_raw_gvar.vcf -selectType INDEL -o ./results/",sample_name,"_raw_gindels.vcf",sep="")
	system(cmd)
	fexp = "\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\""
	fname = "\"my_indel_filter\""
	cmd = paste(gatk," -T VariantFiltration -R ",ref," -V ./results/",sample_name,"_raw_gindels.vcf --filterExpression ",fexp," --filterName ",fname," -o ./results/",sample_name,"_filt_gindels.vcf", sep="")
	system(cmd)
	file = paste("./results/",sample_name,"_re_alignment.bam",sep="")
        gzip(file)
}

#function to call variatns with varscan
varscan = function(sample_name, ref){
	varscan = "java -jar ~/VarScan.v2.3.9.jar"
	cmd = paste("samtools mpileup -A -f ",ref," ./results/",sample_name,"_re_alignment.bam >./results/",sample_name,".mpileup", sep="")
	system(cmd)
	cmd = paste(varscan," mpileup2snp ./results/",sample_name,".mpileup --output-vcf 1 >./results/",sample_name,"_pileup_snps.vcf", sep="")
	system(cmd)
	cmd = paste(varscan," mpileup2indel ./results/",sample_name,".mpileup --output-vcf 1 >./results/",sample_name,"_pileup_indels.vcf", sep="")
	system(cmd)
	cmd = paste(varscan," filter ./results/",sample_name,"_pileup_snps.vcf --min-read-depth 20 --min-var-feq 0.8 --output-file ./results/",sample_name,"_pileup_snps_filtr.vcf", sep="")
	system(cmd)
	cmd = paste(varscan," filter ./results/",sample_name,"_pileup_indels.vcf --min-read-depth 20 --min-var-freq 0.8 --output-file ./results/",sample_name,"_pileup_indels_filtered.vcf", sep="")
	system(cmd)
}

#function filter GATK variants using R (inbetween of GATK and Varscan)
filteR = function(sample_name, ref){
	library(VariantAnnotation)
	vcf_list = c(readVcf(paste("./results/",sample_name,"_raw_gsnps.vcf",sep=""), ref), readVcf(paste("./results/",sample_name,"_raw_gindels.vcf",sep=""),ref))
	names(vcf_list) = c(paste(sample_name,"_snps", sep=""), paste(sample_name, "_indels",sep=""))
	for(vcf in vcf_list){
	#read depth = geno(vcf)$DP
	# check if read_depth >= 20
		p1 = vcf[which(geno(vcf)$DP >=20)]
#alt allele frequence geno(vcf[i])$AD[[1]][2]
#check if alt_feq /read_depth >= 0.8

		p2 = sapply(1:dim(p1)[1], function(x) p1[which(geno(p1[x])$AD[[1]][2]/geno(p1[x])$DP[1] >= 0.8)])
		index = c()
		limit = dim(p1)[1]
		for( i in 1:limit){
        		altd = geno(p1[i])$AD[[1]][2]
#remove .vcf then add _filteRed.vcf
        		totd = geno(p1[i])$DP[1]
        		if(altd / totd >= 0.5){
                		index = c(index, i)
       			}
		}
		p2 = p1[index]
#need to figure out what the softfilter matrix in girke program does
		out_file = paste("./results/",sample_name,"_filteRed.vcf", sep="")
		writeVcf(p2, out_file)
	}

}
#function for using SNPeff to get the predicted impact use TAIR10.29 database
snpeff = function(sample_name){
	
	cmd = paste("java -Xmx4g -jar ~/Downloads/snpEff/snpEff.jar -csvStats ./results/",sample_name,"_snp_ann_stats.csv TAIR10.29 ./results/",sample_name,"_unique_snps.vcf >./results/",sample_name,"_Usnps_ann.vcf", sep="")
	system(cmd)
	
	##rename summary file and out put as csv with -csv option for further analysics later on	
	
	cmd = paste("java -Xmx4g -jar ~/Downloads/snpEff/snpEff.jar -csvStats ./results/",sample_name,"_indel_ann_stats.csv TAIR10.29 ./results/",sample_name,"_unique_indels.vcf >./results/",sample_name,"_Uindels_ann.vcf", sep="")
	system(cmd)
	
	##repeat renaming so that results do not get overwritten for indels
}


line_compare =function(v1, v2, t){
        line1_string = names(rowRanges(v1))
        line2_string = names(rowRanges(v2))
        res = c()
        if(t == "u1"){
        res = setdiff(line1_string, line2_string)
        }
        if(t == "i"){
        res = intersect(line1_string, line2_string)
        }
        if(t == "u2"){
        res = setdiff(line2_string, line1_string)
        }
        return(res)
}

num_chrmChng = function(v){
        chrom = c("1","2","3","4","5","mitochondria","chloroplast")
        num_changes = c()
        for( i in chrom){
                changes = v[seqnames(rowRanges(v)) == i]
                num_changes = c(num_changes, dim(changes)[1])
        }
        names(num_changes) = c("1","2","3","4","5","M","C")
        return(num_changes)
}

get_alt = function(vcf, ref_letter){
        nuc = c("A","C","G","T")
        alt = do.call(c, alt(vcf))
        al_let = c()
        freq = c()
        for ( i in nuc) {
                if(i == ref_letter){
                        next
                }
                freq = c(freq, length(alt[alt == DNAString(i)]))
                al_let = c(al_let, i)
        }

        names(freq) = al_let
        return(freq)
}

nuc_changes = function(vcf){
	ref_a = vcf[ref(vcf) == "A"]
        ref_t = vcf[ref(vcf) == "T"]
        ref_c = vcf[ref(vcf) == "C"]
        ref_g = vcf[ref(vcf) == "G"]
	nuc_chg = list(ref_a, ref_t, ref_c, ref_g)
	names(nuc_chg) = c("A","T", "C", "G")
	res = lapply(names(nuc_chg), function(x) get_alt(nuc_chg[[x]],x))
	names(res) = names(nuc_chg)
	return(res)
}

gene_overlaps = function(vcf, gff_sub){
        chrom_names = c("1", "2", "3", "4", "5", "chloroplast", "mitochondria")
        gff_chrom = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
        ol = c()
        for( i in 1:7){
                gff_region = gff_sub[seqnames(gff_sub) == gff_chrom[i]]
                vcf_region = vcf[seqnames(rowRanges(vcf)) == chrom_names[i]]
                #should I be using subsetbyoverlap instead of find overlaps
                ol = c(ol, findOverlaps(ranges(vcf_region),ranges(gff_region)))
        }
        return(ol)
}

#Function that takes a vcf file as input and breaks down the changes in to type, size and includes 
#the actual variant change
indel_breakdown = function(v){
        # get size of ref 
        #r = ref(v)
        #a = do.call(c, alt(v))
        n = names(rowRanges(v))
        d = c()
        i = c()
	dc = c()
	ic = c()
        changes = sub(".*_","",n)
        changes = strsplit(changes, "/")
        for ( x in changes){

                if(nchar(x[1]) > nchar(x[2])){
                        size = (nchar(x[1]) - nchar(x[2]))
                        d = c(d, size)
			dc = c(dc, paste(x[1],"*",x[2]))
                }
                else{
                        size = (nchar(x[2]) - nchar(x[1]))
                        i = c(i, size)
			ic = c(ic, paste(x[1],"*",x[2]))
                }
        }
        n = c(rep("insert", length(i)),rep("deletion", length(d)))
        id = c(i,d)
	ch = c(ic,dc)
        inde = data.frame(n, id, ch)
        names(inde) = c("type", "size", "change")
        return(inde)
}

#function for determing the base content of indels in a specific vcf files can be combined with 
#things like lappy to do multiple vcfs at once
base_content = function(v){
        n = names(rowRanges(v))
        d = c()
        i = c()
        changes = sub(".*_","",n)
        changes = strsplit(changes, "/")
        for ( x in changes){

                if(nchar(x[1]) > nchar(x[2])){
                        d = c(d, x[2])
                }
                else{
                        i = c(i, x[2])

                }
        }

        n = c(rep("insert", length(i)),rep("deletion", length(d)))
        id = c(i,d)
        dat =c()
        dcg =c()
        for( y in id){
                dat = c(dat,(count_occ(y,"A")+count_occ(y,"T")))
                dcg = c(dcg,(count_occ(y,"C")+count_occ(y,"G")))
        }
        inde = data.frame(n, id,dat,dcg)
        names(inde) = c("type", "change", "AT","GC")
        return(inde)

}
#counts the number of occurances for a character in a specific string
count_occ = function(string, char){
        string2 = gsub(char,"",string)
        return(nchar(string)-nchar(string2))
}
#function for creating a list of line changes common to multiple vcf files and returns a list of 
#changes that they have in common for filtering purposes
line_changes = function(change_list){
	list_changes = list()
	for(i in 1:length(change_list)){
		int = c()
		for(j in 1:length(change_list)){
			#print(i,j)
			if(i == j){
				next
			}
			else{
				int = c(int, intersect(change_list[[i]],change_list[[j]]))
				#print(int)
			}
		}
		changes = unique(int)
	list_changes[[length(list_changes)+1]] = changes
	}
	return(list_changes)

}

#allows the user to select from a list of brokendown indels a specific type and 
indel_select = function( indel_list, t, s){
		res = lapply(indel_list, function(x) x[(x$type == t &x$size == s),])
	
	return(res)
}

#Takes the unique filtered breakdown of indels and creates a numerical table to for export
indel_metrics = function(breakdown){
	t_list = lapply(breakdown, function(x) table(x$size))
	names = lapply(t_list, function(x) names(x))
	uNames = sort(as.integer(unique(unlist(names))))
	#write out table showing distribution of indels in table format
	m = matrix(0, length(t_list), length(uNames))
	rownames(m) = names(t_list)
	colnames(m) = uNames
	for(i in 1:length(t_list)){
        	m[i,intersect(names(t_list[[i]]), colnames(m))] = t_list[[i]]

	}

	return(m)
}


breakdown_table = function(bp_bd, s=0) {
        #subroutine takes the values from the breakdown values of the variants
        if(s == 0){
        #check if default value has changes
                bd_change = lapply(bp_bd, function(x) x$change)
        }
        else{
        #given a specific size only perform table on those items
                bd_change = lapply(bp_bd, function(x) x[x$size == s,]$change)
        }
        #create each list of variants as its own table
        dfl = lapply(bd_change, function(x) data.frame(table(x)))
        #remove any that are entirely zero
        dfl = lapply(dfl, function(x) x[x$Freq >0,])
        #change the names of the table so they can be combined together
        dfl = lapply(names(dfl), function(x) setNames(dfl[[x]], c("change",x)))
        #reduce the tables to combine together into one larger table
        dfl = Reduce(function(...) merge(..., all=TRUE), dfl)
        #remove NAs from table (variants that one sample had and another did not)
        dfl[is.na(dfl)]=0
        #return table to user
        return(dfl)
} 


box_graph = function(bp_table, s=0){
	if(s == 1){
		dfm = melt(t(bp_table))
		names(dfm) = c("Change", "Al", "Counts")
		dfm$Change = gsub("$","bp",dfm$Change)
		dfm$Change = factor(dfm$Change, levels = unique(dfm$Change), ordered=TRUE)
	}
	
	else{
		dfm = melt(bp_table)
		names(dfm) = c("Change", "Al", "Counts")
	}
	
	dfm$Al = gsub(".*_A", "", dfm$Al)
	dfm$Al = as.factor(dfm$Al)
	
	g = ggplot(dfm, aes(Change,Counts) ) + geom_boxplot(aes(color= Al))+theme(axis.text.x=element_text(angle=90))
	return(g)
}


pois_stat = function(id_bd, m = "a"){
	dfm = melt(id_bd)
	names(dfm) = c("Change", "Al", "Counts")
	dfm$Al = gsub(".*_A", "", dfm$Al)
        dfm$Al = as.factor(dfm$Al)
	#dfm$Change = as.factor(dfm$Change)
	if(m =="i"){
		pois = glm(Counts ~ Change*Al, dfm, family = poisson)
	}
	if(m == "s"){
		dfm$Change = factor(dfm$Change)
		pois = glm(Counts ~ Change*Al, dfm, family = poisson)
	}
	else{
		pois = glm(Counts ~ Change+Al, dfm, family = poisson)
	}
	return(pois)

}

f_test = function(id_bd, type = "a"){
	dfm = melt(id_bd)
        names(dfm) = c("Change", "Al", "Counts")
        dfm$Al = gsub(".*_A", "", dfm$Al)
        dfm$Al = as.factor(dfm$Al)
	#dfm$Change = as.factor(dfm$Change)
        if(type =="i"){
                nova = aov(Counts ~ Change*Al, dfm)
        }
        else{
                nova = aov(Counts ~ Change+Al, dfm)
        }
        return(nova)


}




