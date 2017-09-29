## GATK Best Practices

# 1. mapping and mark duplicates

# sort by genomic coordinate
java	–jar	picard.jar SortSam	INPUT=unsorted.sam	OUTPUT=sorted.sam	\	
SORT_ORDER=coordinate	

# mark	duplicates
java	–jar	picard.jar MarkDuplicatesWithMateCigar	\	
INPUT=unmarked.sam	OUTPUT=marked.sam

# split N trim?

# realign indels: running the indel realigner only at known sites

# Call variants in your seq data
java -jar GenomeAnalysisTK.jar \ 
-T HaplotypeCaller \ 
-R reference.fa \ 
-I preprocessed_reads.bam \  
-L 20 \ 
--genotyping_mode DISCOVERY \ 
-stand_emit_conf 10 \ 
-stand_call_conf 30 \ 
-o raw_variants.vcf 
