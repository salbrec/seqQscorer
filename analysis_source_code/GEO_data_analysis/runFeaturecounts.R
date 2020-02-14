# given a bam file, this script runs featureCounts from the Rsubread package
# and wites the results into file based on the given "out_file" path

library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)
bam_file = args[1]
out_file = args[2]

fc <- featureCounts(files=bam_file, annot.inbuilt="hg38", GTF.attrType="gene_name", GTF.featureType="exon", isPairedEnd=FALSE)

write.table(fc$counts,file=out_file, sep="\t", row.names = TRUE, quote=FALSE)


