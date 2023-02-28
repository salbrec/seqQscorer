library(GenomicFeatures)
library(ChIPpeakAnno)
library(ChIPseeker)

args <- commandArgs(trailingOnly = TRUE)

bed_file = args[1]
gtf_file = args[2]
out_file_base = args[3]

# create TxDb from given gtf file
custom_TxDb = makeTxDbFromGFF(gtf_file)
if (seqlevelsStyle(custom_TxDb)[1] != 'UCSC') {
  seqlevelsStyle(custom_TxDb) <- "UCSC"
}

# run the TSS annotations to derive the TSS features
reads <- toGRanges(bed_file, format="BED", header=FALSE)
if (seqlevelsStyle(reads)[1] != 'UCSC') {
  seqlevelsStyle(reads) <- "UCSC"
}

annoDataTSS <- toGRanges(custom_TxDb, feature="gene")

tssAnno = binOverFeature(reads, annotationData=annoDataTSS, radius=5000, nbins=5, select='nearest', FUN=length)
tssAnno = DataFrame(tssAnno)

tssAnno$tss_dist = as.numeric(rownames(tssAnno))
tssAnno$perc = tssAnno$reads / length(reads) * 100

tss_out = paste(out_file_base, ".TSS", sep="")
write.table(tssAnno,file=tss_out, sep="\t", row.names = FALSE)


# run the peak annotation to derive the LOC features
readsAnno = annotatePeak(reads, tssRegion=c(-1000, 150), TxDb=custom_TxDb, verbose=TRUE)

percentages = data.frame(show(readsAnno))

loc_out = paste(out_file_base, ".LOC", sep="")
write.table(percentages,file=loc_out, sep="\t")




