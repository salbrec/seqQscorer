library(ChIPseeker)

args <- commandArgs(trailingOnly = TRUE)
bed_file = args[1]
assembly = args[2]
out_file = args[3]

txdb = NULL
if (assembly == 'GRCh38') {
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
}
if (assembly == 'GRCm38') {
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
}
peaks = readPeakFile(bed_file)
# added a check, if the levels of seqnames contain "chr" 
# (onbly the first is checks since they are expected to be uniform)
# if chr is not present, is is added wit a helperfunction of the packages diffloop
if (!grepl("chr", levels(peaks@seqnames), fixed=TRUE)[1]){
	library(diffloop)
	peaks = addchr(peaks)
}

readsAnno = annotatePeak(peaks, tssRegion=c(-1000, 150), TxDb=txdb, verbose=TRUE)
percentages = data.frame(show(readsAnno))

write.table(percentages,file=out_file, sep="\t")



