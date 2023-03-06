library(GenomicFeatures)
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

reads = readPeakFile(bed_file)
if (seqlevelsStyle(reads)[1] != "UCSC") {
  seqlevelsStyle(reads) <- "UCSC"
}

readsAnno = annotatePeak(reads, tssRegion=c(-1000, 150), TxDb=txdb, verbose=TRUE)
percentages = data.frame(show(readsAnno))

write.table(percentages,file=out_file, sep="\t")



