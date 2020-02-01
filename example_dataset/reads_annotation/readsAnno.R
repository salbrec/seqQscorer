library(ChIPseeker)

args <- commandArgs(trailingOnly = TRUE)
bed_file = args[1]
assembly = args[2]
out_file = args[3]

txdb = NULL
if (assembly == 'hg38') {
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
}
if (assembly == 'mm10') {
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
}
if (assembly == 'dm6') {
	library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
	txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
}
if (assembly == 'ce10') {
	library(BSgenome.Celegans.UCSC.ce10)
	txdb = BSgenome.Celegans.UCSC.ce10
}

readsAnno = annotatePeak(bed_file, tssRegion=c(-1000, 150), TxDb=txdb)
percentages = data.frame(show(readsAnno))

write.table(percentages,file=out_file, sep="\t")


