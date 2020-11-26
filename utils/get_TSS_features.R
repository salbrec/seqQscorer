library(ChIPpeakAnno)

args <- commandArgs(trailingOnly = TRUE)
bed_file = args[1]
assembly = args[2]
out_file = args[3]

reads <- toGRanges(bed_file, format="BED", header=FALSE)

annoDataTSS = NULL
if (assembly == 'GRCh38') {
	data(TSS.human.GRCh38)
	annoDataTSS = TSS.human.GRCh38
}
if (assembly == 'GRCm38') {
	data(TSS.mouse.GRCm38)
	annoDataTSS = TSS.mouse.GRCm38
}

tss_anno = binOverFeature(reads, annotationData=annoDataTSS, radius=5000, nbins=5, select='nearest', FUN=length)

tss_anno = DataFrame(tss_anno)

tss_anno$tss_dist = as.numeric(rownames(tss_anno))
tss_anno$perc = tss_anno$reads / length(reads) * 100

write.table(tss_anno,file=out_file, sep="\t", row.names = FALSE)



