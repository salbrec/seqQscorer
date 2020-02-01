# Quality Feature Preprocessing


This directory includes subfolders for the different tools that are applied. The FastQ files can be downloaded from ENCODE.

### RAW Features (FastQC)

change directory: `./example_dataset/FastQC`:
```
mkdir ./reports
ls ../fastq_files/* | awk -F"/" '{print $3}' | awk -F"." '{print "mkdir ./reports/"$1}' > make_dirs.sh
sh make_dirs.sh

ls ../fastq_files/* | awk -F"/" '{print $3}' | awk -F"." '{print "fastqc ../fastq_files/"$1".fastq.gz --extract '-d' ./reports/"$1"/ -o ./reports/"$1"/"}' > run_FastqQC.sh
sh run_FastqQC.sh

mkdir ./report_summaries
ls ./reports/ | awk '{print "cp ./reports/"$1"/"$1"_fastqc/summary.txt ./report_summaries/"$1".txt"}' > copy_summaries.sh
sh copy_summaries.sh
```
Depending on the number of FastQ files this can still run on a usual laptop or desktop PC. However, we recommend to perform this on a high performance cluster using the appropriate scheduler.

### MAP Features (Bowtie2)

#### Mapping

change directory: `./example_dataset/bowtie`:
```
python wrapper_mapping_paired_end.py ../fastq_files/ENCFF137DWP.fastq.gz ../fastq_files/ENCFF033TZD.fastq.gz ~/path_to_human_genome_assembly_hg38/bowtie_index ./sam_files/ENCFF137DWP.sam ./mapping_statistics/ENCFF137DWP.txt
python wrapper_mapping_paired_end.py ../fastq_files/ENCFF165NJF.fastq.gz ../fastq_files/ENCFF226TLO.fastq.gz ~/path_to_human_genome_assembly_hg38/bowtie_index ./sam_files/ENCFF165NJF.sam ./mapping_statistics/ENCFF165NJF.txt
python wrapper_mapping_paired_end.py ../fastq_files/ENCFF677ULQ.fastq.gz ../fastq_files/ENCFF410PAM.fastq.gz ~/path_to_human_genome_assembly_hg38/bowtie_index ./sam_files/ENCFF677ULQ.sam ./mapping_statistics/ENCFF677ULQ.txt
python wrapper_mapping_paired_end.py ../fastq_files/ENCFF895UGJ.fastq.gz ../fastq_files/ENCFF910QTP.fastq.gz ~/path_to_human_genome_assembly_hg38/bowtie_index ./sam_files/ENCFF895UGJ.sam ./mapping_statistics/ENCFF895UGJ.txt

```
We highly recommend to perform the mapping on a high performance cluster using multiple cores per mapping. Also a specific genome assembly has to be provided for bowtie in particular.

#### SAM to BAM conversion
change directory: `./example_dataset/bowtie`:
```
python wrapper_sam2bam.py ./sam_files/ENCFF137DWP.sam ./bam_files/ENCFF137DWP.bam
python wrapper_sam2bam.py ./sam_files/ENCFF165NJF.sam ./bam_files/ENCFF165NJF.bam
python wrapper_sam2bam.py ./sam_files/ENCFF677ULQ.sam ./bam_files/ENCFF677ULQ.bam
python wrapper_sam2bam.py ./sam_files/ENCFF895UGJ.sam ./bam_files/ENCFF895UGJ.bam
```

### LOC Features (ChIPseeker)
change directory: `./example_dataset/read_annotation`:
```
python wrapper_reads_annotation.py ../bowtie/bam_files/ENCFF137DWP.bam human ./statistics/ENCFF137DWP.tsv
python wrapper_reads_annotation.py ../bowtie/bam_files/ENCFF165NJF.bam human ./statistics/ENCFF165NJF.tsv
python wrapper_reads_annotation.py ../bowtie/bam_files/ENCFF677ULQ.bam human ./statistics/ENCFF677ULQ.tsv
python wrapper_reads_annotation.py ../bowtie/bam_files/ENCFF895UGJ.bam human ./statistics/ENCFF895UGJ.tsv
```

### TSS Features (ChIPpeakAnno)
change directory: `./example_dataset/tss_annotation`:
```
python wrapper_tss_annotation.py ../bowtie/bam_files/ENCFF137DWP.bam human ./statistics/ENCFF137DWP.tsv
python wrapper_tss_annotation.py ../bowtie/bam_files/ENCFF165NJF.bam human ./statistics/ENCFF165NJF.tsv
python wrapper_tss_annotation.py ../bowtie/bam_files/ENCFF677ULQ.bam human ./statistics/ENCFF677ULQ.tsv
python wrapper_tss_annotation.py ../bowtie/bam_files/ENCFF895UGJ.bam human ./statistics/ENCFF895UGJ.tsv
```

