"""Provides the full preprocessing needed for the blocklist seqQscorer extension

"python deriveBLfeatures.py --help" will display a formatted help text on
the console. A comprehensive description is provided in the GitHub README
that includes examples as well. In short, this script derives the quality
features used by seqQscorer to perform automatic NGS quality control.

date:	2025-08-26
author:	Steffen Albrecht

"""

import argparse
import pandas as pd
from os.path import exists
from os import system, stat, makedirs
from sys import argv
import utils.utils as utils

script_dir = './'
if argv[0].find('/') >= 0:
	script_dir = argv[0][: - argv[0][::-1].find('/')]

parser = argparse.ArgumentParser(description='Blocklist features Preprocessing - derives feature sets needed leveraged by seqQscorer')
parser.add_argument('--fastq', '-i', type=str, required=True, help='Input fastq file. For paired-end sequencing, use Read 1.')
parser.add_argument('--assembly', '-a', type=str, required=True, choices=['hg38', 'mm10'], help='Species assembly defines the blocklist to be used. (Has to be consistent with the species used in for Bowtie2)')
parser.add_argument('--btidx','-ix', type=str, default=None, help='Filename prefix for Bowtie2 Genome Index (minus trailing .X.bt2). If not used, the blocklist-restricted Bowtie2 index is used from this repository.')
parser.add_argument('--outdir', '-o', type=str, default='./features_BL/', help='Output directory. Default: "./features_BL/"')
parser.add_argument('--cores', '-c', type=int, default=1, help='Defines the number of processors (CPUs) to be used by bowtie2 and samtools. (decreases runtime)')
parser.add_argument('--name', '-n', type=str, default=None, help='By default the output files are named by the file name of --fastq. In order to change this to a custom name, use this option.')

# parse and pre-process command line arguments
args = parser.parse_args()
outdir = args.outdir if args.outdir[-1] == '/' else args.outdir + '/'
makedirs(outdir, exist_ok=True)
out_file_name = utils.getFileName(args.fastq)
if args.name != None:
	out_file_name = args.name

# check if a Bowtie2 index has been provided.
# If not, use the blocklist-restricted Bowtie2 index provided in this repository.
idx_Bowtie2 = args.btidx
if args.btidx == None:
	idx_Bowtie2 = script_dir + 'utils/idxBowtie2/BL/' + args.assembly

# read in the right blocklist
bl_file = '%sutils/blocklists/%s.bed'%(script_dir, args.assembly)
blocklist = utils.read_blocklist(bl_file)
blocklist_df = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','ID'])

# first read in the chomosome sizes. the sizes are not needed, 
# however, it provides a list of relevant chromosomes.
fp_sizes = '%sutils/chromosome_sizes/%s.tsv'%(script_dir, args.assembly)
chrom_sizes = pd.read_csv(fp_sizes, sep='\t', names=['chr', 'size'])
chrom_sizes = chrom_sizes.loc[ [not c in ['chrX', 'chrY'] for c in chrom_sizes['chr']] ]
chrom_size_map = dict(zip(chrom_sizes['chr'], chrom_sizes['size']))

# run Bowtie2 to map reads against the genome assembly
mapping_dir = '%smapping_data/'%(outdir)
mapping_stats = '%s%s_stats.txt'%(mapping_dir, out_file_name)
makedirs(mapping_dir, exist_ok=True)
bam_file = '%s%s.bam'%(mapping_dir, out_file_name)
if not exists(mapping_stats) or not exists(bam_file):
	print('Running the mapping now with Bowtie2...')
	bowtie2  = 'bowtie2 -p %d -x %s '%(args.cores, idx_Bowtie2)
	bowtie2 += '-U %s 2> %s | '%(args.fastq, mapping_stats)
	bowtie2 += 'samtools view -@ %d -Sb -o %s'%(args.cores, bam_file)
	system(bowtie2)

# TODO: check if mapping was successful

# convert bam to bed
bed_file = bam_file.replace('.bam', '.bed')
if not exists(bed_file):
	print('Converting BAM file to BED using BEDtools...')
	bamtobed = 'bedtools bamtobed -i %s > %s'%(bam_file, bed_file)
	system(bamtobed)

# TODO: check if bed file was created successfully

# read in bed file regions as summits
# the summits describe the center of the mapped read locations in the full BED
# important: by default a restricted genome assembly is used from this repository
# to speed up the mapping. If the user provides the path to a full Bowtie2 index, 
# the chromosome names must be in the format "chr1", "chr2", ..., "chr19", "chr20"
print('Reading in BED file to derive features...')
summits = dict((chrom, []) for chrom in chrom_size_map.keys())
if args.btidx == None:
	summits = dict((blID, []) for blID in blocklist_df['ID'])
with open(bed_file, 'r') as f:
	for line in f:
		line = line.strip().split('\t')
		chrom_or_BLreg = line[0] if args.btidx == None else 'chr' + line[0]
		summits[chrom_or_BLreg].append( int((int(line[1])+int(line[2]))/2.0) )
# count the reads overlapping with the blocklisted regions
count_BL_reg = utils.count_reads_in_regions(summits, blocklist, 
									chrom_size_map, args.btidx == None)	
bl_features_file = outdir + out_file_name + '.csv'
count_BL_reg.to_csv(bl_features_file, index=False)








