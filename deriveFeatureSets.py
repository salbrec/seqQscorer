"""Provides the full preprocessing needed for seqQscorer

"python deriveFeatureSets.py --help" will display a formatted help text on 
the console. A comprehensive description is provided in the GitHub README 
that includes examples as well. In short, this script derives the quality
features used by seqQscorer to perform automatic NGS quality control.

date:	2020-11-26
author:	Steffen Albrecht

"""

import os
import sys
import argparse
from subprocess import call
import subprocess
import locale

import utils.Exceptions as myExceptions

def getSystemCall(call):
	process = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	out = out.decode(locale.getdefaultlocale()[1])
	err = err.decode(locale.getdefaultlocale()[1])
	if process.returncode:
		print("call failed, call was: %s" % ' '.join(call), file=sys.stderr)
		print("Message was: %s" % str(out), file=sys.stderr)
		print("Error code was %s, stderr: %s" % (process.returncode, err), file=sys.stderr, end='')
		raise Exception('runSystemCall Exception') 
	return out, err

def getFileName(file_path):
	# TODO: use endswith string functions.
	file_name = file_path[ -file_path[::-1].find('/') : ]
	if file_name[-3:] == '.gz':
		file_name = file_name[:-3]
	if file_name[-6:] == '.fastq':
		file_name = file_name[:-6]
	elif file_name[-3:] == '.fq':
		file_name = file_name[:-3]
	return file_name

script_dir = './'
if sys.argv[0].find('/') >= 0:
	script_dir = sys.argv[0][: - sys.argv[0][::-1].find('/')]

parser = argparse.ArgumentParser(description='seqQscorer Preprocessing - derives feature sets needed leveraged by seqQscorer')
parser.add_argument('--fastq1', '-f1', type=str, required=True, help='Input fastq file. Either the fastq file for a single-end sample or the fastq file for read 1 of a paired-end sample. Using this default destination "./feature_sets/", all feature sets are computed and saved. The file names define the sampleID while the file endings define the feature sets RAW, MAP, LOC and TSS.')
parser.add_argument('--fastq2', '-f2', type=str, default=None, help='In case of a paired-end sample, the fastq file for read 2. When the preprocessing is applied on paired-end samples, FastQC is applied to the read 1 by default. Also the sampleID for the output files are named by the file name of read1. These two aspects can be cahnged by using --fastqc respectively --name.')
parser.add_argument('--btidx','-ix', type=str, required=True, help='Filename prefix for Bowtie2 Genome Index (minus trailing .X.bt2).')
parser.add_argument('--outdir', '-o', type=str, default='./feature_sets/', help='Output directory. Default: "./feature_sets/"')
parser.add_argument('--cores', '-c', type=int, default=None, help='Defines the number of processors (CPUs) to be used by bowtie2 and samtools. (decreases runtime)')
parser.add_argument('--fastqc', '-f', type=int, default=1, choices=[1, 2], help='The fastq on which FastQC is applied on. Can optionally be selected for paired-end samples.')
parser.add_argument('--assembly', '-a', type=str, default='GRCh38', choices=['GRCh38', 'GRCm38'], help='Species assembly needed to define the gene structure / annotation used by the bioconductor functions. (has to be consistent with the species used in for Bowtie2)')
parser.add_argument('--gtf', '-g', type=str, default=None, help='File path for a gtf file to be used to get the LOC and TSS features. (--assembly will be ignored then)')
parser.add_argument('--name', '-n', type=str, default=None, help='By default the output files are named by the file name of --fastq1. In order to change this to a custom name, use this option.')

# parse and pre-process command line arguments
args = parser.parse_args()
outdir = args.outdir if args.outdir[-1] == '/' else args.outdir + '/'
if not os.path.exists(outdir):
	os.mkdir(outdir)
out_file_name = getFileName(args.fastq1)
if args.name != None:
	out_file_name = args.name

# run FastQC to derive the RAW features
feaures_RAW = '%s%s.RAW'%(outdir, out_file_name)
if not os.path.exists(feaures_RAW):
	input_fastq = args.fastq1
	if args.fastqc == 2:
		if args.fastq2 == None:
			raise myExceptions.WrongSettingException('Read for FastQC cannot be specified without providing a read 2 file.')
		input_fastq = args.fastq2
	fastqc_file_name = getFileName(input_fastq)
	print('Running FastQC on %s (read %s)...'%(fastqc_file_name, args.fastqc))
	
	reports_FastQC = '%sFastQC_report_%s/'%(outdir, out_file_name)
	os.mkdir(reports_FastQC) if not os.path.exists(reports_FastQC) else None

	fastqc_call = ['fastqc', input_fastq, '--extract', '-d', reports_FastQC, '-o', reports_FastQC]

	out, err = getSystemCall(fastqc_call)

	copy_summary = 'cp %s%s_fastqc/summary.txt %s'%(reports_FastQC, fastqc_file_name, feaures_RAW)
	os.system(copy_summary)
	print('FastQC report created!\n')
else:
	print('RAW features for %s exist already...'%(out_file_name))
	print('.. therefore the FastQC call is skipped.\n')

# run Bowtie2 to derive the MAP features
features_MAP = '%s%s.MAP'%(outdir, out_file_name)
mapping_dir = '%smapping_data_%s/'%(outdir, out_file_name)
bam_file_path = '%s%s.bam'%(mapping_dir, out_file_name)
if not os.path.exists(features_MAP) or not os.path.exists(bam_file_path):
	print('Running the mapping now with Bowtie2...')
	os.mkdir(mapping_dir) if not os.path.exists(mapping_dir) else None
	sam_file_path = bam_file_path.replace('.bam', '.sam')
	bowtie2_call = ['bowtie2', '-x', args.btidx]
	mapping_stats = ''
	if args.fastq2 == None:
		# run single-ended mapping
		bowtie2_call += ['-q', args.fastq1]
	else:
		# run paired-ended mapping
		bowtie2_call += ['-1', args.fastq1, '-2', args.fastq2]
	bowtie2_call += ['-S', sam_file_path]
	if args.cores != None:
		bowtie2_call += ['-p', str(args.cores)]
	
	out, err = getSystemCall(bowtie2_call)
	mapping_stats = err
	
	open(features_MAP, 'w').write(mapping_stats)
	print('Mapping is done!\n')
	print('Now creating the sorted bam file...')
	unsorted_bam = sam_file_path.replace('.sam', '.unsorted.bam')
	multi_thread = ['-@', str(args.cores)] if args.cores != None else []
	samtools_view = ['samtools', 'view', '-bS', sam_file_path, '-o', unsorted_bam] + multi_thread
	output, error = getSystemCall(samtools_view)
	samtools_sort = ['samtools', 'sort', unsorted_bam, '-o', bam_file_path] + multi_thread
	output, error = getSystemCall(samtools_sort)
	print('bam file processed!\n')
	
	os.system('rm ' + unsorted_bam)
	os.system('rm ' + sam_file_path)
else:
	print('MAP features and bam file already exist!')
	print('... therefore the Bowtie2 call is skipped.\n')

# preprocess bam file for ChIPseeker and ChIPpeakAnno
subsampled_1Mb_bed = bam_file_path.replace('.bam','_1Mb.bed')
if not os.path.exists(subsampled_1Mb_bed):
	print('Some preprocessing is needed before LOC and TSS features are derived.')
	print('Subsampling 1.000.000 mapped reads...')
	
	# converting bam to bed
	mapping_bed_file_path = bam_file_path.replace('.bam', '.bed')
	bedtools = ['bedtools', 'bamtobed', '-i', bam_file_path, '>', mapping_bed_file_path]
	os.system(' '.join(bedtools))
	
	# subsampling 1.000.000 reads
	reads_sampling = ['shuf', '-n', '1000000', mapping_bed_file_path, '>', subsampled_1Mb_bed]
	os.system(' '.join(reads_sampling))
	print('Preprocessing is done!\n')

if args.gtf == None:
	# run ChIPseeker to derive the LOC features
	features_LOC = '%s%s.LOC'%(outdir, out_file_name)
	if not os.path.exists(features_LOC):
		print('Running ChIPseeker to derive the LOC features...')
		ChIPseeker_call = ['Rscript', '%sutils/get_LOC_features.R'%(script_dir), subsampled_1Mb_bed, args.assembly, features_LOC]
		os.system(' '.join(ChIPseeker_call))
		print('LOC features done!\n')
	else:
		print('LOC features already exist!')
		print('... therefore the ChIPseeker call is skipped.\n')

	# run ChIPpeakAnno to derive the TSS features
	features_TSS = '%s%s.TSS'%(outdir, out_file_name)
	if not os.path.exists(features_TSS):
		print('Running ChIPpeakAnno to derive the TSS features...')
		ChIPpeakAnno_call = ['Rscript', '%sutils/get_TSS_features.R'%(script_dir), subsampled_1Mb_bed, args.assembly, features_TSS]
		os.system(' '.join(ChIPpeakAnno_call))
		print('TSS features done!\n')
	else:
		print('TSS features already exist!')
		print('... therefore the ChIPpeakAnno call is skipped.\n')
else:
	out_file_base = '%s%s'%(outdir, out_file_name)
	if not os.path.exists(out_file_base + '.LOC') or not os.path.exists(out_file_base + '.TSS'):
		print('Running ChIPseeker and ChIPpeakAnno with the given gtf file...')
		bioconductor_calls = ['Rscript', '%sutils/gtf_LOC_TSS_features.R'%(script_dir), subsampled_1Mb_bed, args.gtf, out_file_base]
		os.system(' '.join(bioconductor_calls))
		print('LOC and TSS features done!\n')
	else:
		print('LOC and TSS features already exist!')
		print('... therefore the ChIPseeker and ChIPpeakAnno calls are skipped.\n')

# this PDF is automatically created by one of the bioconductor functions
# it is not needed so it is automatically removed here
if os.path.exists('./Rplots.pdf'):
	os.system('rm ./Rplots.pdf')





