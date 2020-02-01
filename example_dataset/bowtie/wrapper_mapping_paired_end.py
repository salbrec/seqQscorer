from sys import *
import utils

fastq1 = argv[1]
fastq2 = argv[2]
genome_index = argv[3]
sam_file = argv[4]
output_file = argv[5]

n_cores = 62
for i in range(len(argv)):
	if argv[i] == '-cores':
		n_cores = int(argv[i+1])

# create call for FastQC
call = ['bowtie2', '-x', genome_index, '-1', fastq1, '-2', fastq2, '-S', sam_file, '-p', str(n_cores)]
output, error = utils.getSystemCall(call)

with open(output_file, 'w') as summary_file:
	summary_file.write(error)










