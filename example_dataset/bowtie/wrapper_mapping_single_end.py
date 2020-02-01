from sys import *
import utils

fastq = argv[1]
genome_index = argv[2]
sam_file = argv[3]
output_file = argv[4]

n_cores = 62
for i in range(len(argv)):
	if argv[i] == '-cores':
		n_cores = int(argv[i+1])

# create call for FastQC
call = ['bowtie2', '-x', genome_index, '-q', fastq, '-S', sam_file, '-p', str(n_cores)]
output, error = utils.getSystemCall(call)

with open(output_file, 'w') as summary_file:
	summary_file.write(error)







