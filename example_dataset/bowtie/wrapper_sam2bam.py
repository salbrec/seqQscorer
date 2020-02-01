from sys import *
import utils
import os

sam_file_path = argv[1]
bam_file_path = argv[2]

sam_file_basename = os.path.basename(sam_file_path)
unsorted_bam = './%s.unsorted.bam'%(sam_file_basename)

view = ['samtools', 'view', '-bS', sam_file_path, '-o', unsorted_bam]
output, error = utils.getSystemCall(view)
print('view is done!')

sort = ['samtools', 'sort', unsorted_bam, '-o', bam_file_path]
output, error = utils.getSystemCall(sort)
print('sort is done!')

rm = ['rm', unsorted_bam]
output, error = utils.getSystemCall(rm)
print('removed unsorted bam file')
