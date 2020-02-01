from sys import *
import os

assemblies = {'human': 'hg38', 'mouse': 'mm10', 'drosophila': 'dm6', 'celegans': 'ce10'}

bam_file_path = argv[1]
organism = argv[2]
output_file = argv[3]

# use the name of bam file for the temporary bed files
bam_file_basename = os.path.basename(bam_file_path)
mapping_bed_file_path = ('./%s.mapping.bed'%(bam_file_basename))
mapping_1M_bed_file_path = ('./%s.mapping_1M.bed'%(bam_file_basename))


print('Converting bam to bed ...')
bedtools = ['bedtools', 'bamtobed', '-i', bam_file_path, '>', mapping_bed_file_path]
bedtools = ' '.join(bedtools)
os.system(bedtools)

print('Sampling 1M reads from mapping.bed ...')
reads_sampling = ['shuf', '-n', '1000000', mapping_bed_file_path, '>', mapping_1M_bed_file_path]
reads_sampling = ' '.join(reads_sampling)
os.system(reads_sampling)

print('\nRunning Bioconda now ...')
annots = ['Rscript', './readsAnno.R', mapping_1M_bed_file_path, assemblies[organism], output_file]
annots = ' '.join(annots)
os.system(annots)

print('\nDeleting the temporary bed-files ...')
rm = ['rm', mapping_bed_file_path]
os.system(' '.join(rm))
rm = ['rm', mapping_1M_bed_file_path]
os.system(' '.join(rm))


