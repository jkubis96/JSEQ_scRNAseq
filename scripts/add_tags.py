import pysam 
import re
import csv
from collections import defaultdict
import sys
import os


barcodes_path=sys.argv[1]
umis_path=sys.argv[2]
input_bam=sys.argv[3]
output_bam=sys.argv[4]
barcode_len=sys.argv[5]
umi_len=sys.argv[6]

barcodes = list(open(barcodes_path, 'r'))
barcodes = [barcodes[n][5:int(barcode_len)+5] for n, x in enumerate(barcodes)]

umis = list(open(umis_path, 'r'))
umis = [umis[n][5:int(umi_len)+5] for n, x in enumerate(umis)]



save = pysam.set_verbosity(0)	



infile_bam = pysam.AlignmentFile(input_bam, "rb")
pysam.set_verbosity(save)



outfile = pysam.AlignmentFile(output_bam, "wb", template=infile_bam)
pysam.set_verbosity(save)

current_barcodes = dict()
current_barcodes['XC'] = barcodes
current_barcodes['XM'] = umis



for n, bam_read in enumerate(infile_bam):
    tags = bam_read.get_tags()
    tags.extend([
        ('XC', current_barcodes['XC'][n],'Z'),
        ('XM', current_barcodes['XM'][n],'Z')])
    bam_read.set_tags(tags)
    outfile.write(bam_read)