import pysam
import re
import csv
from Bio import SeqIO
import gzip
from collections import defaultdict
import sys
import os

cwd = os.getcwd()
R1_FILE=sys.argv[1]
project_name_mode=sys.argv[2]


discard_secondary_alignements = True

save = pysam.set_verbosity(0)	


barcodes_struct = {
	'BC_start':0,
	'BC_end':12,
	'UMI_start':13,
	'UMI_end':20
	}

def parse_barcodes(fastq_parser, query_name, read_barcodes, barcodes_struct):
	for fastq_R1 in fastq_parser:
		if '/' in fastq_R1.id:
			R1_id = fastq_R1.id[:fastq_R1.id.find("/")]
		else:
			R1_id = fastq_R1.id
		read_barcodes[R1_id]['XC'] = str(fastq_R1.seq)[barcodes_struct['BC_start']:barcodes_struct['BC_end']]
		read_barcodes[R1_id]['XM'] = str(fastq_R1.seq)[barcodes_struct['UMI_start']:barcodes_struct['UMI_end']]
		if(read_barcodes[R1_id]['XM']==''):
			sys.SystemExit('UMI empty for read {}.\n The barcode is: {}.\nWhole entry is:{}'.format(R1_id, fastq_R1.seq,fastq_R1))
		if (R1_id == query_name):
			return(fastq_parser,read_barcodes)
	return(fastq_parser,read_barcodes)
    


infile_bam = pysam.AlignmentFile(cwd + '/projects/' + project_name_mode + '/tmp/Aligned.out.bam', "rb")
pysam.set_verbosity(save)

fastq_parser = SeqIO.parse(gzip.open(R1_FILE, "rt"), "fastq")

outfile = pysam.AlignmentFile(cwd + '/projects/' + project_name_mode + '/tmp/R1_joined_map.out.bam', "wb", template=infile_bam)
pysam.set_verbosity(save)

read_barcodes = defaultdict(lambda :{'XC':'','XM':''})

for bam_read in infile_bam:
	if(discard_secondary_alignements & bam_read.is_secondary):
		continue
	if (bam_read.query_name) in read_barcodes:
		current_barcodes = read_barcodes.pop(bam_read.query_name)
		tags = bam_read.get_tags()
		tags.extend([
			('XC', current_barcodes['XC'],'Z'),
			('XM', current_barcodes['XM'],'Z')])
		bam_read.set_tags(tags)
	else:
		fastq_parser,read_barcodes = parse_barcodes(fastq_parser, bam_read.query_name, read_barcodes, barcodes_struct)
		if (bam_read.query_name) not in read_barcodes:
			raise SystemExit('Read {} from mapped file is missing in reference fastq file!'.format(bam_read.query_name))
			os.remove(cwd + '/projects/' + project_name_mode + '/tmp/R1_joined_map.out.bam')
		current_barcodes = read_barcodes.pop(bam_read.query_name)
		tags = bam_read.get_tags()
		tags.extend([
			('XC', current_barcodes['XC'],'Z'),
			('XM', current_barcodes['XM'],'Z')])
		bam_read.set_tags(tags)
	outfile.write(bam_read)