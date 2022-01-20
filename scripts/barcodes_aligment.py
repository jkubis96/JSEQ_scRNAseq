import pickle
import pysam
import os
import sys

cwd = os.getcwd()
project_name_mode=sys.argv[1]

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

save = pysam.set_verbosity(0)	
infile_bam = pysam.AlignmentFile(cwd + '/projects/' + project_name_mode + '/tmp/R1_joined_map.out.bam', "rb")
pysam.set_verbosity(save)

outfile = pysam.AlignmentFile(cwd + '/projects/' + project_name_mode + '/tmp/Map_barcodes.out.bam', "wb", template=infile_bam)

mapping = load_obj(cwd + '/projects/' + project_name_mode + '/tmp/barcode_mapping')
barcode_ref = load_obj(cwd + '/projects/' + project_name_mode + '/tmp/barcode_ref')
barcode_ext_ref = load_obj(cwd + '/projects/' + project_name_mode + '/tmp/barcode_ext_ref')
unknown_barcodes = set()

for bam_read in infile_bam:
	barcode = bam_read.get_tag('XC')
	if barcode in barcode_ref:
		mapping[0][barcode]['count'] += 1
		outfile.write(bam_read)
		continue
	elif barcode in barcode_ext_ref:
		reference_barcode = mapping[1][barcode]['ref']
		mapping[1][barcode]['count'] += 1
		bam_read.set_tag('XC',reference_barcode,value_type='Z',replace=True)
		outfile.write(bam_read)
		continue
	else:
		if barcode in unknown_barcodes:
			mapping['unknown'][barcode]['count'] += 1
		else:
			mapping['unknown'][barcode] = {'count':1}
			unknown_barcodes.add(barcode)
		outfile.write(bam_read)

save_obj(obj=mapping, name=cwd + '/projects/' + project_name_mode + '/tmp/barcode_mapping_count')