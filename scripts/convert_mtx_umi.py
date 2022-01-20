import os 
import subprocess
import sys

cwd = os.getcwd()
path=sys.argv[1]
tmp=sys.argv[2]


barcodes = {}
features = {}

os.mkdir(path)

out_barcodes = path + '/barcodes.tsv'
out_features = path + '/genes.tsv'
mtx = path + '/matrix.mtx'
temp_mtx = tmp + '/temp_umi.mtx'
header_mtx = tmp + '/header.mtx'
expression = tmp + '/umi_expression.long'
n_lines = 0
barcode_index = 1
feature_index = 1

with open(temp_mtx,'w') as mtx_stream:
		with open(expression,'r') as input_file:
			next(input_file) # skip first line
			for line in input_file:
				barcode,feature,count = line.strip().split('\t')
				if(barcode not in barcodes):
					barcodes[barcode] = barcode_index
					barcode_index += 1
				if(feature not in features):
					features[feature] = feature_index
					feature_index += 1
				mtx_stream.write('{} {} {}\n'.format(features[feature],barcodes[barcode],count))
				n_lines +=1

with open(out_barcodes,'w') as barcodes_outfile:
	for barcode in barcodes:
		barcodes_outfile.write('{}\n'.format(barcode))

with open(out_features,'w') as features_outfile:
	for feature in features:
		features_outfile.write('{}\n'.format(feature))
        
with open(header_mtx,'w') as header_outfile:
        header_outfile.write("%%MatrixMarket matrix coordinate real general\n")
        header_outfile.write('{} {} {}\n'.format(len(features), len(barcodes), n_lines))

subprocess.call("cat {} {} > {}".format(header_mtx, temp_mtx, mtx), shell=True)

os.remove(temp_mtx)
os.remove(header_mtx)