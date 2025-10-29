import pysam
import sys
from tqdm import tqdm

# Arguments
barcodes_path = sys.argv[1]
umis_path = sys.argv[2]
input_bam = sys.argv[3]
output_bam = sys.argv[4]
barcode_len = int(sys.argv[5])
umi_len = int(sys.argv[6])

count = 0

# Open files in streaming mode
with open(barcodes_path, 'r') as barcodes_file, \
     open(umis_path, 'r') as umis_file, \
     pysam.AlignmentFile(input_bam, "rb", check_sq=False) as infile_bam, \
     pysam.AlignmentFile(output_bam, "wb", template=infile_bam) as outfile:

    # Iterate line by line with progress bar
    for bam_read, barcode_line, umi_line in tqdm(
            zip(infile_bam.fetch(until_eof=True), barcodes_file, umis_file),
            desc="Processing BAM", unit="reads", leave=True):

        barcode = barcode_line[5:5+barcode_len]
        umi = umi_line[5:5+umi_len]

        tags = bam_read.get_tags()
        tags.extend([
            ('XC', barcode, 'Z'),
            ('XM', umi, 'Z')
        ])
        bam_read.set_tags(tags)
        outfile.write(bam_read)
        count += 1

    # Check for extra lines in barcode/UMI files
    barcode_extra = sum(1 for _ in barcodes_file)
    umi_extra = sum(1 for _ in umis_file)

    if barcode_extra > 0 or umi_extra > 0:
        sys.stderr.write(
            f"[ERROR] Mismatch in number of records!\n"
            f"Processed: {count} reads.\n"
            f"Remaining: barcodes={barcode_extra}, UMIs={umi_extra}\n"
        )
        sys.exit(1)

print(f"[OK] Processed {count} reads successfully, all counts match.")
