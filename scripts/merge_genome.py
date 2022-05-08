import datetime 
import re
import os
import sys


input_mice = sys.argv[1]
input_human = sys.argv[2]
output_mix = sys.argv[3]
		
		

with open(input_mice, 'r') as annotation1:
    with open(input_human, 'r') as annotation2:
        with open(output_mix, 'w') as outfile:
            for line in annotation1:
                if(not line.startswith('#')):
                    outfile.write(re.sub('^','mice_',line))
            for line in annotation2:
                if(not line.startswith('#')):
                    outfile.write(re.sub('^','human_',line))