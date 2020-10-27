import datetime
import re
import os
import sys


input_mice = sys.argv[1]
input_human = sys.argv[2]
output_mix = sys.argv[3]
		
		
		
		
header1 = "#!Mixed reference of mice and human\n"
header2 = "#!genome-builds GRC mice 38 GRC human 38\n"
header3 = "#!genome-releases 91 91\n"
header4 = "#!genome-date {}\n".format(str(datetime.date.today()))
header=[header1,header2,header3,header4]
with open(input_mice) as annotation1:
    with open(input_human) as annotation2:
        with open(output_mix, 'w') as outfile:
            outfile.writelines(header)
            for line in annotation1:
                if(not line.startswith('#!')):
                    outfile.write(re.sub('^','mice',line))
            for line in annotation2:
                if(not line.startswith('#!')):
                    outfile.write(re.sub('^','human',line))