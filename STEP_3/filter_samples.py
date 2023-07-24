#!/usr/bin/env ./venv/bin/python3
##test script to itirate through all trees in ExomeDepth directory, check all CNV file and create summary based on files and their location to exclude 
##remember to activate conda environment first
##conda activate *name_of_env*

import os
import sys
import csv
import re

##create sample list array
ori_samples = []

with open("/path/to/list/of/samples.txt") as file:
    for line in file:
        ori_samples.append(line.rstrip())


samples_csv = [sample + ".csv" for sample in ori_samples]  

samples_bam = [sample + "_sorted_unique_recalibrated.bam" for sample in ori_samples]  

## pull out genes of interest from sample of interest's CNV files
directory = r'/path/to/working/directory/'
aligned_directory = r'/path/to/alignment/directory/'


original_stdout = sys.stdout 

with open("samples_to_exclude.csv", "w+") as out:
    sys.stdout = out
    print ('sample', 'location', 'number of cnvs', 'size_of_bam', sep = ',')
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename in samples_csv:
                f = os.path.join(root,filename)
                f2 = f.rstrip()
                with open(f) as myfile:
                    csvreader = pd.read_csv(myfile, sep = ',')
                    #csvreader['file_name'] = filename
                    #csvreader['file_location'] = f2
                    number_of_cnvs = len(csvreader)
                    if number_of_cnvs > 250:
                        print (filename, f2, number_of_cnvs, sep = ',')
    for root2, dirs2, files2 in os.walk(aligned_directory):
        for filename2 in files2:
            if filename2 in samples_bam:
                f_bam = os.path.join(root2,filename2)
                f2_bam = f_bam.rstrip()
                file_stats = os.stat(f_bam)
                if file_stats.st_size < 2563000000:
                    print (filename2, f2_bam, '-', file_stats.st_size, sep = ',')

##messy samples are excluded from analysis pipeline

sys.stdout = original_stdout


## find samples alignment location 
print ('Samples to remove due to high CNV counts or low quality BAM file:')

for line in open("samples_to_exclude.csv", "r"):
    output = line.split('.')
    samples = output[0]
    for root, dirs, files in os.walk(aligned_directory):
        for filename in dirs:
            if filename in samples:
                f = os.path.join(root,filename)
                f2 = f.rstrip()
                print(f2)
                
exit()
