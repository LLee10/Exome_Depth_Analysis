#!/usr/bin/env bash

##Directory Set Up for ExomeDepth CNV query
##Li Ling

work_location="query_file" #edit
exomedepth_work_dir="/path/to/working/directory/${work_location}/"
aligned_work_dir="/path/to/alignment/directory/${work_location}/"

#Step1: Go to ExomeDepth directory
cd /path/to/working/directory/ 

#Step2: Create working directory
mkdir ${work_location}

#Step3: Copy ExomeDepth.Rmd file into directory
cp ExomeDepth.Rmd ./${work_location}/

#Step4: Create BAMlocation file in working directory by pulling out sample names in chosen Aligned directory
cd ${work_location}
echo list_of_bam_files >config_location.txt

# Go to Aligned directory to read file names
cd ${aligned_work_dir}
find -name "*_recalibrated.bam" >> ${exomedepth_work_dir}config_location.txt

# Back to working directory
cd ${exomedepth_work_dir}

#Step5: Check location file and create config.csv file to read into ExomeDepth, make sure it has header and list of samples
more config_location.txt | sed 's%./%%' > config.csv

#Step6: Create samplename.csv and bamname.csv file to read into ExomeDepth for CNV calling loop
more config_location.txt | sed 1d | awk -F'/' '{print $2}' > samplename.csv
more config_location.txt | sed 1d | awk -F'/' '{print $3}' > bamname.csv

#Step7: Attach list of samples onto main list of samples for gene/location query later
more samplename.csv > ../list_of_samples.txt

echo DONE, check ${exomedepth_work_dir}bamname.csv 
