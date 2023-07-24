---
title: "CNV Calling With Exome Depth"
author: "Li Ling Lee"
date: "31/05/2023"
output: html_document
---

## R Version 4.0.0

```{r}
#install.packages("ExomeDepth")
#install.packages("GenomeInfoDb")
#install.packages("Rsamtools")
library(ExomeDepth)
library(GenomeInfoDb)
library(Rsamtools)
```

###Set working directory to file directory
##Copy path from Files>More>Copy Folder Path
```{r}
setwd("/path/to/working/directory/")
```

##Set definition files location
```{r}
FILElocation="/path/to/working/directory/" 
```

## Import Modified hg38 exons, conrad, omim and gnomad SV definitions
#  The exons have ~5% of the total exons removed due to poor mappability and repetitive sequences
```{r}
exons.hg38 <- as.data.frame(read.csv2(file = file.path(FILElocation,"exons.hg38.csv"), header = TRUE, sep = ","))
```

```{r}
conrad <- as.data.frame(read.csv2(file = file.path(FILElocation,"conrad.hg38.csv"), sep = ","))
```

```{r}
conrad.GRanges <- GenomicRanges::GRanges(seqnames = conrad$chromosome, 
											IRanges::IRanges(start=conrad$start,end=conrad$end), 
											names = conrad$cnv)
```

```{r}
omim <- as.data.frame(read.csv2(file = file.path(FILElocation,"omim_phenotypes.hg38.csv"), sep = ","))
```

```{r}
omim.GRanges <- GenomicRanges::GRanges(seqnames = omim$chromosome, 
											IRanges::IRanges(start=omim$start,end=omim$end), 
											names = omim$OMIM)
```

```{r}
gnomad <- as.data.frame(read.csv2(file = file.path(FILElocation,"gnomad.hg38.csv"), sep = ","))
```

```{r}
gnomad.GRanges <- GenomicRanges::GRanges(seqnames = gnomad$chromosome, 
											IRanges::IRanges(start=gnomad$start,end=gnomad$end), 
											names = gnomad$freq)
```

##Set bam location; most likely Mimir2 or Mimir, copy the windows filepath by clicking on it and change the backslashes to forward slashes \ ---> /
```{r}
BAMlocation="/path/to/alignment/directory/" 
```

# Read in the config.csv file containg our bam sample list
```{r}
analysisConfig <- read.csv('config.csv', 
							              header = TRUE, 
							              fill = TRUE)

list_of_samples <- as.vector(analysisConfig$list_of_bam_files)
```

#collate file path and get name of where our bams are located
```{r}
filenames <- file.path(BAMlocation, 
						paste0(analysisConfig$list_of_bam_files))

file.exists(filenames)
```
# Read in the samplename.csv and bamname.csv file containg our sample list
```{r}
samplename <- read.csv('samplename.csv', header = FALSE)

bamname <- read.csv('bamname.csv', header = FALSE)
```

# Perform read counting of exonic intervals - modified exon list
```{r}
my.counts <- getBamCounts(bed.frame = exons.hg38,
                          bam.files = filenames,
                          include.chr = FALSE)

print(head(my.counts))
```

# Save count table just in case
```{r}
write.table(my.counts,
			file="read_counts.txt")

##If R crash after reading counts (boo) load saved count file
#my.counts <- read.table("fold_congenica_read_counts.txt", header = TRUE)
```

# Create a dataframe and also a matrix of the exonic read counts
```{r}
ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

#check the fields and columns in dataframe:
print(head(ExomeCount.dafr))

# Create matrix of the bam counts only
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), 
							pattern = '*.bam')])

##check the fields and columns in matrix:
print(head(ExomeCount.mat))
```

# Create sets of 10 samples for each CNV calling iteration 
```{r}
##define number of subgroups we want to create
for (number_of_samples in bamname) {
  groups = round(length(number_of_samples)/10)
}

##split bamname and samplename
##create BamSet and SampleSet for each analysis run
for (i in samplename) {
  SampleSet = split(i, cut(seq_along(i), breaks = groups, labels = FALSE))
}

for (j in bamname) {
  BamSet = split(j, cut(seq_along(j), breaks = groups, labels = FALSE))
}
```

#Create range to allow each set of analysis to run 
```{r}
for (k in 1:groups){
  my_range = 1:k
}
```

# Perform ExomeDepth CNV calling 
```{r}
for(i in my_range){

  ## Create matrix and dataframe for each set
  my_sampleset = paste0(SampleSet[[i]])
  my_matrix = as.matrix(ExomeCount.dafr[, BamSet[[i]]])
  my_dataframe = as.data.frame(ExomeCount.dafr[, c("chromosome","start","end","exon", BamSet[[i]])])
  
  
  nsamples <- ncol(my_matrix)
  
  print(head(nsamples))
  
  message('Now looping over all the samples innit')
  
  # Perform ExomeDepth CNV calling on each sample with a loop
  for (l in 1:nsamples) {
    
    my.test.data <- as.matrix(my_matrix[, l])
    
    my.reference.set <- as.matrix(my_matrix[, -l])
  
    my.choice <- select.reference.set(test.counts = my.test.data,
                                        reference.counts = my.reference.set,
                                        bin.length = (my_dataframe$end - my_dataframe$start)/1000,
                                        n.bins.reduced = 10000)
    
    
    head(my.choice)
    
    my.matrix <- as.matrix(my_dataframe[, my.choice$reference.choice, drop = FALSE])
    
    my.reference.set <- apply(X = my.matrix,
                                    MAR = 1,
                                    FUN = sum)
    
    
    ## CNV calling
    
    all.exons <- new('ExomeDepth',
                      test = my_matrix[,l],
                      reference = my.reference.set,
                      formula = 'cbind(test, reference) ~ 1')
    
    #We can now call the CNV by running the underlying hidden Markov model:
    
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = my_dataframe$chromosome,
                          start = my_dataframe$start,
                          end = my_dataframe$end,
                          name = my_dataframe$exon)
    
    #check output
    head(all.exons@CNV.calls)
    
    
    print(head(all.exons@CNV.calls))
    
    #now annotating with exon/gene level information / conrad et al. and OMIM 
    
    
     exons.hg38.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg38$chromosome, 
    											IRanges::IRanges(start=exons.hg38$start,end=exons.hg38$end), 
    											names = exons.hg38$name)
    
    #here the minimum overlap should be very close to 0  
    all.exons <- AnnotateExtra(x = all.exons, 
    						reference.annotation = exons.hg38.GRanges, 
    						min.overlap = 0.0001, 
    						column.name = 'exons.hg38')
    
    all.exons <- AnnotateExtra(x = all.exons, 
    						reference.annotation = conrad.GRanges, 
    						min.overlap = 0.5, 
    						column.name = 'conrad.cnv')
    
    all.exons <- AnnotateExtra(x = all.exons, 
    						reference.annotation = omim.GRanges, 
    						min.overlap = 0.0001, 
    						column.name = 'OMIM')
    
    all.exons <- AnnotateExtra(x = all.exons, 
    						reference.annotation = gnomad.GRanges, 
    						min.overlap = 0.0001, 
    						column.name = 'gnomAD_0.0001')
    
    all.exons <- AnnotateExtra(x = all.exons, 
    						reference.annotation = gnomad.GRanges, 
    						min.overlap = 0.1, 
    						column.name = 'gnomAD_0.1')
    
    #Now save it in an easily readable format
    
    output.file <- paste(my_sampleset[l], 
    					          '.csv', sep = '')
    
    write.csv(file = output.file, 
    		      x = all.exons@CNV.calls, 
    		      row.names = FALSE)
    
    }
}

```
