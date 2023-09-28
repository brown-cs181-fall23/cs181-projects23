#!/bin/bash

# get a list of all 6 .fq files
cd fastqs
FASTQS=$(ls *.fq)
cd ..

for i in $FASTQS
do
	## STEP 1, run bowtie2: 
	# Usage: bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]

	## STEP 2, run samtools sort:
	# Usage: samtools sort [options...] [in.sam]

done 

## STEP 3, run featurecounts
# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...