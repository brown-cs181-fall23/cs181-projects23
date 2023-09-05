#!/bin/bash

cd fastqs
FASTQS=$(ls *.fq)
cd ..

for i in $FASTQS
do
	bowtie2 -x ~/genomes/dmel_all_chromosome -U fastqs/${i} -S SAMS/${i}.sam

	samtools sort -o BAMS/${i}.sorted.bam SAMS/${i}.sam
done

featureCounts -t gene -g gene_id -C \
	-a ~/genomes/dmel-all-r6.42.gtf \
	-o abridged_counts.txt \
	BAMs/*.sorted.bam