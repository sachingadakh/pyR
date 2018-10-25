for a in 2.4 2.5 2.6 2.7 2.8
do 
	/home/ngs/bin/STAR-master/bin/Linux_x86_64/STAR --runThreadN 4 --runMode alignReads --twopassMode Basic --sjdbGTFfile ../starIndexdata/dm6.gtf --genomeDir ../starIndexdata/ --outSAMstrandField intronMotif --readFilesIn ${a}-E_R1trimmed.fastq ${a}-E_R2trimmed.fastq --outSAMtype BAM Unsorted SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ./star/remaining/$a
done
