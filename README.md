# Wightman-et-al-2020-script
Crude Shell Scripts that was used to separate Synthetic Yeast Illumina DNA reads from a heterologous diploid strains


Bioinformatics pipeline to filter sequencing reads of heterozygous diploid S. cerevisiae
# Do once to build an “index” of reference (in index folder)
bowtie2-build -f ~/data/finalscripts/references/JB11consensus369nogaps.fasta ~/data/finalscripts/index/JB11consensus369nogaps.index &&
bowtie2-build -f ~/data/finalscripts/references/MH1000consensus369nogaps.fasta ~/data/finalscripts/index/MH1000consensus369nogaps.index

# Do every time in relevant strain folder:

# Creates a logfile and prints the total number of reads to this file
echo Total no. reads for strain > redwhitelogfile.txt &&
cat *_1.fq | wc -l >> redwhitelogfile.txt &&
echo plus >> redwhitelogfile.txt &&
 cat *_2.fq | wc -l >> redwhitelogfile.txt &&

# Uses bowtie2 aligner to map reads to the reference genome 
bowtie2 -x ~/data/finalscripts/index/MH1000consensus369nogaps.index -q -1 *_1.fq -2 *_2.fq -S RW2MH.sam.temp &&
bowtie2 -x ~/data/finalscripts/index/JB11consensus369nogaps.index -q -1 *_1.fq -2 *_2.fq -S RW2JB.sam.temp &&

# ‘Sorts’ the bowtie output 
samtools sort RW2JB.sam.temp -o RW2JB.sorted.sam.temp -O sam &&
samtools sort RW2MH.sam.temp -o RW2MH.sorted.sam.temp -O sam &&

# Pulls out all reads that map perfectly to either parent reference and adds the sam ‘header’. Grep 100 for BGI sequenced strains and grep 101 for parent
cat RW2MH.sorted.sam.temp | head -5 > RW2MH.MD100.sam.temp && cat RW2MH.sorted.sam.temp | grep MD:Z:100 >> RW2MH.MD100.sam.temp &&
cat RW2JB.sorted.sam.temp | head -5 > RW2JB.MD100.sam.temp && cat RW2JB.sorted.sam.temp | grep MD:Z:100 >> RW2JB.MD100.sam.temp &&

# Prints the number of 100% matching reads, to either parent genome, to the logfile
echo No. reads 100% matched to JB11 ie MD100 >> redwhitelogfile.txt &&
cat RW2JB.MD100.sam.temp | grep MD:Z:100 | wc -l >> redwhitelogfile.txt &&
echo No. reads 100% matched to MH1000 ie MD100 >> redwhitelogfile.txt &&
cat RW2MH.MD100.sam.temp | grep MD:Z:100 | wc -l >> redwhitelogfile.txt &&

#calculates the read depth of 100% matched reads
samtools depth RW2JB.MD100.sam.temp > RW2JB.MD100.depth.temp &&
samtools depth RW2MH.MD100.sam.temp > RW2MH.MD100.depth.temp &&
samtools depth RW2JB.sorted.sam.temp > RW2JB.depth.temp &&
samtools depth RW2MH.sorted.sam.temp > RW2MH.depth.temp && 

# Prints to logfile the depth information
echo number of nucleotides in JB11 reference >> redwhitelogfile.txt &&
cat RW2JB.depth.temp | wc -l >> redwhitelogfile.txt &&
echo number of nucleotides in MH1000 reference >> redwhitelogfile.txt &&
cat RW2MH.depth.temp | wc -l >> redwhitelogfile.txt &&
echo Sum of read depths mapped to JB11 >> redwhitelogfile.txt &&
awk '{sum+=$3;} END{print sum;}' RW2JB.depth.temp >> redwhitelogfile.txt &&
echo Sum of read depths mapped to MH1000 >> redwhitelogfile.txt &&
awk '{sum+=$3;} END{print sum;}' RW2MH.depth.temp >> redwhitelogfile.txt &&
echo Sum of read depths mapped to JB11 with MD filter >> redwhitelogfile.txt &&
awk '{sum+=$3;} END{print sum;}' RW2JB.MD100.depth.temp >> redwhitelogfile.txt &&
echo Sum of read depths mapped to MH1000 with MD filter >> redwhitelogfile.txt &&
awk '{sum+=$3;} END{print sum;}' RW2MH.MD100.depth.temp >> redwhitelogfile.txt &&

# Using the fact that paired reads have the same name, creates separate lists of read names that are either paired or ‘single’
cat RW2JB.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==1) {print $2}}' > RW2JB.MD100.single.name.temp &&
cat RW2MH.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==1) {print $2}}' > RW2MH.MD100.single.name.temp &&
cat RW2JB.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==2) {print $2}}' > RW2JB.MD100.paired.name.temp &&
cat RW2MH.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==2) {print $2}}' > RW2MH.MD100.paired.name.temp &&

# Prints to logfile the number of single and paired reads 
echo No. single reads mapping 100% to JB11 >> redwhitelogfile.txt &&
cat RW2JB.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==1) {print $2}}' | wc -l >> redwhitelogfile.txt &&
echo No. single reads mapping 100% to MH1000 >> redwhitelogfile.txt &&
 cat RW2MH.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==1) {print $2}}' | wc -l >> redwhitelogfile.txt &&
echo No. pairs mapping 100% to JB11, x2 to get no. reads >> redwhitelogfile.txt &&
 cat RW2JB.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==2) {print $2}}' | wc -l >> redwhitelogfile.txt &&
echo No. pairs mapping 100% to MH1000, x2 to get no. reads >> redwhitelogfile.txt &&
 cat RW2MH.MD100.sam.temp | tail -n+6 | awk '{print $1}' | sort | uniq -c | sort | awk '{if($1==2) {print $2}}' | wc -l >> redwhitelogfile.txt &&

# Uses the names of paired and single reads to pull out all information of the read from the sam file 
grep -wf RW2JB.MD100.single.name.temp RW2JB.MD100.sam.temp > RW2JB.MD100.single.temp &&
grep -wf RW2MH.MD100.single.name.temp RW2MH.MD100.sam.temp > RW2MH.MD100.single.temp &&
grep -wf RW2JB.MD100.paired.name.temp RW2JB.MD100.sam.temp > RW2JB.MD100.paired.temp &&
grep -wf RW2MH.MD100.paired.name.temp RW2MH.MD100.sam.temp > RW2MH.MD100.paired.temp &&

# Creates file with the names and sequences of paired and single reads
cat RW2JB.MD100.paired.temp | awk '{print $1 "\t" $10}' > RW2JB.MD100.paired.nameseq.temp &&
cat RW2MH.MD100.paired.temp | awk '{print $1 "\t" $10}' > RW2MH.MD100.paired.nameseq.temp &&
cat RW2JB.MD100.single.temp | awk '{print $1 "\t" $10}' > RW2JB.MD100.single.nameseq.temp &&
cat RW2MH.MD100.single.temp | awk '{print $1 "\t" $10}' > RW2MH.MD100.single.nameseq.temp &&

# Separates single and paired reads that are universal, ie map to both parents, or unique, ie map to only one parent reference
cat RW2JB.MD100.paired.temp RW2MH.MD100.paired.temp | sort | uniq -c | awk '{print $2,$11}' | sort | uniq -c | sort | awk '{if($1==2) {print $2 "\t" $3}}' > RW_univ.paired.nameseq.temp &&
cat RW2JB.MD100.paired.temp RW2MH.MD100.paired.temp | sort | uniq -c | awk '{print $2,$11}' | sort | uniq -c | sort | awk '{if($1==1) {print $2 "\t" $3}}' > RW_uniq.paired.nameseq.temp &&

# Finds universal single reads
cat RW2JB.MD100.single.temp RW2MH.MD100.single.temp | awk '{print $1,$10}' | sort | uniq -c |  awk '{if($1==2) {print $2 "\t" $3}}' > RW_singleVsingle.univ.temp &&
cat RW2JB.MD100.single.temp RW2MH.MD100.paired.temp | awk '{print $1,$10}' | sort | uniq -c |  awk '{if($1==2) {print $2 "\t" $3}}' > RW_singleVspaired.univ.temp &&
cat RW_singleVsingle.univ.temp RW_singleVspaired.univ.temp | sort | uniq -c | awk '{if($1==1) {print $2 "\t" $3}}' > RW_univ.single.nameseq.temp &&

# compares the universal single reads to ALL single reads belonging to JB11, takes out the universal reads, and leaves the unique reads
cat RW_univ.single.nameseq.temp RW2JB.MD100.single.nameseq.temp | sort | uniq -c | awk '{if($1==1) {print $2 "\t" $3}}' > RW_uniq.single.nameseq.temp &&

# Prints to log file the number of paired or single unique and universal reads 
echo No. paired universal reads >> redwhitelogfile.txt &&
cat RW_univ.paired.nameseq.temp | wc -l >> redwhitelogfile.txt &&
echo No. paired unique, to JB11 OR MH1000 reads  >> redwhitelogfile.txt &&
cat RW_uniq.paired.nameseq.temp | wc -l >> redwhitelogfile.txt &&
echo No. single universal reads >> redwhitelogfile.txt &&
cat RW_univ.single.nameseq.temp | wc -l >> redwhitelogfile.txt &&
echo No. single unique, to JB11 OR MH1000 reads >> redwhitelogfile.txt &&
cat RW_uniq.single.nameseq.temp | wc -l >> redwhitelogfile.txt &&

# Creates file with the names and sequences of paired and single reads 
cat RW2JB.MD100.paired.temp | awk '{print $1 "\t" $10}' > RW2JB.MD100.paired.nameseq.temp &&
cat RW2MH.MD100.paired.temp | awk '{print $1 "\t" $10}' > RW2MH.MD100.paired.nameseq.temp &&
cat RW2JB.MD100.single.temp | awk '{print $1 "\t" $10}' > RW2JB.MD100.single.nameseq.temp &&
cat RW2MH.MD100.single.temp | awk '{print $1 "\t" $10}' > RW2MH.MD100.single.nameseq.temp &&

# Specifies which reads are unique to which parent genome 
cat RW2JB.MD100.paired.nameseq.temp RW_uniq.paired.nameseq.temp | sort | uniq -c | awk '{if($1==2) {print $2 "\t" $3}}' > RW_JB11.uniq.paired.nameseq.temp &&
cat RW2MH.MD100.paired.nameseq.temp RW_uniq.paired.nameseq.temp | sort | uniq -c | awk '{if($1==2) {print $2 "\t" $3}}' > RW_MH.uniq.paired.nameseq.temp &&
cat RW2JB.MD100.single.nameseq.temp RW_uniq.single.nameseq.temp | sort | uniq -c | awk '{if($1==2) {print $2 "\t" $3}}' > RW_JB11.uniq.single.nameseq.temp &&
cat RW2MH.MD100.single.nameseq.temp RW_uniq.single.nameseq.temp | sort | uniq -c | awk '{if($1==2) {print $2 "\t" $3}}' > RW_MH.uniq.single.nameseq.temp &&

# Prints to log file the number of paired and single unique reads 
echo No. reads that are: paired unique to JB11, paired unique to MH1000, all paired unique reads >> redwhitelogfile.txt &&
cat RW_JB11.uniq.paired.nameseq.temp | wc -l >> redwhitelogfile.txt && cat RW_MH.uniq.paired.nameseq.temp | wc -l >> redwhitelogfile.txt && cat RW_uniq.paired.nameseq.temp | wc -l >> redwhitelogfile.txt &&
echo No. reads that are: single unique to JB11, single unique to MH1000, all single unique reads >> redwhitelogfile.txt &&
cat RW_JB11.uniq.single.nameseq.temp | wc -l >> redwhitelogfile.txt && cat RW_MH.uniq.single.nameseq.temp | wc -l >> redwhitelogfile.txt && cat RW_uniq.single.nameseq.temp | wc -l >> redwhitelogfile.txt &&

# Uses the names of single unique reads and pulls that read, plus it’s pair mate from the sam file. Uses the names of paired reads to pull out all information from the sam file.
cat RW_JB11.uniq.single.nameseq.temp | awk '{print $1}' > RW_JB11.uniq.single.name.temp &&
grep -wf RW_JB11.uniq.single.name.temp RW2JB.sam.temp > RW_JB11.pulledmate.paired.sam.temp &&
cat RW_JB11.uniq.paired.nameseq.temp | awk '{print $1}' | sort | uniq -c | awk '{if($1==2) {print $2 }}' > RW_JB11.uniq.paired.name.temp &&
grep -wf RW_JB11.uniq.paired.name.temp RW2JB.sam.temp > RW_JB11.uniq.paired.sam.temp &&

# Combines these reads into a single file containing all the final, filtered reads
cat RW_JB11.pulledmate.paired.sam.temp RW_JB11.uniq.paired.sam.temp > RW_JB11_finalreads.temp &&

# The final section involves converting our ‘final read’ file into a ‘fastq’ file that can be imported into other programs such as Geneious.
paste -d ' ' - - - - <*_1.fq > RW_1.fastq.1read1line.temp &&
paste -d ' ' - - - - <*_2.fq > RW_2.fastq.1read1line.temp &&
cat RW_JB11_finalreads.temp | awk '{print $1}' | sort | uniq -c | awk '{print $2}' > RW_JB11_finalreads.names.temp &&
grep -wf RW_JB11_finalreads.names.temp RW_1.fastq.1read1line.temp > RW_finalreads1_1read1line.fastq.temp &&
grep -wf RW_JB11_finalreads.names.temp RW_2.fastq.1read1line.temp > RW_finalreads2_1read1line.fastq.temp  &&
cat RW_finalreads1_1read1line.fastq.temp RW_finalreads2_1read1line.fastq.temp | sort > RW_finalreads_read1line.fastq.temp &&
awk '{print $1 " \n" $2 "\n" $3 "\n" $4}' < RW_finalreads_read1line.fastq.temp > RW_finalreads.fastq  &&

# Prints to log file the final number of reads 
echo Final number of reads after all filtering >> redwhitelogfile.txt &&
cat RW_finalreads.fastq | wc -l >> redwhitelogfile.txt  &&

# Determines read depth of final reads 
cat RW2JB.sam.temp | head -5 > RW_JB11_finalreads_withheader.temp && cat RW_JB11_finalreads.temp >> RW_JB11_finalreads_withheader.temp  && 
samtools sort RW_JB11_finalreads_withheader.temp -o RW_JB11_finalreads_withheader_sorted.temp -O sam &&
samtools depth RW_JB11_finalreads_withheader_sorted.temp > RW_JB11_finalreads.depth.temp &&

# Prints to log file the information on the read depth of final reads
echo Sum of final read depths on JB11 >> redwhitelogfile.txt &&
awk '{sum+=$3;} END{print sum;}' RW_JB11_finalreads.depth.temp >> redwhitelogfile.txt
