# make sure you use python 3
# python RNAexome_preprocessing.py <other args>
# e.g. python RNAexome_preprocessing.py -t pe -b /Users/test_pipeline -o /Users/test_pipeline_res -m RNA0 --repeat no --subs no -s yes --dedup picard --adapter truseq -u eva.hulstaert@ugent.be

import subprocess
import os
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Make RNA exome preprocessing scripts (and run on command line or submit to slurm server)')

# Read arguments
parser.add_argument('-t', nargs=1, choices=['se','pe'], required=True, help='Single end (se) or paired end (pe) data')
parser.add_argument('-b', nargs=1, required=True, help='Base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-o', nargs=1, required=True, help='Directory where output should be created', metavar='output_dir')
parser.add_argument('-m', nargs=1, required=True, help='String match to select sample folders in base directory e.g. RNA0', metavar='string_match')
parser.add_argument('-s', nargs=1, choices=['yes','no'], required=True, help='(reverse) stranded sequencing or not (unstranded)?')
parser.add_argument('--dedup', nargs=1, choices=['no', 'picard'], required=True, help='Method of duplicate removal [no, picard]')
parser.add_argument('--polyA', nargs=1, choices=['yes','no'], default=['no'], help='PolyA trimming needed? [default: no]')
parser.add_argument('--adapter', nargs=1, choices=['truseq','no'], default=['no'], help='Adapter trimming needed? If so, which one [default: no]')
parser.add_argument('--subs', nargs=1, choices=['no','start'], default=['no'], help='Subsampling needed? [default: no] If so, when? At start or after clumpify duplicate removal')
#parser.add_argument('--seq', nargs=1, choices=['novaseq','hiseq','nextseq'], default=['nextseq'], help='Sequencing machine? [default: nextseq] Is important for duplicate removal', metavar=seq_machine)
parser.add_argument('-n', nargs=1, default=['0'], help='Nr of reads to subsample to. If you need the analysis to stop in order to determine the subsampling level first, do not enter a number here (or -n 0)')
parser.add_argument('--repeat', nargs=1, choices=['yes','no'], default=['no'], help='Repeated analysis? [default: no] If yes: general fastq copying and filtering will not be repeated')
parser.add_argument('-u', nargs=1, required=True, help='Submitter email address (used for error reporting)', metavar='user_email')

# Parse arguments
args = parser.parse_args()
print(args)
data_type = args.t[0]
base_dir = args.b[0].rstrip("/")
string_match = args.m[0]
email = args.u[0]
deduplication = args.dedup[0]
stranded = args.s[0]
adaptercontam = args.adapter[0]
polyAtrim = args.polyA[0]
#print(picoV2)
subs_timing = args.subs[0]
repeat_analysis = args.repeat[0]
subsample_to_nr = args.n[0] # nr of reads you want to subsample to (based on floor of wc -l divided by 4), will be ignored if no subsampling asked for in command line
#seq_machine = args.seq[0]

# python submit_circRNA_slurm.py -p CIRCexplorer2_STAR -t se -d clumpify -b $VSC_DATA_VO/RNA_seq_pipeline/data/KVB/NSQ_Run219b-32853823 -g hg38 -o $VSC_DATA_VO/RNA_seq_pipeline/data/KVB/NSQ_Run219b-32853823/output_circRNA -s H -u jasper.anckaert@ugent.be
star_index = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/STAR_index/Genome_hg38_spikes_chrIS_MTr45S_star_index"
gtf = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes.gtf"
chrom = "/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.chr_spikes.txt"
exon_bed = "/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/ensembl_bedregions/Homo_sapiens.GRCh38.91_exons_sorted_merged2.bed"
tasks = 1
mem = 40


def sbatch(job_name, command, index, mem = mem, tasks = tasks, workdir = '.', time='4:00:00', dep=''): 
	if dep != '':
		dep = ' --dependency=afterok:{} --kill-on-invalid-dep=yes'.format(dep) 
	printlines = [
		"#!/bin/bash",
		"",
		##uncomment lines below to do SBATCH submission on server       
		#"#SBATCH -J {}".format(job_name+'_'+dedup_mode),
		#"#SBATCH -D {}".format(workdir),
		#"#SBATCH --mem={}G".format(mem),
		#"#SBATCH --cpus-per-task={}".format(tasks), #nr of processors per task needed
		#"#SBATCH -t {}".format(time),
		#"#SBATCH --mail-user={}".format(email),
		#"#SBATCH --mail-type=FAIL",
		#"#SBATCH -o {1}/{2}/logs/{3:02d}_{0}.out".format(job_name, workdir, dedup_mode, index),
		#"#SBATCH -e {1}/{2}/logs/{3:02d}_{0}.err".format(job_name, workdir, dedup_mode, index),
		#"#SBATCH --test-only", #Uncomment this to test instead of directly submitting job
		"cd {}".format(workdir),
		""       
	]
	
	printlines.extend(command.split("; "))
	printjobfile('{0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, dedup_mode, index, job_name),printlines)
	
	##uncomment lines below to do SBATCH submission on server
	#sbatch_command = 'sbatch{4} {0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, dedup_mode, index, job_name, dep)
	#sbatch_response = subprocess.getoutput(sbatch_command) 
	#job_id = sbatch_response.split(' ')[3].strip() 
	#return job_id

def printjobfile(filename, printlines):
        with open(filename, 'w') as the_file:
                for line in printlines:
                        the_file.write(line+"\n")

### Copy from DATA to SCRATCH
def combinecopy(sampleID,dep='',jobsuffix=''):
	""" Copy from original directory to working directory + unzip
	Check if the fastq files per lane are already combined ([[ -s file ]]: True if file has a Size greater than zero),
	if not concatenate them
	"""
	if data_type == 'se':
		command = '[[ -s {0}/{2}.fastq.gz ]] && cp {0}/{2}.fastq.gz {1}/{2}.fastq.gz || cat {0}/*R1_*.fastq.gz > {1}/{2}.fastq.gz; '.format(input_dir, output_dir, sampleID)
	elif data_type == 'pe':
		command = '[[ -s {0}/{2}_1.fastq.gz ]] && cp {0}/{2}_1.fastq.gz {1}/{2}_1.fastq.gz || cat {0}/*R1_*.fastq.gz > {1}/{2}_1.fastq.gz; '.format(input_dir, output_dir, sampleID)
		command = command + '[[ -s {0}/{2}_2.fastq.gz ]] && cp {0}/{2}_2.fastq.gz {1}/{2}_2.fastq.gz || cat {0}/*R2_*.fastq.gz > {1}/{2}_2.fastq.gz; '.format(input_dir, output_dir, sampleID)
	command = command + 'gunzip {0}/{1}*.fastq.gz; '.format(output_dir, sampleID)
	job_id = sbatch('combinecopy_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)

### If not first analysis: unzip the fastq files for further use (if they are gzipped)
def preprepeat(sampleID,dep='',jobsuffix=''):
	""" If the quality filtered files are already present, unzip them (if still needed) for further analysis
	Test if the file is gzipped [[ -s file.gz ]]: True if file has a Size greater than zero -> gunzip
	Else: test whether the unzipped file exists (if it doesn't, the function will throw an error)
	"""
	if data_type == 'se':
		command = '[[ -s {0}/{1}_qc.fastq.gz ]] && gunzip {0}/{1}_qc.fastq.gz || head -1 {0}/{1}_qc.fastq; '.format(output_dir, sampleID)
	elif data_type == 'pe':
		command = '[[ -s {0}/{1}_1_qc.fastq.gz ]] && gunzip {0}/{1}_1_qc.fastq.gz || head -1 {0}/{1}_1_qc.fastq; '.format(output_dir, sampleID)
		command = command + '[[ -s {0}/{1}_2_qc.fastq.gz ]] && gunzip {0}/{1}_2_qc.fastq.gz || head -1 {0}/{1}_2_qc.fastq; '.format(output_dir, sampleID)
	job_id = sbatch('preprepeat_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep)
	
### FastQC
def fastqc(sampleID,suffix='',outDIR='.',subDIR='',dep='',jobsuffix=''): 
	# Build the command for fastqc 
	command = 'module purge; module load FastQC/0.11.8-Java-1.8; mkdir {}/{}{}; '.format(output_dir, subDIR, outDIR)
	if data_type == 'se':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}{4}.fastq; '.format(output_dir, subDIR, outDIR, sampleID, suffix)
	elif data_type == 'pe':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_1{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_2{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
	job_id = sbatch('fastqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)

def removebadqc(sampleID,suffix='',percentage=str(80),quality=str(19),dep='',jobsuffix=''):
	"""
	Keep only read pairs where each read has at least p% of bases with a phred score > q
	-q: minimum quality phred score
	-p: percentage that has the minimum quality phred score
	FASTQC_filter.py: Keep only read pairs where each read has at least p% of bases with a phred score > q
	fastq_quality_filter: remove reads with lower quality, keeping only reads that have at least p% of bases with a quality score of q or more
	-i: FASTA/Q input file
	-o: FASTA/Q output file
	-Q33: If Illumina encoding is >= 1.8, you need to provide this option to fastq_quality_filter. For Illumina encoding <=1.5, option is not required
	
	"""
	command = 'module purge; module load Biopython/1.72-foss-2018b-Python-3.6.6; '
	#command = command + 'mkdir {0}/{1}/{2:02d}_removebadqcout; '.format(output_dir, dedup_mode, index)
	if data_type == 'se':
		command = 'module purge; module load FASTX-Toolkit/0.0.14-intel-2018a; '
		command = command + 'fastq_quality_filter -p {3} -q {4} -Q33 -i {0}/{1}{2}.fastq -o {0}/{1}_qc.fastq; '.format(output_dir, sampleID, suffix, percentage, str(int(quality)+1))
		command = command + 'qc_lines=`wc -l {0}/{1}_qc.fastq`; echo {1} $qc_lines > {0}/qc_lines.txt; '.format(output_dir,sampleID)
	elif data_type == 'pe':
		#split fastq files per 2.5 M reads in order to make it more feasible to process the different subparts
		command = command + 'split -l 10000000 --additional-suffix=.fastq {0}/{1}_1{2}.fastq {0}/temp_1_; '.format(output_dir, sampleID, suffix)
		command = command + 'rm {0}/{1}_1{2}.fastq; '.format(output_dir, sampleID, suffix)
		command = command + 'split -l 10000000 --additional-suffix=.fastq {0}/{1}_2{2}.fastq {0}/temp_2_; '.format(output_dir, sampleID, suffix)
		command = command + 'rm {0}/{1}_2{2}.fastq; '.format(output_dir, sampleID, suffix)
		command = command + 'for suffix in $(find {0}/. -type f -name "*temp_1*" | awk -F \'temp_1_\' \'{{print $2}}\' | cut -d \'.\' -f1); '.format(output_dir, index)
		command = command + 'do; '
		command = command + '   python FASTQC_filter.py -r1 {0}/temp_1_$suffix.fastq -r2 {0}/temp_2_$suffix.fastq -p {1} -q {2} -o1 {0}/temp_qc_1_$suffix.fastq -o2 {0}/temp_qc_2_$suffix.fastq; '.format(output_dir, percentage, quality)
		command = command + '   rm {0}/temp_1_$suffix.fastq;    rm {0}/temp_2_$suffix.fastq; '.format(output_dir)
		command = command + 'done; '
		command = command + 'cat {0}/temp_qc_1_* > {0}/{1}_1_qc.fastq; rm {0}/temp_qc_1*; '.format(output_dir, sampleID)
		command = command + 'cat {0}/temp_qc_2_* > {0}/{1}_2_qc.fastq; rm {0}/temp_qc_2*; '.format(output_dir, sampleID)
		command = command + 'qc_lines=`wc -l {0}/{1}_1_qc.fastq`; echo {1} $qc_lines > {0}/qc_lines.txt; '.format(output_dir, sampleID)
	job_id = sbatch('removebadqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep, tasks=4)

def trimadapter(sampleID,suffix='',adapterR1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',adapterR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',dep='',jobsuffix=''):
	""" Remove adapters and trim reads
	-a: regular 3' adapter R1 (To reduce the number of falsely trimmed bases, the alignment algorithm requires that, by default, at least three bases match between adapter and read)
	-A: regular 3' adapter R2
	--minimum-length 20: discard reads shorter than 20
	-U: remove x bases from beginning (positive x) or end (negative x) of R2 read (-U 3 needed for picoV2 data)
	-l: shorten each read down to a certain length, use the --length option or the short version -l (we want to trim last nt off 76nt reads)
	-q x,y: quality-trim read with a threshold of x from 3' and with threshold of y from 5'
	-o & -p: output file for first read (-o) and second read (-p) of pair
	--pair-filter=any: (default) read pair is discarded if one of the reads (R1 or R2) fulfills the filtering criterion (e.g. one read is shorter than 20)
	"""
	command = 'module purge; module load cutadapt/1.18-intel-2018b-Python-3.6.6; '
	command = command + 'mkdir {0}/{1:02d}_trimout; '.format(output_dir, index)
	if data_type == 'se':
		command = command + 'cutadapt -a {4} --minimum-length=35 -o {0}/{1:02d}_trimout/{2}_adapter.fastq {0}/{2}{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
	elif data_type == 'pe':
		command = command + 'cutadapt --pair-filter=any -a {4} -A {5} --minimum-length=35 -o {0}/{1:02d}_trimout/{2}_1_adapter.fastq -p {0}/{1:02d}_trimout/{2}_2_adapter.fastq {0}/{2}_1{3}.fastq {0}/{2}_2{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
	### PolyA adapter trimming necessary?
	if polyAtrim=='yes':
		if data_type == 'se':
			command = command + 'cutadapt -a T{{100}} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_adapterpA.fastq {0}/{1:02d}_trimout/{2}_adapter.fastq; '.format(output_dir,index,sampleID) #remove polyA tails (adapter sequence is considered to be a sequence of 100A for the R2 read and a sequence of 100T for the R1 read)
			command = command + 'mv {0}/{1:02d}_trimout/{2}_adapterpA.fastq {0}/{2}{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
		elif data_type == 'pe':
			command = command + 'cutadapt --pair-filter=any -a T{{100}} -A A{{100}} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_1_adapterpA.fastq -p {0}/{1:02d}_trimout/{2}_2_adapterpA.fastq {0}/{1:02d}_trimout/{2}_1_adapter.fastq {0}/{1:02d}_trimout/{2}_2_adapter.fastq; '.format(output_dir,index,sampleID) #remove polyA tails (adapter sequence is considered to be a sequence of 100A for the R2 read and a sequence of 100T for the R1 read)
			command = command + 'mv {0}/{1:02d}_trimout/{2}_1_adapterpA.fastq {0}/{2}_1{3}_trim.fastq; mv {0}/{1:02d}_trimout/{2}_2_adapterpA.fastq {0}/{2}_2{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
	else: #no polyA trimmming necessary
		if data_type == 'se':
			command = command + 'mv {0}/{1:02d}_trimout/{2}_adapter.fastq {0}/{2}{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
		elif data_type == 'pe':
			command = command + 'mv {0}/{1:02d}_trimout/{2}_1_adapter.fastq {0}/{2}_1{3}_trim.fastq; mv {0}/{1:02d}_trimout/{2}_2_adapter.fastq {0}/{2}_2{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
	command = command + 'rm -r {0}/{1:02d}_trimout; '.format(output_dir, index) #rm temp dir
	job_id = sbatch('trimadapter_'+shortname+jobsuffix, command, index,workdir=output_dir,dep=dep,tasks=4)

### Subsampling
def subsample(sampleID,nreads,suffix='',subDIR='',dep='',jobsuffix=''):
	""" Seqtk subsampling
	Downsample fastq files to x number of reads (take min of all samples)
	-s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of code
	"""
	command = 'module purge; module load seqtk/1.3-foss-2018b; '
	#if subsample_to_nr == '0': #if no subsampling nr given
	#	command = command + 'echo \"Determine the subsampling level first! (i.e. floored minimum of wc-l fastq divided by 4)\"; exit 1; '
	if data_type == 'se':
		command = command + 'seqtk sample -s100 {0}/{1}{2}{3}.fastq {4} > {0}/{5}/{2}{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
	elif data_type == 'pe':
		command = command + 'seqtk sample -s100 {0}/{1}{2}_1{3}.fastq {4} > {0}/{5}/{2}_1{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
		command = command + 'seqtk sample -s100 {0}/{1}{2}_2{3}.fastq {4} > {0}/{5}/{2}_2{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
	job_id = sbatch('subsample_'+shortname+jobsuffix,command,index,workdir=output_dir,mem=40, time='8:00:00', dep=dep,tasks=2)

### Duplicate removal
def picarddupremcoord(sampleID,bamIN='Aligned.sortedByCoord.out.bam',bamOUT='Aligned.sortedByCoord.picard.bam',subdir='.',dep='',jobsuffix=''):
	""" Duplicate removal with Picard
	(When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates)
	ASSUME_SORT_ORDER=coordinate
	REMOVE_DUPLICATES=true: remove all duplicates (optical and sequencing)
	VALIDATION_STRINGENCY=SILENT: (default: STRICT) Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags)
	M: file to write duplication metrics to
	"""
	command = 'module purge; module load picard/2.18.5-Java-1.8.0_162; '
	command = command + 'mkdir {0}/{1}/{2:02d}_picardout; '.format(output_dir, dedup_mode, index)
	command = command + 'java -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates I={0}/{1}/{2:02d}_srout/{3} O={0}/{1}/{5:02d}_picardout/{4} ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M={0}/{1}/{5:02d}_picardout/{6}_picard_dup.metrics; '.format(output_dir,dedup_mode, alignindex,bamIN,bamOUT,index, sampleID)
	#'rm {}_srout/{}'.format(sampleID,bamOUT) #not needed anymore, just interesting to see how many duplicates picard still finds (picard_dup.metrics)
	command = command + 'mv {0}/{1}/{2:02d}_picardout/{3} {0}/{1}/{4:02d}_srout/{3}; '.format(output_dir, dedup_mode, index,bamOUT,alignindex)
	job_id = sbatch('removedupPicard_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep, time='8:00:00', tasks=4)

### Mapping and counting (STAR/HTSeq/Kalllisto)
def alignment(sampleID,suffix='',subDIR='',dep='', jobsuffix=''): 
	""" Paired end STAR alignment with output as bam (coordinate sorted)
	--outFileNamePrefix *: output folder
	--outSAMtype BAM SortedByCoordinate: outputs BAM file that is sorted by coordinate ("Aligned.sortedByCoord.out.bam")
	--outReadsUnmapped Fastx: unmapped reads output in separate fastq files (Unmapped.out.mate1/2)
	--twopassMode Basic: run STAR 2-pass mapping for each sample separately
	--outMultimapperOrder Random: outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments
	--outSAMmultNmax -1: max number of multiple alignments for a read that will be output to the SAM/BAM files (default -1: all alignments (up to outFilterMultimapNmax) will be output)
	--outFilterMultimapNmax 10: (default 10) max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped (counted as "mapped to too many loci")
	--outSAMprimaryFlag AllBestScore: For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in the FLAG (column 2 of the SAM).
	The unmarked alignment is selected from the best ones (i.e. highest scoring). outSAMprimaryFlag AllBestScore option will change default behavior and output all alignments with the best score as primary alignments (i.e. 0x100 bit in the FLAG unset)
	--outFilterScoreMinOverLread 0.66: alignment will be output only if its score normalized over read length (sum of mate lengths for PE reads) is higher than or equal to this value (default 0.66)
	--outFilterMatchNminOverLread 0.66: alignment will be output only if the number of matched bases normalized over read length (sum of mate lengths for PE reads) is higher than or equal to this value (default 0.66)
	--outFilterMatchNmin 20: alignment will be output only if the number of matched bases is higher than or equal to this value (default 0)
	"""
	command = 'module purge; module load STAR/2.6.0c-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_srout; '.format(output_dir, dedup_mode, index)
	#command = command + 'if [ {0} = . ] ; then echo {0}; else mkdir {0}; fi; '.format(subdir) #if subdir is specified, make directory for it (if it does not already exist)
	#'STAR --runThreadN 10 --outFileNamePrefix {}_srout/ --readFilesIn {}.fastq {}.fastq --genomeDir {} --sjdbGTFfile {} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag OneBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20'.format(sampleID,sampleID1,sampleID2,star_index,gtf)
	if data_type == 'se':
		command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/{3:02d}_srout/{4}. --readFilesIn {1}/{7}{4}{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,sampleID,suffix,gtf,subDIR)
	elif data_type == 'pe':
		command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/{3:02d}_srout/{4}. --readFilesIn {1}/{7}{4}_1{5}.fastq {1}/{7}{4}_2{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,sampleID,suffix,gtf,subDIR)
	command = command + 'gzip {0}/{1}/{2:02d}_srout/*Unmapped.*; '.format(output_dir, dedup_mode, index)
	job_id = sbatch('alignSTAR_'+shortname+jobsuffix, command, index, workdir=output_dir, time='8:00:00', mem=40,dep=dep,tasks=8) 

def sortbambyname(bamIN, bamOUT, dep='', jobsuffix=''):
	""" Some algorithms only run on name sorted bam files instead of coordinate sorted ones, this function makes the conversion from coo to name sorted.
	-n: sort by read names instead of chromosomal coordinates
	-o: output (sam, bam, cram format is deduced from filename extension -> make sure it ends on .bam for bam file
	"""
	command = 'module purge; module load SAMtools/1.8-intel-2018a; '
	command = command + 'samtools sort -o {0}/{1}/{2:02d}_srout/{3} -n {0}/{1}/{2:02d}_srout/{4}; '.format(output_dir,dedup_mode,alignindex,bamOUT, bamIN)
	job_id = sbatch('sortbambyname_'+shortname+jobsuffix,command, index, workdir=output_dir, dep=dep, mem=40,tasks=2)

def countshtseq(bamIN, dep='', jobsuffix=''):
	""" HTSeq quantification of STAR mapped bam files.
	--order name: (default) needs name sorted bam files
	--nonunique none: (default) if the union of features for each position in the read is > 1, the read (pair) is counted as ambiguous and not counted for any features +
	if read (pair) aligns to more than one location in reference, it is scored as alignment_not_unique (for each location)
	--stranded reverse: for PE reads, the second read has to be on same strand and the first read has to be on opposite strand (for stranded=yes it is vice versa)
	"""
	command = 'module purge; module load HTSeq/0.11.0-foss-2018b-Python-2.7.15; '
	command = command + 'mkdir {0}/{1}/{2:02d}_htout; '.format(output_dir, dedup_mode, index)
	if stranded == 'yes':
		command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1}/{2:02d}_srout/{3} {4} > {0}/{1}/{6:02d}_htout/{5}_htseq_counts.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,gtf,shortname,index)
	else:
		command = command + 'htseq-count --format bam --order name --nonunique none --stranded no {0}/{1}/{2:02d}_srout/{3} {4} > {0}/{1}/{6:02d}_htout/{5}_htseq_counts.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,gtf,shortname,index)
	#'rm {0}_srout/{1}'.format(sampleID,bamIN) #remove name sorted bam for memory reasons (not needed anymore in pipeline)
	job_id = sbatch('countshtseq_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep, tasks=2)

def idxstat(bamIN,dep='', jobsuffix=''):
	""" Perform samtools idxstats
	Retrieve and print stats in the index file corresponding to the input file. The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
	Before calling idxstats, the input BAM file should be indexed by samtools index.
	"""
	command = 'module purge; module load SAMtools/1.8-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_idxout; '.format(output_dir, dedup_mode, index)
	command = command + 'samtools index {0}/{1}/{2:02d}_srout/{3}; '.format(output_dir,dedup_mode,alignindex,bamIN) #index bam file
	command = command + 'samtools idxstats {0}/{1}/{2:02d}_srout/{3} > {0}/{1}/{4:02d}_idxout/{5}_idxstat.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,index,shortname)
	job_id = sbatch('idxstat_'+shortname+jobsuffix, command, index,workdir=output_dir, dep=dep)

def strandedness(sampleID, bamIN,dep='',jobsuffix=''):
	""" RSeQC to retrieve % correct strandedness.
	Grep the line that shows what fraction of reads is explained by fr-firststrand (reverse in htseq)
	(1+-,1-+,2++,2-- category, e.g. 1+- read 1 '+' mapped to + strand while gene is on '-' strand is what we expect for reverse stranded)
	"""
	command = 'module purge; module load Python/2.7.14-intel-2018a; module load bx-python/0.8.1-intel-2018a-Python-2.7.14; module load RSeQC/2.6.4-intel-2018a-Python-2.7.14; '
	# strandedness on BAM file (without duplicate removal) that was already created with STAR
	command = command + 'infer_experiment.py -r {4} -i {0}/{1}/{2:02d}_srout/{3} > {0}/{1}/{5}_RSeQC_output_all.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exon_bed,sampleID) #output all metrics
	command = command + 'out=`cat {0}/{1}/{2}_RSeQC_output_all.txt | grep "1+-" | cut -d":" -f2`; '.format(output_dir,dedup_mode,sampleID) #grep the percentage we are interested it (how many are on correct strand?)
	command = command + 'echo {2} $out > {0}/{1}/{2}_RSeQC_output.txt; '.format(output_dir,dedup_mode,sampleID)
	job_id = sbatch('strandedness_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep)
	
def zipfastq(dep='',jobsuffix=''):
	""" Gzip fastq files and change permissions """
	command = 'gzip {0}/{1}/*.fastq; gzip {0}/*.fastq; '.format(output_dir, dedup_mode) #gzip all fastq files
	command = command + 'chmod -R 774 {0}/{1}; chgrp -R gvandesompele_lab {0}/{1}; '.format(output_dir, dedup_mode)
	job_id = sbatch('zipfastq_'+shortname+jobsuffix,command,index, workdir=output_dir, tasks=2, time='8:00:00', dep=dep)

### Select the functions you want to use + make sure the dependencies match!

## Make sure you have a file with the samples you want (e.g. by command line: ls | grep "RNA" > listsamples.txt)
#samples = [line.rstrip() for line in open(origdir+"/"+sys.argv[2])]
#for samplename in samples:
for samplename in os.listdir(base_dir): ##alternative approach if you want all samples with RNA in the name immediately
	if os.path.isdir(os.path.join(base_dir,samplename)):
		if re.search(string_match,samplename):
			print(samplename)
			output_dir = args.o[0].rstrip("/")
			input_dir = base_dir+"/"+samplename
			output_dir = output_dir+"/"+samplename
			shortname=str.split(samplename,"-")[0] #retrieve only the first part of the name to use in job name
			
			if deduplication == 'picard':
				dedup_mode = 'dedup_picard'
			else:
				dedup_mode = 'dedup_none'
			
			# name of subdirectory to put files in
			if subs_timing == 'start':
				dedup_mode += '-subs_atstart'			
			else:
				dedup_mode += '-subs_none'
				
			#make subdir for each sample in working directory (if it does not exist yet)
			os.makedirs(output_dir, exist_ok=True)
			os.makedirs(output_dir+"/"+dedup_mode+"/scripts", exist_ok=True)
			
			index = 0
			sub_dir = ''
			
			#### General fastq generation (if needed)
			# copy and unzip everything to this working directory
			if repeat_analysis == 'no':
				# copy and unzip everything to this working directory
				copycombine_jobid = combinecopy(samplename)
				index += 1
				fastqc_jobid = fastqc(samplename,'','FASTQC_original', dep=copycombine_jobid, jobsuffix='_orig')
				index += 1
				fastq_suffix = ''
				
				#lentrim_jobid = cutadapt_l75(samplename, suffix=fastq_suffix,dep=copycombine_jobid)
		
				#index += 1
				if adaptercontam != 'no':
					if adaptercontam == 'truseq':
						adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
						adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
					trimadapt_jobid = trimadapter(samplename, suffix=fastq_suffix, adapterR1=adapter1, adapterR2=adapter2, dep=copycombine_jobid, jobsuffix='')
					fastq_suffix += '_trim'
					index += 1
				else:
					trimadapt_jobid = copycombine_jobid
				
				removebadqc_jobid = removebadqc(samplename, suffix=fastq_suffix, dep=trimadapt_jobid)
				index += 1
				fastq_suffix = '_qc'
				fastqcfilt_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, dep=removebadqc_jobid, jobsuffix='_filt') #fastq file not in subdir yet
				index +=1
				### end of general fastq generation
				
				# if deduplication == 'clumpify': #already determine strandedness on original fastq (for others, it will be determined after running STAR)
				# 	# Determine strandedness (first STAR will map fastq without duplicate removal, then strandedness determined)
				# 	strandcoordsorted_jobid = strandedness('Aligned.sortedByCoord.out.bam',dep=removebadqc_jobid,jobsuffix='')
				# 	index += 1
			
			else: #not first analysis (repeat_analysis = 'yes')
				#skip the steps above and immediately start with the other tasks (by putting the dependency as '')
				index = 20 #start from a new index (to make clear it is a repeat)
				fastq_suffix = '_qc'
				preprepeat_jobid = preprepeat(samplename)
				#preprepeat_jobid = ''
				index += 1
				removebadqc_jobid = preprepeat_jobid
			####
			
			if subs_timing == 'start':
				if subsample_to_nr == '0':
					continue #go to next samplename
				else:
					#sub_dir = ''
					subsstart_jobid = subsample(samplename, str(subsample_to_nr), suffix=fastq_suffix, subDIR=sub_dir, dep=removebadqc_jobid, jobsuffix='')
					fastq_suffix += '_subs'
					sub_dir = dedup_mode+'/'
					index += 1
			else:
				#sub_dir = ''
				subsstart_jobid = removebadqc_jobid
			
			#names of bam files
			bam_coord = samplename+'.Aligned.sortedByCoord.out.bam'
			bam_name = samplename+'.Aligned.sortedByName.out.bam'
			
			# Mapping with STAR of low quality reads filtered FASTQ
			star_jobid = alignment(samplename, suffix=fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, jobsuffix='')
			alignindex = index
			index +=1
			
			endfastq_jobid=star_jobid
			
			# RSeQC: to retrieve % correct strandedness (based on original BAM - no duplicate removal)i
			strandcoordsorted_jobid = strandedness(samplename, bam_coord,dep=star_jobid,jobsuffix='')
			index += 1
		
			if deduplication == 'picard':
				# update bam file names for further use
				bam_coord_nodedup = bam_coord
				bam_coord = samplename+'.Aligned.sortedByCoord.picard.bam'
				bam_name = samplename+'.Aligned.sortedByName.picard.bam'
				# duplicate removal based on coordinate sorted bam (Picard)
				picardcoord_jobid = picarddupremcoord(samplename, bamIN=bam_coord_nodedup, bamOUT=bam_coord, dep=star_jobid, jobsuffix='_picard') #dep=starsubs_jobid
				index += 1
				namesortedbam_jobid = sortbambyname(bamIN=bam_coord, bamOUT=bam_name, dep=picardcoord_jobid, jobsuffix='_picard')
				index += 1
				
			else: #no duplicate removal
				namesortedbam_jobid = sortbambyname(bamIN=bam_coord, bamOUT=bam_name, dep=star_jobid, jobsuffix='_nodedup')
				index += 1
			
			# HTSeq quantification (bam needs to be sorted by name)
			htseq_jobid = countshtseq(bamIN=bam_name, dep=namesortedbam_jobid, jobsuffix='')
			index += 1
			# Run idxstats (samtools)
			idxstat_jobid = idxstat(bamIN=bam_coord, dep=namesortedbam_jobid, jobsuffix='')
			index += 1

			# Gzip all fastq files
			zipfastq_jobid = zipfastq(dep=endfastq_jobid)

if (subs_timing == 'start') & (subsample_to_nr == '0'):
	sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/${sample}_1_qc.fastq; done > lines_qcfilt_fastq.txt")
