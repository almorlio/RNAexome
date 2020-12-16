from Bio import SeqIO
import argparse
import itertools
import os

#import multiprocessing as mp


def filter_PEreads(fastqfile1,fastqfile2,out1,out2,qcutoff,pcutoff):
	fastq1 = SeqIO.parse(open(fastqfile1),"fastq")
	fastq2 = SeqIO.parse(open(fastqfile2),"fastq")
	outfile1 = open(out1,"w")
	outfile2 = open(out2,"w")
	#num_cores = multiprocessing.cpu_count()
	for read1 in fastq1:
		read2 = next(fastq2)
		read1q = read1.letter_annotations["phred_quality"]
		read2q = read2.letter_annotations["phred_quality"]
		qread1 = 100 * len([i for i in read1q if i > qcutoff]) / len(read1q)
		qread2 = 100 * len([i for i in read2q if i > qcutoff]) / len(read2q)
		if qread1 >= pcutoff and qread2 >= pcutoff:
			SeqIO.write(read1, outfile1, "fastq")
			SeqIO.write(read2, outfile2, "fastq")
	outfile1.close()
	outfile2.close()

if __name__ == '__main__':
	ap = argparse.ArgumentParser()
	ap.add_argument("-q","--quality", default=20, help="The minimum quality phred score")
	ap.add_argument("-p","--percentage", default=80, help="The percentage that has the minimum quality phred score")
	ap.add_argument("-r1","--read1", required=True, help="Fastq for read pair 1")
	ap.add_argument("-r2","--read2", required=True, help="Fastq for read pair 2")
	ap.add_argument("-o1","--outread1", required=True, help="Output filtered Fastq for read pair 1")
	ap.add_argument("-o2","--outread2", required=True, help="Output filtered Fastq for read pair 2")
	args = vars(ap.parse_args())
	filter_PEreads(args["read1"],args["read2"],args["outread1"],args["outread2"],int(args["quality"]),int(args["percentage"]))

