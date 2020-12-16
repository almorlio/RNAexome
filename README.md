# RNA exome preprocessing pipeline

## Scripts
### Preparation
#### BasespaceRunDownloader.py
##### Set-up
You first have to get an accession token from `https://developer.basespace.illumina.com/dashboard` > `Create a new Application` > Fill in the fields (doesn't matter what you fill in) > Go to your new project > Credentials > Your Access Token.

##### Downloading sequencing runs.
BaseSpaceRunDownloader works with python 2.7, on the HPC e.g. `ml Python/2.7.15-intel-2018b`.
```
$ BaseSpaceRunDownloader.py -h

  Usage: BaseSpaceRunDownloader.py [-h] -p PROJECT -a TOKEN [-d DIRECTORY]

  Options:
    -h, --help            show this help message and exit
    -p PROJECT, --project=PROJECT
                          Project Name: required
    -a ACCESSTOKEN, --accesstoken=ACCESSTOKEN
                          Access Token: required if BASESPACE_API_TOKEN env
                          variable is not set
    -d DIRECTORY, --directory=DIRECTORY
                          Where the project directory is created, [default: ./]
```

You can find the Project Name in the BaseSpace: e.g. NSQ_Run532.

Fill in the project name, accesstoken and directory and the download should automatically start.


### Main script
#### Usage
Make sure to run the pipeline with Python3
This pipeline makes use of several other tools that need to be installed or loaded first (change module loading in script to the version you have in case of submission to server)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - for quality control
* [FASTX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) - for quality filtering
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - for adapter trimming
* [Seqtk](https://github.com/lh3/seqtk/blob/master/README.md) - for downsampling FASTQ files
* [STAR](https://github.com/alexdobin/STAR) - for mapping
* [HTSeq](https://htseq.readthedocs.io/en/master/) - for quantification
* [Picard](https://broadinstitute.github.io/picard/) - for PCR duplicate removal
* [SAMtools](http://samtools.sourceforge.net) - for sorting BAM file and index statistics
* [RSeQC](http://rseqc.sourceforge.net) - for strandedness info

```bash
$ python RNAexome_preprocessing.py -h

usage: RNAexome_preprocessing.py [-h] -t {se,pe} -b base_dir -o output_dir -m
                                 string_match -s {yes,no} --dedup {no,picard}
                                 [--polyA {yes,no}] [--adapter {truseq,no}]
                                 [--subs {no,start}] [-n N]
                                 [--repeat {yes,no}] -u user_email

Make RNA exome preprocessing scripts (and run on command line or submit to
slurm server)

optional arguments:
  -h, --help            show this help message and exit
  -t {se,pe}            Single end (se) or paired end (pe) data
  -b base_dir           Base directory where the sample subdirectories are
                        located
  -o output_dir         Directory where output should be created
  -m string_match       String match to select sample folders in base
                        directory e.g. RNA0
  -s {yes,no}           (reverse) stranded sequencing or not (unstranded)?
  --dedup {no,picard}   Method of duplicate removal [no, picard]
  --polyA {yes,no}      PolyA trimming needed? [default: no]
  --adapter {truseq,no}
                        Adapter trimming needed? If so, which one [default:
                        no]
  --subs {no,start}     Subsampling needed? [default: no] If so, when? At
                        start or after clumpify duplicate removal
  -n N                  Nr of reads to subsample to. If you need the analysis
                        to stop in order to determine the subsampling level
                        first, do not enter a number here (or -n 0)
  --repeat {yes,no}     Repeated analysis? [default: no] If yes: general fastq
                        copying and filtering will not be repeated
  -u user_email         Submitter email address (used for error reporting)
```

This will create, for each sample in your base directory, individual jobscripts to run subparts of the pipeline which can be run on command line.

```bash
e.g. bash 00_combinecopy_RNA0001.sh
```
If you have a server that accepts slurm, you can comment out lines in the script to directly submit jobs with dependencies.

#### Output
In the output directory, there will be a number of files and folders created. The directory will look like this:
```bash
.
└── <samplename>                                                  # Main sample directory
    ├── dedup_picard-subs_none                                    # Subdirectory according to chosen deduplication and subsampling
    │   ├── 05_srout                                              # Directory that contains all BAM files
    │   │   ├── <samplename>.Aligned.sortedByCoord.out.bam        # Initial BAM file (mapping with STAR)
    │   │   ├── <samplename>.Aligned.sortedByCoord.picard.bam     # BAM after duplicate removal
    │   │   ├── <samplename>.Aligned.sortedByCoord.picard.bam.bai
    │   │   ├── <samplename>.Aligned.sortedByName.picard.bam
    │   │   ├── <samplename>.Log.final.out
    │   │   ├── <samplename>.Log.out
    │   │   ├── <samplename>.Log.progress.out
    │   │   ├── <samplename>.SJ.out.tab
    │   │   ├── <samplename>._STARgenome
    │   │   │   ├── ...
    │   │   ├── <samplename>._STARpass1
    │   │   │   ├── Log.final.out
    │   │   │   └── SJ.out.tab
    │   │   ├── <samplename>.Unmapped.out.mate1
    │   │   ├── <samplename>.Unmapped.out.mate1.gz
    │   │   ├── <samplename>.Unmapped.out.mate2
    │   │   └── <samplename>.Unmapped.out.mate2.gz
    │   ├── 07_picardout            
    │   │   └── <samplename>_picard_dup.metrics                     # Picard duplicate removal stats
    │   ├── 09_htout
    │   │   └── <samplename>_htseq_counts.txt                       # HTSeq quantification output
    │   ├── 10_idxout
    │   │   └── <samplename>_idxstat.txt                            # SAMtools index stats
    │   ├── <samplename>_RSeQC_output_all.txt                       # RSeQC strandedness information
    │   ├── <samplename>_RSeQC_output.txt
    │   └── scripts                                                 # Folder that contains all scripts to run
    │       ├── 00_combinecopy_<samplename>.sh
    │       ├── 01_fastqc_<samplename>_orig.sh
    │       ├── 02_trimadapter_<samplename>.sh
    │       ├── 03_removebadqc_RNA010101.sh
    │       ├── 04_fastqc_<samplename>_filt.sh
    │       ├── 05_alignSTAR_<samplename>.sh
    │       ├── 06_strandedness_<samplename>.sh
    │       ├── 07_removedupPicard_<samplename>_picard.sh
    │       ├── 08_sortbambyname_<samplename>_picard.sh
    │       ├── 09_countshtseq_<samplename>.sh
    │       ├── 10_idxstat_<samplename>.sh
    │       └── 11_zipfastq_<samplename>.sh
    ├── FASTQC_original                                               # FASTQC of original FASTQ files
    │   ├── <samplename>_1_fastqc.html
    │   ├── <samplename>_1_fastqc.zip
    │   ├── <samplename>_2_fastqc.html
    │   └── <samplename>_2_fastqc.zip
    ├── FASTQC_qc                                                     # FASTQC of adapter trimmed and filtered FASTQ files
    │   ├── <samplename>_1_qc_fastqc.html
    │   ├── <samplename>_1_qc_fastqc.zip
    │   ├── <samplename>_2_qc_fastqc.html
    │   └── <samplename>_2_qc_fastqc.zip
    ├── qc_lines.txt                                                  # number of lines remaining after adapter trimming and filtering (divide by 4 to get number of reads)
    ├── <samplename>_1.fastq.gz                                       # Original FASTQ file R1
    ├── <samplename>_1_qc.fastq.gz                                    # Adapter trimmed and filtered FASTQ file R1
    ├── <samplename>_2.fastq.gz                                       # Original FASTQ file R2
    └── <samplename>_2_qc.fastq.gz                                    # Adapter trimmed and filtered FASTQ file R2
```
