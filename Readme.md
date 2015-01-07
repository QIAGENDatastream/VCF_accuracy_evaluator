###INSTALLATION

To install:

* Type "make" in the root to build the suite of programs used in this tool
* type pip install -r requirements.txt to install needed modules
* Python dependencies:filemagic, numpy, matplotlib
* Ubuntu Issues:
- numpy: "sudo apt-get install python-numpy" should give a precompiled version 
- matplolib: "sudo apt-get install python-matplotlib" should give the precompiled version
- can't find zlib.h: "sudo apt-get install zlib1g-dev"
* OSX issues
- filemagic: "brew install libmagic"
- installing Homebrew for Mac is outside the scope of this document [HomeBrew.sh](http://brew.sh/)
 
Requirements (shipped alongside and compiled automatically with the Makefile- but can be manually specified if already installed):

* bedtools
* bam-readcount
* vt
* tabix
* samtools to reference the fasta (not done automatically by tool yet)
* ftp://ftp-trace.ncbi.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.18/ provides suitable vcf and bed file     

###SKIPPING MAKE if you MEET MOST OR ALL REQUIREMENTS ALREADY:

There are static config variables at the top of accuracy_evaluator.py to specify the locations of various utilities.

    usage: accuracy_evaluator.py [-h] [--bed_file BED_FILE] [--bam_file BAM_FILE]
                                 [--log_level LOG_LEVEL] --ref_vcf REF_VCF
                                 --eval_vcf EVAL_VCF --ref_fasta REF_FASTA
                                 [--no-deep-compare] [--gatk] [--rough-and-graph]
    Take in a high quality snp/indel set and evaluate a call set for concordance
    optional arguments:
      -h, --help            show this help message and exit
      --bed_file BED_FILE   bed file to limit comparisons to certain regions-
                            optional (default: None)
      --bam_file BAM_FILE   bam file to use to get readcounts for missing sites
                            (default: None)
      --log_level LOG_LEVEL
                            INFO for basic, DEBUG for (very) detailed debug
                            output) (default: 20)
      --ref_vcf REF_VCF     file of SNPs or INDELs you know to be true in your
                            callset (default: None)
      --eval_vcf EVAL_VCF   file of SNPs and/or INDELs you want to compare to the
                            known true callset (default: None)
      --ref_fasta REF_FASTA
                            reference fasta used to call the bam (default: None)
      --no-deep-compare     Default, most detailed mode:break things down by
                            allele, snp, and indel (default: True)
      --gatk                run GATK GenotypeConcordance for a further check
                            (default: False)
      --rough-and-graph     used bedtools for positional/missing site bed file
                            creation then interrogate bam to generate graph if bam
                            is supplied (default: False)


###EXAMPLE OUTPUT:

    MIXED   POSITONAL_MATCH_WRONG_ALLELES   53
    MIXED   EXACT_MATCH 229
    MIXED   ONE_ALLELE_CORRECT  744
    MIXED   NOT_FOUND_IN_EVALUATION_CALLSET 458
    MIXED   TOTAL   1484
    MIXED   PERCENTAGE_EXACTLY_CORRECT  15.43
    MIXED   PERCENTAGE_COMPLETELY_MISSED    30.86
    INDEL   POSITONAL_MATCH_WRONG_ALLELES   1684
    INDEL   EXACT_MATCH 134674
    INDEL   ONE_ALLELE_CORRECT  6678
    INDEL   NOT_FOUND_IN_EVALUATION_CALLSET 30161
    INDEL   TOTAL   173197
    INDEL   PERCENTAGE_EXACTLY_CORRECT  77.76
    INDEL   PERCENTAGE_COMPLETELY_MISSED    17.41
    SNP POSITONAL_MATCH_WRONG_ALLELES   185
    SNP EXACT_MATCH 2599397
    SNP ONE_ALLELE_CORRECT  4086
    SNP NOT_FOUND_IN_EVALUATION_CALLSET 137734
    SNP TOTAL   2741402
    SNP PERCENTAGE_EXACTLY_CORRECT  94.82
    SNP PERCENTAGE_COMPLETELY_MISSED    5.02
    OVERALL PERCENTAGE_EXACTLY_CORRECT  93.77
    OVERALL EXACT_MATCH 2734300
    OVERALL TOTAL_SITES 2916083
