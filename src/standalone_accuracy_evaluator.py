#!/usr/bin/env python

import sys,re,os
import matplotlib, 

GATK_LOCATION =  ""
BAM_READCOUNT_LOCATION = ""
VT_LOCATION = ""
BEDTOOLS_LOCATION = ""
def main():

    











if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed_file", action="store", dest="bed_file", default=None)
    parser.add_argument("--bam_file", action="store", dest="bam_file")
    parser.add_argument("--ref_snp_vcf", action="store", dest="ref_snp_vcf")
    parser.add_argument("--eval_snp_vcf", action="store", dest="eval_snp_vcf")
    args=parser.parse_args()
    main(args.bed_file,args.bam_file, args.ref_snp_vcf, args.eval_snp_vcf);

