#!/usr/bin/env python

import sys,re,os
import matplotlib, 

GATK_LOCATION =  ""
BAM_READCOUNT_LOCATION = ""
VT_LOCATION = ""
BEDTOOLS_LOCATION = ""
def main(bed_file, bam_file, ref_snp_vcf, eval_snp_vcf, ref_fasta):
    #run VT on eval file to make sure any indels are left shifted
    normalized_vcf_output = "temp.normalized.vcf.gz"
    cmd = [ VT_LOCATION, "-r", ref_fasta, eval_snp_vcf, "-o", normalized_vcf_output ]
    #run GATK evaluator and store output


    #run bed intersect by position and get the disjoint set only in reference
    missed_sites_file = "sites_only_in_ref.vcf.gz"
    cmd = [BEDTOOLS_LOCATION, "-v", "-a", ref_snp_vcf, "-b", eval_snp_vcf, ">", missed_sites_file]
    #prepare sites into bam-readcount format CHR START STOP
    bam_readcount_sites_file = prepare_sites_file_from_vcf(missed_sites_file)
    #run bam-readcount
    bam_readcount_output = "bam_readcount.output"
    cmd = [BAM_READCOUNT_LOCATION, "-l", bam_readcount_sites_file, "-f", ref_fasta, ">", bam_readcount_output]
    #generate graphs
    missing_sites_coverage_quality_png = generate_coverage_and_quality_graph(bam_readcount_output)











if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed_file", action="store", dest="bed_file", default=None)
    parser.add_argument("--bam_file", action="store", dest="bam_file")
    parser.add_argument("--ref_snp_vcf", action="store", dest="ref_snp_vcf")
    parser.add_argument("--eval_snp_vcf", action="store", dest="eval_snp_vcf")
    parser.add_argument("--ref_fasta", action="store", dest="ref_fasta")
    args=parser.parse_args()
    #FIXME add nice error messages if not enough files are supplied
    main(args.bed_file,args.bam_file, args.ref_snp_vcf, args.eval_snp_vcf, args.ref_fasta);

