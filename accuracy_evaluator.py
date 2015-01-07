#!/usr/bin/env python

import sys,re,os, subprocess, logging
import magic, gzip
import argparse 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from collections import Counter as mset

#GATK will crash unless you use a particular version of java. you could change this variable to point to a different one 
#if your default system java is problematic
JAVA_LOCATION = "java"
GATK_LOCATION =  "./GATK/GenomeAnalysisTK.jar"
BAM_READCOUNT_LOCATION = "./bam-readcount/bin/bam-readcount"
VT_LOCATION = "./vt/vt"
BEDTOOLS_LOCATION = "./bedtools/bin/bedtools"
TABIX_LOCATION = "tabix"
logger = None

def rough_positional_intersect(ref_snp_vcf, eval_snp_vcf, bed_file=None):
    missed_sites_file = "sites_only_in_ref.vcf"
    cmd = [BEDTOOLS_LOCATION, "intersect", "-v", "-header", "-a", ref_snp_vcf, "-b", eval_snp_vcf, ">", missed_sites_file]
    logger.debug("Running:%s" % " ".join(cmd))
    subprocess.call(" ".join(cmd), shell=True)
    extra_sites_file = "sites_only_in_eval.vcf"
    cmd = [BEDTOOLS_LOCATION, "intersect", "-v", "-header", "-a", eval_snp_vcf, "-b", ref_snp_vcf, ">", extra_sites_file]
    logger.debug("Running:%s" % " ".join(cmd))
    subprocess.call(" ".join(cmd), shell=True)
    missed_sites = count_sites(missed_sites_file)
    extra_sites = count_sites(extra_sites_file)
    total_ref_sites = count_sites(ref_snp_vcf)
    total_eval_sites = count_sites(eval_snp_vcf)
    print "[POSITION_BASED]Total Sites in Reference/Canonical VCF: %s" % total_ref_sites
    print "[POSITION_BASED]Total Sites in Supplied/Evaluation VCF: %s" % total_eval_sites
    print "[POSITION_BASED]Missed Sites: %s" % missed_sites
    print "[POSITION_BASED]Extra Sites: %s" % extra_sites
    if bed_file != None:
        missed_sites_region_limited_file="sites_only_in_ref.region_limited.vcf"
        limit_cmd = [BEDTOOLS_LOCATION, "intersect", "-wa", "-u", "-header", "-a", missed_sites_file, "-b", bed_file, ">", missed_sites_region_limited_file]
        logger.debug("Running:%s" % " ".join(limit_cmd))
        subprocess.call(" ".join(limit_cmd), shell=True)
        extra_sites_region_limited_file="sites_only_in_eval.region_limited.vcf"
        limit_cmd = [BEDTOOLS_LOCATION, "intersect", "-wa", "-u", "-header", "-a", extra_sites_file, "-b", bed_file, ">", extra_sites_region_limited_file]
        logger.debug("Running:%s" % " ".join(limit_cmd))
        subprocess.call(" ".join(limit_cmd), shell=True)
        missed_sites_in_region = count_sites(missed_sites_region_limited_file)
        extra_sites_in_region = count_sites(extra_sites_region_limited_file)
        print "[POSITION_BASED]Missed Sites (Limited by supplied bed file): %s" % missed_sites_in_region
        print "[POSITION_BASED]Extra Sites (Limited by supplied bed file): %s" % extra_sites_in_region



def calculate_average_mapping_quality(per_base_counts):
    total_depth = 0
    total_mapping_quality = 0
    for basecount in per_base_counts:
        #base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end
        fields = basecount.split(":")
        total_depth+=int(fields[1])
        total_mapping_quality=(int(fields[1])*float(fields[2]))
    if total_depth ==0:
        return 0
    return total_mapping_quality/total_depth

def xlabel_ticks(max_value, num_ticks):
    #    tick_spacing = max_value / num_ticks 
    #tick_labels = list()
    #tick_positions = list()
    #current_tick_position=0
    #for x in range(num_ticks):
    #    current_tick_position+=tick_spacing
    #    tick_positions.append(current_tick_position)
    #    tick_labels.append(str(current_tick_position))
    tick_positions = [0, 10, 100, 1000]
    tick_labels = ["0", "10", "100", "1000"]
    return (tick_positions, tick_labels)

def calculate_color(fields):
    ref_base = fields[2].upper()
    ref_count=0
    non_ref_count=0
    for basecount in fields[-5:]:
        basecount_fields = basecount.split(":")
        if (basecount_fields[0] == ref_base) or (basecount_fields[0]=="="):
            #reference base
            ref_count+=int(basecount_fields[1]) #add coverage to total
        else:
            non_ref_count+=int(basecount_fields[1]) 
    ref_ratio = float(ref_count) / (non_ref_count + ref_count) #for red yellow green, we want 100% ref points to be "green" / good
    return plt.cm.RdYlGn(ref_ratio)


def generate_coverage_and_quality_graph(bam_readcount_output, mapping_qual_cutoff=20,coverage_cutoff=15):
    """ this function relies on output of the Wash. U program bam-readcount to produce a matplotlib PNG 
    that describes variants not called by the VCF we are evaluating, but present in the reference"""
    fh = open(bam_readcount_output, "r")
    x=list()
    y=list()
    colors = list()
    max_coverage = 0
    max_quality = 0
    passing_sites = 0
    while(True):
        line = fh.readline()
        if not line:
            break
        fields = line.split()
        coverage = fields[3] # bam readcount format: chr pos base depth [per base counts]
        if(max_coverage < int(coverage)):
            max_coverage = int(coverage)
        average_mapping_quality = calculate_average_mapping_quality(fields[-5:])
        color_value = calculate_color(fields)
        if(average_mapping_quality > max_quality):
            max_quality=average_mapping_quality
        if(average_mapping_quality > mapping_qual_cutoff and int(coverage) > coverage_cutoff):
            passing_sites+=1
        x.append(int(coverage))
        y.append(average_mapping_quality)
        colors.append(color_value)
    fig = plt.figure(figsize=(10,7.5), dpi=80)
    plt.scatter(x, y, color=colors, marker=".")
    plt.axhline(y=20, color="red")
    plt.axvline(x=15, color="red")
    plt.xscale('log')
    (axis_positions, axis_labels) = xlabel_ticks(max_coverage+10, 4)
    plt.xticks(axis_positions, axis_labels)
    plt.xlim(0, max_coverage+10)
    plt.ylim(0, max_quality+5)
    plt.xlabel("Coverage (log scale)")
    plt.ylabel("Average Mapping Quality")
    plt.suptitle("Gold SNP missed sites: Mapping Quality vs Coverage")
    plt.text(15, max_quality+2, "Coverage Cutoff:%s" % coverage_cutoff, color="red", fontsize="10")
    plt.text(max_coverage-100, 21, "Quality Cutoff:%s" % mapping_qual_cutoff, color="red", fontsize="10", horizontalalignment="right")
    ax = plt.gca()
    plt.text(.8, .8, "HQ sites missed:%s" % passing_sites, horizontalalignment="right", color="red", transform=ax.transAxes)
    #plt.show()
    scalar_colormap = plt.cm.ScalarMappable(cmap="RdYlGn")
    scalar_colormap.set_array(colors)
   # axins = inset_axes(ax, width ="2.5%", height="50%", loc=3, bbox_to_anchor=(1.02, 0., 1,1), bbox_transform=ax.transAxes, borderpad=0)
    axins = inset_axes(ax, width ="2.5%", height="30%", loc=2)
    cbar = fig.colorbar(scalar_colormap, cax=axins, ticks=[0, .5,1])
    cbar.ax.set_yticklabels(["100% Variant", "50% Ref", "100% Ref"])
    for label in cbar.ax.get_yticklabels():
        label.set_fontsize(8)
    plt.savefig("test.png")




def prepare_sites_file_from_vcf(missed_sites_file):
    """ bam readcount takes a sub-optimal input format you have to construct to supply to it """
    fh = open(missed_sites_file, "r")
    ofh = open("bam_readcount.input_sites", "w")
    while(True):
       line = fh.readline()
       if not line:
           break
       fields = line.split()
       #CHR POS POS for bam readcount, pulling from a vcf
       ofh.write("\t".join([fields[0], fields[1], fields[1]]) + "\n")
    ofh.close()
    return "bam_readcount.input_sites"

def count_sites(vcf_file):
    """ count lines in vcf file """
    fh = open(vcf_file, "rb")
    line_count = 0
    while(1):
        line = fh.readline()
        if not line:
            return line_count
        else:
            if line[0]!='#':
                line_count+=1

def tabix_file(vcf):
    """ index a vcf file with tabix for random access"""
    cmd = [TABIX_LOCATION, "-p" , "vcf", vcf]
    logger.debug("Tabix command: %s" % " ".join(cmd))
    rv = subprocess.call(cmd)

def get_tabix_chrom(vcf, chromosome, output=None):
    """ ask tabix to give us one chromosome in a file to work with """
    if not output:
        output = chromosome + "." + os.path.basename(vcf)
    cmd = [TABIX_LOCATION, "-h", vcf, chromosome, ">", output]
    logger.debug("Tabix command: %s" % " ".join(cmd))
    subprocess.call(" ".join(cmd), shell=True)
    return output

def bedtools_loj(ref_vcf, eval_vcf, output=None):
    """ use left outer join to produce a combined-vcf format including reference vcfs lines that 
    do not intersect a position in the vcf to be evaluated """
    if not output:
        output="left_outer_join_temp"
    cmd = [BEDTOOLS_LOCATION, "intersect", "-loj", "-f 1", "-a", ref_vcf, "-b", eval_vcf, ">" , output]
    logger.debug("Bedtools command: %s" % " ".join(cmd))
    subprocess.call(" ".join(cmd), shell=True)
    return output

def make_vcf_dict(vcf_line):
    """ convenience function to make a python dictionary from a VCF line
    also populates a "GT" key in dict for even more convenience"""
    keys = ["CHR", "POS", "ID", "REF", "ALT" , "SCORE", "FILTER", "INFO", "FORMAT", "SAMPLE"] 
    vcf=dict(zip(keys, vcf_line))
    sample_fields = vcf["SAMPLE"].split(":")
    vcf["GT"]=sample_fields[0]
    self.debug("VCF Dict Created: %s" % vcf)
    return vcf 

def determine_mutation_type(vcf):
    """ categorize mutation we are evaluating according to refence and alt lengths as SNP, indel, or BOTH """
    ref_allele = vcf["REF"]
    alt_allele = vcf["ALT"]
    alt_alleles = alt_allele.split(",")
    types = list()
    for alt in alt_alleles:
        if len(ref_allele) == len(alt) and len(ref_allele) ==1:
            types.append("SNP")
        elif len(ref_allele) > 1 or len(alt) > 1:
            types.append("INDEL")
    #set function in python will condense duplicate entries so we just check if both SNP and INDEL are here
    if(len(set(types))==1):
        return types[0]
    else:
        return "MIXED"

def compare_genotype(ref_vcf, eval_vcf):
    """ the meat of the comparator tool.  properly compare GT's across two files
    need to reconstruct the [REF, ALT[0]...ALT[n]] index and turn genotypes back into real alleles 
    then bin the site according to how well the evaluation vcf has done """
    ref_alt = ref_vcf["ALT"].split(",")
    ref_gt = re.split("[/|]", ref_vcf['GT'])
    #create list to use ref GT alleles and indices into
    ref_alt.insert(0,ref_vcf['REF'])
    #use the GT indices as array indices as VCF intends - Ref = index 0, all alts = indices 1..n
    ref_genotype = [ref_alt[int(ref_gt[0])],ref_alt[int(ref_gt[1])]]
    eval_alt = eval_vcf["ALT"].split(",")
    eval_gt = re.split("[/|]", eval_vcf['GT'])
    #create list to use eval GT alleles as indices into
    eval_alt.insert(0, eval_vcf['REF'])
    #use the GT indices as array indices as VCF intends - Ref = index 0, all alts = indices 1..n
    eval_genotype = [eval_alt[int(eval_gt[0])],eval_alt[int(eval_gt[1])]]
    #mset DOES allow sets to have non-unique members e.g. '["A","A"]' will remain uncondensed
    #so taking the intersection of two of these i
    ##will tell us how many of the alleles we got right
    shared_allele_counts = mset(ref_genotype) & mset(eval_genotype)
    number_of_alleles_matched = sum(shared_allele_counts.values())
    logger.debug("computed genotypes: %s, %s, shared_allele_count: %s, number_of_alleles_matched: %s" % (ref_genotype, eval_genotype, shared_allele_counts, number_of_alleles_matched))
    return number_of_alleles_matched 






def compare_vcf_lines(ref_line, eval_line):
    """based on the output of our core comparison function, classify the correctness of a call
    using a few different bins"""

    #CHROM, POS, ID, REF, ALT, SCORE, FILTER, INFO, FORMAT, SAMPLE
    ref_vcf = make_vcf_dict(ref_line)
    mutation_type = determine_mutation_type(ref_vcf)
    if eval_line[0] == ".":
        return ("NOT_FOUND_IN_EVALUATION_CALLSET", mutation_type, None)
    eval_vcf = make_vcf_dict(eval_line)
    number_of_alleles_matched = compare_genotype(ref_vcf, eval_vcf)
    if number_of_alleles_matched==2:
        #are genotype an exact match
        return ("EXACT_MATCH", mutation_type, None)
    elif number_of_alleles_matched==1:
        return  ("ONE_ALLELE_CORRECT", mutation_type, [ref_vcf["GT"], eval_vcf["GT"]])
    elif number_of_alleles_matched==0:
        return ("POSITONAL_MATCH_WRONG_ALLELES", mutation_type, [ref_vcf["ALT"], eval_vcf["ALT"]])


def parse_bedtools_intersection(loj_bedtools_file):
    """after a bedtools loj on 2 single sample vcf files, we should have a 20 column file. here we break
    it apart into two vcf dicts line-by-line for easy comparison from this function on """
    #expecting a  left outer joined bedtools file here 
    #should be 20 columns for two one sample vcf files
    fh = open_file(loj_bedtools_file)
    stats_dict = dict()
    prev_pos = None
    prev_line = None
    doublings=0
    while(1):
        line = fh.readline()
        if not line:
            break
        fields = line.split()
        if(prev_pos == fields[1]):
            doublings+=1 
        if len(fields) > 20: 
            print >>sys.stderr, "Expected <=20 fields when joining two single sample vcf files, problem line: %s" % line
            sys.exit(1)
        ref_vcf_line = fields[0:10]
        eval_vcf_line = fields[10:20]
        prev_pos = fields[1]
        prev_line = line
        (status, mutation_type, differences) = compare_vcf_lines(ref_vcf_line, eval_vcf_line)
        if mutation_type not in stats_dict:
            stats_dict[mutation_type]=dict()
        if status not in stats_dict[mutation_type]:
            stats_dict[mutation_type][status]=1
        else:
            stats_dict[mutation_type][status]+=1
      #  if(differences):
            # print >>sys.stderr, differences
    return (stats_dict, doublings)

def open_file(filename):
    """convenience method to handle gz or not transparently based on magic number of file"""
    with magic.Magic(flags=magic.MAGIC_MIME_TYPE) as m:
        if(m.id_filename(filename).find("gzip") != -1):
            fh = gzip.open(filename, "rb")
        else:
            fh = open(filename, "rb")
        return fh


def determine_chrom_list(ref_vcf, eval_vcf):
    """check to determine at least the big 24 chromosomes overlap 
    set chromosome intersection list and most importantly
    give the user an idea if the references both include unplaced contigs or not"""
    chrom_dict = dict()
    for vcf in [ref_vcf, eval_vcf]:
        fh = open_file(vcf)
        while(1):
            line = fh.readline()
            if not line:
                break
            if line[0]!="#":
                fields = line.split()
                if fields[0] not in chrom_dict:
                    chrom_dict[fields[0]]=vcf
                elif fields[0] == vcf:
                    True
                else:
                    chrom_dict[fields[0]]=2
    shared_chroms = list()
    for (key, value) in chrom_dict.items():
        if(value ==2):
            shared_chroms.append(key)
    print >>sys.stderr, "Found %d chromosomes in common between the two files: %s" % ( len(shared_chroms), ", ".join(shared_chroms))
    return shared_chroms

def pretty_print_stats(stats, chrom=None):
    """take an stat dict (example in comments in code) and display it in an excel/user friendly way"""
    # {'MIXED': {'DIFFERENT ALTS': 214, 'NOT_FOUND_IN_EVALUATION_CALLSET': 1174}, 'INDEL': {'DIFFERENT ALTS': 2874, 'EXACT_MATCH': 1147, 'DIFFERENT_GENOTYPES': 208, 'NOT_FOUND_IN_EVALUATION_CALLSET': 169375}, 'SNP': {'DIFFERENT ALTS': 107, 'EXACT_MATCH': 2715136, 'DIFFERENT_GENOTYPES': 4688, 'NOT_FOUND_IN_EVALUATION_CALLSET': 21429}}
    totals = dict()
    for k in stats.iterkeys():
        for j in stats[k].iterkeys():
            if k not in totals:
                totals[k]=stats[k][j]
            else:
                totals[k]+=stats[k][j]
    if(chrom):
        print "Chromosome %s stats:" % chrom
    total_sites=0
    total_exact_match_sites=0

    for mut_type in stats.iterkeys():
        for stat in stats[mut_type].iterkeys():
            print "%s\t%s\t%d"%(mut_type, stat, stats[mut_type][stat])
        print "%s\t%s\t%d"%(mut_type, "TOTAL", totals[mut_type])
        total_sites+=totals[mut_type]
        if "EXACT_MATCH" in stats[mut_type]:
            correct_percent = float(stats[mut_type]["EXACT_MATCH"]) / totals[mut_type] * 100
            total_exact_match_sites+=stats[mut_type]["EXACT_MATCH"]
        else:
            correct_percent=0
        if "NOT_FOUND_IN_EVALUATION_CALLSET" in stats[mut_type]:
            missed_percent = float(stats[mut_type]["NOT_FOUND_IN_EVALUATION_CALLSET"]) / totals[mut_type] * 100
        else:
            missed_percent = 0
        print "%s\t%s\t%s"%(mut_type, "PERCENTAGE_EXACTLY_CORRECT", "%0.2f" % correct_percent)
        print "%s\t%s\t%s"%(mut_type, "PERCENTAGE_COMPLETELY_MISSED", "%0.2f" % missed_percent)
    print "%s\t%s\t%s"%("OVERALL", "PERCENTAGE_EXACTLY_CORRECT", "%0.2f" % (total_exact_match_sites/total_sites * 100))
    print "%s\t%s\t%s"%("OVERALL", "EXACT_MATCH", total_exact_match_sites)
    print "%s\t%s\t%s"%("OVERALL", "TOTAL_SITES", total_sites)


def deep_compare(ref_vcf, eval_vcf):
    """ do a thorough position and allele comparison"""
    chrom_list = determine_chrom_list(ref_vcf, eval_vcf)
    #chrom_list = ["1"]
    tabix_file(ref_vcf)
    tabix_file(eval_vcf)
    total_stats = None 
    total_doublings = 0
    for chromosome in chrom_list:
        ref_chrom = get_tabix_chrom(ref_vcf, chromosome)
        eval_chrom = get_tabix_chrom(eval_vcf, chromosome)
        intersected_output = bedtools_loj(ref_chrom, eval_chrom)
        (stats_dict, doublings) = parse_bedtools_intersection(intersected_output)
        total_doublings+=doublings
       # print >>sys.stderr, "%d doublings in chrom %s" % (doublings, chromosome)
        os.unlink(ref_chrom)
        os.unlink(eval_chrom)
        os.unlink(intersected_output)
      #  pretty_print_stats(stats_dict, chrom=chromosome)
        if not total_stats:
            total_stats=stats_dict
        else:
            for k in stats_dict.iterkeys():
                for j in stats_dict[k].iterkeys():
                    if not j in total_stats[k]:
                        total_stats[k][j]=0
                    total_stats[k][j]+=stats_dict[k][j]
    print "Overall Stats:"
    print "Total Doubled lines(positions that appear twice in the evaluation vcf file potentially giving conflicting genotypes): %d" % total_doublings
    pretty_print_stats(total_stats)
     #indels_by_position, indels_by_position_and_genotype
        #snps_by_position, snps_by_position_and_genotype

def main(bed_file, bam_file, ref_snp_vcf, eval_snp_vcf, ref_fasta, do_deep_compare, do_gatk, do_positional, log_level):
    #run VT on eval file to make sure any indels are left shifted
    #FIXME pass in other logging levels
    global logger
    logger = configure_logging(log_level)
    normalized_vcf_output = "temp.normalized.vcf.gz"
    cmd = [ VT_LOCATION, "normalize", "-r", ref_fasta, eval_snp_vcf, "-o", normalized_vcf_output ]
    logger.info("Running:" + " ".join(cmd))
    subprocess.call(cmd)
    #run GATK evaluator and store output
    if(do_gatk==True):
        logger.info("GATK flag set... Running GATK GenotypeConcordance")
        gatk_results = "gatk_genotype_concordance_output"
        cmd = [JAVA_LOCATION, "-jar", GATK_LOCATION, "-T", "GenotypeConcordance", "--comp", ref_snp_vcf, "--eval", normalized_vcf_output, "-R", ref_fasta, "-o", gatk_results]
        if(bed_file != None):
            logger.info("Restricting GATK comparison to bed file region: %s" % bed_file)
            cmd.extend(["-L", bed_file])
        logger.debug("GATK COMMAND: %s" % " ".join(cmd))
        subprocess.call(cmd)
    #run bed intersect by position and get the disjoint set only in reference
    if(do_positional==True):
        rough_positional_intersect(ref_snp_vcf, normalized_vcf_output, bed_file, );

    if(do_deep_compare==True):
        #FIXME USE NORMALIZED VCF WHEN DONE TESTING
        deep_compare(ref_snp_vcf, normalized_vcf_output)
   #prepare sites into bam-readcount format CHR START STOP
    if(bam_file != None and os.path.exists(bam_file)):
        #FIXME: need to have a missed_sites_file for sure even if user skips positional
        bam_readcount_sites_file = prepare_sites_file_from_vcf(missed_sites_file)
        #run bam-readcount
        bam_readcount_output = "bam_readcount.output"
        cmd = [BAM_READCOUNT_LOCATION, "-w 1", "-l", bam_readcount_sites_file, "-f", ref_fasta, bam_file, ">", bam_readcount_output]
        print " ".join(cmd)
        subprocess.call(" ".join(cmd), shell=True)
        #generate graphs
        missing_sites_coverage_quality_png = generate_coverage_and_quality_graph(bam_readcount_output)
    else:
        print >>sys.stderr, "BAM file not supplied - skipping missing site interrogation!"


def configure_logging(log_level):
    logger = logging.getLogger("QIAGEN VCF Accuracy Eval")
    if log_level=="DEBUG":
        log_level=logging.DEBUG
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    return logger

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take in a high quality snp/indel set and evaluate a call set for concordance")
    parser.add_argument("--bed_file", action="store", dest="bed_file", default=None, help="bed file to limit comparisons to certain regions- optional")
    parser.add_argument("--bam_file", action="store", dest="bam_file", help="bam file to use to get readcounts for missing sites")
    parser.add_argument("--log_level", action="store", dest="log_level", default=logging.INFO, help="INFO for basic, DEBUG for (very) detailed debug output)")
    parser.add_argument("--ref_vcf", action="store", dest="ref_vcf", help="file of SNPs or INDELs you know to be true in your callset", required=True)
    parser.add_argument("--eval_vcf", action="store", dest="eval_vcf", help="file of SNPs and/or INDELs you want to compare to the known true callset", required=True)
    parser.add_argument("--ref_fasta", action="store", dest="ref_fasta", help="reference fasta used to call the bam", required=True)
    parser.add_argument("--deep-compare", action="store_true", dest="deep_compare", default="False", help="break things down by chromsome, snp, and indel")
    parser.add_argument("--gatk", action="store_true", dest="gatk", default="False", help="run GATK GenotypeConcordance for a further check")
    parser.add_argument("--no-rough", action="store_false", dest="rough", default="True", help="skip positional/missing site bed file creation")
    args=parser.parse_args()
    if not args.deep_compare and not args.gatk:
        print >>sys.stderr, "Only doing fast positional intersect, aka --rough"
    #FIXME add nice error messages if not enough files are supplied
    print args.gatk, args.deep_compare, args.rough
    
    main(args.bed_file,args.bam_file, args.ref_vcf, args.eval_vcf, args.ref_fasta, args.deep_compare, args.gatk, args.rough, args.log_level);