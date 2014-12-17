#!/usr/bin/env python

import sys,re,os, subprocess
import magic, gzip
import argparse 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

JAVA_LOCATION = "/Library/Java/JavaVirtualMachines/jdk1.7.0_60.jdk/Contents/Home/bin/java"
GATK_LOCATION =  "../GATK/GenomeAnalysisTK.jar"
BAM_READCOUNT_LOCATION = "../bam-readcount/bin/bam-readcount"
VT_LOCATION = "../vt/vt"
BEDTOOLS_LOCATION = "../bedtools/bin/bedtools"
TABIX_LOCATION = "tabix"
def calculate_average_mapping_quality(per_base_counts):
    print per_base_counts
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
    cmd = [TABIX_LOCATION, "-p" , "vcf", vcf]
    print >>sys.stderr, " ".join(cmd)
    rv = subprocess.call(cmd)
   
def get_tabix_chrom(vcf, chromosome, output=None):
    if not output:
        output = chromosome + "." + os.path.basename(vcf)
    cmd = [TABIX_LOCATION, "-h", vcf, chromosome, ">", output]
 ###   print >>sys.stderr, " ".join(cmd)
    subprocess.call(" ".join(cmd), shell=True)
    return output

def bedtools_loj(ref_vcf, eval_vcf, output=None):
    if not output:
        output="left_outer_join_temp"
    cmd = [BEDTOOLS_LOCATION, "intersect", "-loj", "-f 1", "-a", ref_vcf, "-b", eval_vcf, ">" , output]
   # print >>sys.stderr, " ".join(cmd)
    subprocess.call(" ".join(cmd), shell=True)
    return output

def make_vcf_dict(vcf_line):
    keys = ["CHR", "POS", "ID", "REF", "ALT" , "SCORE", "FILTER", "INFO", "FORMAT", "SAMPLE"] 
    vcf=dict(zip(keys, vcf_line))
    sample_fields = vcf["SAMPLE"].split(":")
    vcf["GT"]=sample_fields[0]
    return vcf 

def determine_mutation_type(vcf):
    ref_allele = vcf["REF"]
    alt_allele = vcf["ALT"]
    alt_alleles = alt_allele.split(",")
    types = list()
    for alt in alt_alleles:
        if len(ref_allele) == len(alt) and len(ref_allele) ==1:
            types.append("SNP")
        elif len(ref_allele) > 1 or len(alt) > 1:
            types.append("INDEL")
    if(len(set(types))==1):
        return types[0]
    else:
        return "MIXED"
            
def compare_vcf_lines(ref_line, eval_line):
    #CHROM, POS, ID, REF, ALT, SCORE, FILTER, INFO, FORMAT, SAMPLE
    ref_vcf = make_vcf_dict(ref_line)
    mutation_type = determine_mutation_type(ref_vcf)
    if eval_line[0] == ".":
        return ("NOT_FOUND_IN_EVALUATION_CALLSET", mutation_type, None)
    eval_vcf = make_vcf_dict(eval_line)
    if(ref_vcf["ALT"] == eval_vcf["ALT"]):
        #are genotype an exact match
        if cmp(eval_vcf["GT"].split("[/|]"),ref_vcf["GT"].split("[/|]"))==0:
             return ("EXACT_MATCH", mutation_type, None)
        else:
            return  ("DIFFERENT_GENOTYPES", mutation_type, [ref_vcf["GT"], eval_vcf["GT"]])
    else:
        return ("DIFFERENT ALTS", mutation_type, [ref_vcf["ALT"], eval_vcf["ALT"]])
        

def parse_bedtools_intersection(loj_bedtools_file):
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
        if len(fields) != 20: 
            print >>sys.stderr, "Expected 20 fields when joining two single sample vcf files, problem line: %s" % line
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
     with magic.Magic(flags=magic.MAGIC_MIME_TYPE) as m:
        if(m.id_filename(filename).find("gzip") != -1):
            fh = gzip.open(filename, "rb")
        else:
            fh = open(filename, "rb")
     return fh


def determine_chrom_list(ref_vcf, eval_vcf):
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
    for mut_type in stats.iterkeys():
        for stat in stats[mut_type].iterkeys():
            print "%s\t%s\t%d"%(mut_type, stat, stats[mut_type][stat])
        print "%s\t%s\t%d"%(mut_type, "TOTAL", totals[mut_type])
        if "EXACT_MATCH" in stats[mut_type]:
            correct_percent = float(stats[mut_type]["EXACT_MATCH"]) / totals[mut_type] * 100
        else:
            correct_percent=0
        if "NOT_FOUND_IN_EVALUATION_CALLSET" in stats[mut_type]:
            missed_percent = float(stats[mut_type]["NOT_FOUND_IN_EVALUATION_CALLSET"]) / totals[mut_type] * 100
        else:
            missed_percent = 0
        print "%s\t%s\t%s"%(mut_type, "PERCENTAGE_EXACTLY_CORRECT", "%0.2f" % correct_percent)
        print "%s\t%s\t%s"%(mut_type, "PERCENTAGE_COMPLETELY_MISSED", "%0.2f" % missed_percent)
    

def deep_compare(ref_vcf, eval_vcf):
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
        print >>sys.stderr, "%d doublings in chrom %s" % (doublings, chromosome)
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
    print "Total Doubled lines: %d" % total_doublings
    pretty_print_stats(total_stats)
     #indels_by_position, indels_by_position_and_genotype
        #snps_by_position, snps_by_position_and_genotype

def main(bed_file, bam_file, ref_snp_vcf, eval_snp_vcf, ref_fasta, do_deep_compare, do_gatk, do_positional):
    #run VT on eval file to make sure any indels are left shifted
    normalized_vcf_output = "temp.normalized.vcf.gz"
    cmd = [ VT_LOCATION, "normalize", "-r", ref_fasta, eval_snp_vcf, "-o", normalized_vcf_output ]
    print " ".join(cmd)
    subprocess.call(cmd)
    #run GATK evaluator and store output
    if(do_gatk==True):
        gatk_results = "gatk_genotype_concordance_output"
        cmd = [JAVA_LOCATION, "-jar", GATK_LOCATION, "-T", "GenotypeConcordance", "--comp", ref_snp_vcf, "--eval", normalized_vcf_output, "-R", ref_fasta, "-o", gatk_results]
        if(bed_file != None):
            cmd.extend(["-L", bed_file])
        #print " ".join(cmd)
        subprocess.call(cmd)
    #run bed intersect by position and get the disjoint set only in reference
    if(do_positional==True):
        missed_sites_file = "sites_only_in_ref.vcf"
        cmd = [BEDTOOLS_LOCATION, "intersect", "-v", "-header", "-a", ref_snp_vcf, "-b", eval_snp_vcf, ">", missed_sites_file]
        #print " ".join(cmd)
        subprocess.call(" ".join(cmd), shell=True)
        extra_sites_file = "sites_only_in_eval.vcf"
        cmd = [BEDTOOLS_LOCATION, "intersect", "-v", "-header", "-a", eval_snp_vcf, "-b", ref_snp_vcf, ">", extra_sites_file]
        #print " ".join(cmd)
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
            subprocess.call(" ".join(limit_cmd), shell=True)
            extra_sites_region_limited_file="sites_only_in_eval.region_limited.vcf"
            limit_cmd = [BEDTOOLS_LOCATION, "intersect", "-wa", "-u", "-header", "-a", extra_sites_file, "-b", bed_file, ">", extra_sites_region_limited_file]
            subprocess.call(" ".join(limit_cmd), shell=True)
            missed_sites_in_region = count_sites(missed_sites_region_limited_file)
            extra_sites_in_region = count_sites(extra_sites_region_limited_file)
            print "[POSITION_BASED]Missed Sites (Limited by supplied bed file): %s" % missed_sites_in_region
            print "[POSITION_BASED]Extra Sites (Limited by supplied bed file): %s" % extra_sites_in_region

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take in a high quality snp/indel set and evaluate a call set for concordance")
    parser.add_argument("--bed_file", action="store", dest="bed_file", default=None, help="bed file to limit comparisons to certain regions- optional")
    parser.add_argument("--bam_file", action="store", dest="bam_file", help="bam file to use to get readcounts for missing sites")
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
    main(args.bed_file,args.bam_file, args.ref_vcf, args.eval_vcf, args.ref_fasta, args.deep_compare, args.gatk, args.rough);
