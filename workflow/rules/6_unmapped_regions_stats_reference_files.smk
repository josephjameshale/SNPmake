
# in an effort to simplify this step, only the bedtools genomecov output is generated here
# this file is eventually used to mask the final consensus output file

# this rule uses bedtools to identify regions in the reference genome that have zero coverage in the aligned reads
# contig_name \t window_start \t window_end \ LOW_COVERAGE
rule bedtools_genomecov:
    input:
        index_sorted_dups_rmvd_bam ="results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam"
    output:
        lowcoverage_bed_temp = temp("results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_temp.bed"),
        lowcoverage_bed_sorted = temp("results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_temp_sorted.bed"),
    params:
        coverage_threshold = config["coverage_masking_threshold"],
        filter_label = "FAIL_LOW_COVERAGE",
    singularity:
        "docker://staphb/bedtools:2.31.1"
    benchmark: 
        "benchmarks/{prefix}/bedtools_coverage/{sample}.benchmark.tsv"
    shell:
        """
        bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk -v threshold={params.coverage_threshold} -v label={params.filter_label} 'BEGIN {{OFS="\t"}} $4 < threshold {{print $1, $2, $3, label}}' > {output.lowcoverage_bed_temp}
        bedtools sort -i {output.lowcoverage_bed_temp} > {output.lowcoverage_bed_sorted}
        """

rule bedtools_index:
    input:
        lowcoverage_bed_sorted = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_temp_sorted.bed",
    output:
        lowcoverage_bed = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz",
        lowcoverage_bed_tbi = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz.tbi",
        lowcoverage_header = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_header.hdr",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    shell:
        """
        bgzip -c {input.lowcoverage_bed_sorted} > {output.lowcoverage_bed}
        tabix -f -p bed {output.lowcoverage_bed}
        printf '%s\n' '##FILTER=<ID=FAIL_LOW_COVERAGE,Description="Variant is in a low-coverage region">' > {output.lowcoverage_header}
        """


























# old version
# # NOTE: this pipeline uses the term 'unmapped regions' to refer to regions in the reference genome that have low or zero coverage in the aligned reads.
# # this is different from the standard definition of unmapped reads, which are positions in the query reads that do not align to the reference genome.

# def bed_intervals_to_positions(bed_interval_file):
#     '''Convert a bed file of intervals into a tsv file of positions.'''
#     positions_array = []
#     with open(bed_interval_file, 'r') as fp:
#         for line in fp:
#             contig_name, interval_start, interval_end, coverage = line.strip().split('\t')
#             lower_index = int(interval_start) + 1
#             upper_index = int(interval_end) + 1
#             for positions in range(lower_index,upper_index):
#                 positions_array.append(contig_name + "\t" + str(positions))
#     positions_file = bed_interval_file.strip('.bed') + "_positions.tsv"
#     with open(positions_file, 'w') as f1:
#         for line in positions_array:
#             _ = f1.write(line + "\n")
#     return positions_file

# # here is the original function for reference, which assumes only one contig is present in the reference genome
# # def parse_bed_file(final_bed_unmapped_file):
# #     unmapped_positions_array = []
# #     with open(final_bed_unmapped_file, 'r') as fp:
# #         for line in fp:
# #             line_array = line.split('\t')
# #             lower_index = int(line_array[1]) + 1
# #             upper_index = int(line_array[2]) + 1
# #             for positions in range(lower_index,upper_index):
# #                 unmapped_positions_array.append(positions)
# #     only_unmapped_positions_file = final_bed_unmapped_file + "_positions"
# #     f1=open(only_unmapped_positions_file, 'w+')
# #     for i in unmapped_positions_array:
# #         p_string = str(i) + "\n"
# #         f1.write(p_string)
# #     return only_unmapped_positions_file


# # determine coverage of bam file       
# rule gatk_coverage_depth_statistics:
#     input:
#         index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
#         intervals = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
#     output:
#         gatk_depthCoverage_summary = "results/{prefix}/stats/{sample}/{sample}_depth_of_coverage.sample_summary"
#     params:
#         outdir = "results/{prefix}/stats/{sample}",
#         ref_genome = config["reference_genome"], 
#         prefix = "{sample}",
#     log:
#         "logs/{prefix}/post_align/{sample}/{sample}_coverage_depth.log"
#     singularity:
#         "docker://broadinstitute/gatk:4.6.2.0" 
#     benchmark: 
#         "benchmarks/{prefix}/gatk_coverage_depth_statistics/{sample}.benchmark.tsv"
#     shell:
#         "gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_sorted_dups_rmvd_bam} \
#         --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 \
#         --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {input.intervals} &> {log}"

# # this rule uses bedtools to identify regions in the reference genome that have zero coverage in the aligned reads
# # contig_name \t window_start \t window_end \t coverage
# # (coverage column is always zero in this output file)
# rule bedtools_extract_coverage:
#     input:
#         index_sorted_dups_rmvd_bam ="results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam"
#     output:
#         unmapped_bed = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed"
#     log:
#        "logs/{prefix}/bedtools_unmapped/{sample}/{sample}.log"
#     singularity:
#         "docker://staphb/bedtools:2.31.1"
#     benchmark: 
#         "benchmarks/{prefix}/bedtools_extract_coverage/{sample}.benchmark.tsv"
#     shell:
#         "bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk '$4==0' > {output.unmapped_bed}"

# # this converts the bed file of intervals with zero coverage into a file of positions with zero coverage
# # originally, this was a one-column file of positions only, which used the assumption that only one contig was present
# # the format is now:
# # contig_name \t position  
# rule parse_bed_file_find_unmapped_regions:
#     input:
#         unmapped_bed = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed"
#     output:
#         #unmapped_bam_positions = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped.bed_positions"
#         unmapped_bam_positions = "results/{prefix}/bedtools_unmapped/{sample}/{sample}_unmapped_positions.tsv"
#     benchmark: 
#         "benchmarks/{prefix}/parse_bed_file_find_unmapped_regions/{sample}.benchmark.tsv"  
#     run:
#         bed_intervals_to_positions(input.unmapped_bed)

# # prepare reference size and window files
# # the .size file contains the size of each contig in the reference genome
# # contig_name \t contig_length
# # the .bed file contains 1000 bp windows for each contig in the reference genome
# # contig_name \t window_start \t window_end 
# rule prepare_reference_windows:
#     output:
#         reference_size_file="results/{prefix}/ref_genome_files/{ref_name}.size",
#         reference_window_file = "results/{prefix}/ref_genome_files/{ref_name}.bed"
#     params:
#         ref_genome = config["reference_genome"]
#     # benchmark: 
#     #     "benchmarks/{prefix}/prepare_reference_windows/benchmark.tsv"
#     wrapper:
#         "file:workflow/wrapper_functions/prepare_reference_files"

# # NOTE: this rule does not appear to be used in the workflow
# # the original command was:
# # "bedtools coverage -abam {input.index_sorted_dups_rmvd_bam} -b {input.reference_window_file} > {output.bedgraph_cov} 2>&1 | tee {log}"
# # the updated command changes the order of the inputs and removes the logging
# rule bedcoverage:
#     input:
#         index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/post_align/{wildcards.sample}/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
#         reference_window_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
#     output:
#         bedgraph_cov = f"results/{{prefix}}/bedtools/{{sample}}/bedgraph_coverage/{{sample}}.bedcov"
#     singularity:
#         "docker://staphb/bedtools:2.31.1"
#     log:
#         "logs/{prefix}/{sample}/bedgraph_coverage/{sample}_bedcov.log"
#     benchmark: 
#         "benchmarks/{prefix}/bedcoverage/{sample}.benchmark.tsv"
#     shell:
#         "bedtools coverage -a {input.reference_window_file} -b {input.index_sorted_dups_rmvd_bam} > {output.bedgraph_cov}"

# # NOTE: this rule does not appear to be used in the workflow either
# rule alignment_stats:
#     input:
#         index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam"
#     output:
#         alignment_stats = "results/{prefix}/stats/{sample}/{sample}_alignment_stats.tsv" 
#     singularity:
#         "docker://staphb/samtools:1.23.1"
#     log:
#         "logs/{prefix}/post_align/{sample}/{sample}_stats.log"
#     benchmark: 
#         "benchmarks/{prefix}/alignment_stats/{sample}.benchmark.tsv"
#     shell:
#         "samtools flagstat {input.index_sorted_dups_rmvd_bam} > {output.alignment_stats} &> {log}" 

