
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


