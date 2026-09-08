
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

# generate a bed file to mask the beginning and end of each scaffold
rule make_scaffold_mask_bed:
    input:
        ref=REF_GENOME,
    output:
        scaffold_mask_bed_temp = temp("results/{prefix}/bedtools/{prefix}_scaffold_mask_temp.bed"),
        scaffold_mask_bed = "results/{prefix}/bedtools/{prefix}_scaffold_mask.bed",
    params:
        scaffold_mask_size = config["scaffold_mask_size"],
        feature_buffer_size = config["feature_buffer_size"],
        filter_label = "SCAFFOLD_MASK",
    conda:
        "envs/bioconda_gffutils.yaml"
    shell:
        """
        python3.13 workflow/scripts/make_scaffold_mask_bed.py --input {input.ref} --output {output.scaffold_mask_bed_temp} \
        --scaffold_mask_size {params.scaffold_mask_size} --feature_buffer_size {params.feature_buffer_size}
        # add a fourth column to the bed file that simply contains the label "SCAFFOLD_MASK"
        awk -v label={params.filter_label} 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, label}}' {output.scaffold_mask_bed_temp} > {output.scaffold_mask_bed}
        """


# use bedtools to sort the scaffold_mask_bed file
rule combine_scaffold_mask_and_lowcoverage_bed:
    input:
        scaffold_mask_bed = "results/{prefix}/bedtools/{prefix}_scaffold_mask.bed",
        lowcoverage_bed_sorted = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_temp_sorted.bed",
    output:
        final_mask_bed = "results/{prefix}/bedtools/{sample}/{sample}_final_mask.bed",
    singularity:
        "docker://staphb/bedtools:2.31.1"
    shell:
        """
        cat {input.scaffold_mask_bed} {input.lowcoverage_bed_sorted} | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct -delim ';' > {output.final_mask_bed}
        """


rule bedtools_index:
    input:
        final_mask_bed = "results/{prefix}/bedtools/{sample}/{sample}_final_mask.bed",
    output:
        final_mask_bed_gz = "results/{prefix}/bedtools/{sample}/{sample}_final_mask.bed.gz",
        final_mask_bed_gz_tbi = "results/{prefix}/bedtools/{sample}/{sample}_final_mask.bed.gz.tbi",
        final_mask_header = "results/{prefix}/bedtools/{sample}/{sample}_final_mask_header.hdr",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    shell:
        """
        bgzip -c {input.final_mask_bed} > {output.final_mask_bed_gz}
        tabix -f -p bed {output.final_mask_bed_gz}
        printf '%s\n' '##FILTER=<ID=FAIL_LOW_COVERAGE,Description="Variant is in a low-coverage region">' > {output.final_mask_header}
        printf '%s\n' '##FILTER=<ID=SCAFFOLD_MASK,Description="Variant is in a scaffold masking region">' >> {output.final_mask_header}
        """


