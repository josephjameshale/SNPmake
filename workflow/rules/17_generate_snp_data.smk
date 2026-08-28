
# use the vcf file with all pass and fail records to build a dataframe of all positions 
# note the the command below only includes SNPs!!
rule snp_data_isolate:
    input:
        vcf_indelprox_lowcov="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz",
        vcf_indexprox_lowcov_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz.tbi",
    output:
        snp_data_isolate="results/{prefix}/consensus/{sample}/{sample}_snp_data.tsv",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        r"""
        {
            printf 'CHROM\tPOS\t%s\n' {wildcards.sample}

            bcftools query --include 'TYPE="snp" && GT="alt"' --format '%CHROM\t%POS\t%FILTER\n' {input.vcf_indelprox_lowcov}
        } > {output.snp_data_isolate}
        """


# after all snp_data files have been generated, merge them into a single file with a python script
rule merge_snp_data:
    input:
        snp_data_isolate=expand(
            "results/{prefix}/consensus/{sample}/{sample}_snp_data.tsv",
            sample=SAMPLES, prefix=PREFIX
        ),
    output:
        merged_snp_data="results/{prefix}/consensus/{prefix}_merged_snp_data.tsv",
        fail_bed="results/{prefix}/consensus/{prefix}_fail.bed",
    params:
        fail_threshold = config["fail_threshold"],
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        """
        python workflow/scripts/merge_snp_data.py --input {input.snp_data_isolate} --output {output.merged_snp_data} \
        --fail_threshold {params.fail_threshold} --fail_bed {output.fail_bed}
        """

# use bedtools to sort the fail_bed file
rule bedtools_sort_fail_bed:
    input:
        fail_bed="results/{prefix}/consensus/{prefix}_fail.bed",
        ref_genome_fai = config["reference_genome"] + ".fai",
    output:
        fail_bed_sorted="results/{prefix}/consensus/{prefix}_fail_sorted.bed",
    singularity:
        "docker://staphb/bedtools:2.31.1"
    shell:
        """
        bedtools sort -i {input.fail_bed} -faidx {input.ref_genome_fai} > {output.fail_bed_sorted}
        """


# remove the positions that overlap with the fail_bed file from the merged vcf file
rule mask_fail_positions:
    input:
        merged_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz",
        merged_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz.tbi",
        fail_bed_sorted="results/{prefix}/consensus/{prefix}_fail_sorted.bed",
    output:
        merged_filtered_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_filtered.vcf.gz",
        merged_filtered_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_filtered.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        """
        bcftools view -T ^{input.fail_bed_sorted} -Oz -o {output.merged_filtered_vcf} {input.merged_vcf}
        bcftools index -f -t {output.merged_filtered_vcf}
        """



