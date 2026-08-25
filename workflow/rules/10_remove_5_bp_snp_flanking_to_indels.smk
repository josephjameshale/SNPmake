# remove snps with 5bp of an indel 
rule flag_indel_prox_snps:
    input:
        vcf_filter = "results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered.vcf.gz",
    output:
        vcf_indelprox_temp=temp("results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf"),
        masked_positions="results/{prefix}/filtered_vcf/{sample}/{sample}_indelprox_masked_positions.txt",
    params:
        indel_prox_window = config["indel_prox_window"],
    benchmark: 
        "benchmarks/{prefix}/remove_5bp_snp_flanking_to_indels/{sample}.benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60
    wrapper:
        "file:workflow/wrapper_functions/remove_5_bp_snp_flanking_to_indels"


rule bgzip_and_index_indelprox_vcf:
    input:
        vcf_indelprox_temp="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf",
    output:
        vcf_indelprox="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf.gz",
        vcf_indexprox_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30
    shell:
        """
        bgzip -f -c {input.vcf_indelprox_temp} > {output.vcf_indelprox}
        tabix -f -p vcf {output.vcf_indelprox}
        """


