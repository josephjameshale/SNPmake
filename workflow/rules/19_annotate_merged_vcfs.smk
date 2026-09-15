
rule annotate_snp_vcf:
    input:
        merged_snp_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz",
        merged_snp_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz.tbi",
        snpeff_config = config["snpeff_config"],
    output:
        vcf_snp_annotated = "results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_annotated.vcf.gz",
        snp_annotation_stats = "results/{prefix}/merged_vcf/{prefix}_snp_annotation_stats.html",
    singularity:
        "docker://staphb/snpeff:5.4c"
    params:
        ref_genome_id = config["reference_name"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        set -euo pipefail
        snpEff ann -c {input.snpeff_config} -v -s {output.snp_annotation_stats} {params.ref_genome_id} {input.merged_snp_vcf} | bgzip -c > {output.vcf_snp_annotated}
        """


rule annotate_indel_vcf:
    input:
        merged_indel_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only.vcf.gz",
        merged_indel_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only.vcf.gz.tbi",
        snpeff_config = config["snpeff_config"],
    output:
        vcf_indel_annotated = "results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only_annotated.vcf.gz",
        indel_annotation_stats = "results/{prefix}/merged_vcf/{prefix}_indel_annotation_stats.html",
    singularity:
        "docker://staphb/snpeff:5.4c"
    params:
        ref_genome_id = config["reference_name"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        set -euo pipefail
        snpEff ann -c {input.snpeff_config} -v -s {output.indel_annotation_stats} {params.ref_genome_id} {input.merged_indel_vcf} | bgzip -c > {output.vcf_indel_annotated}
        """


