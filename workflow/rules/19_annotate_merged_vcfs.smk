
rule annotate_snp_vcf:
    input:
        merged_snp_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz",
        merged_snp_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz.tbi",
        snpeff_config = config["snpeff_config"],
    output:
        vcf_snp_annotated_temp = temp("results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_annotated.vcf"),
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
        snpEff ann -c {input.snpeff_config} -v -s {output.snp_annotation_stats} {params.ref_genome_id} {input.merged_snp_vcf} > {output.vcf_snp_annotated_temp}
        """

rule index_snp_vcf:
    input:
        vcf_snp_annotated_temp = "results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_annotated.vcf",
    output:
        vcf_snp_annotated = "results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_annotated.vcf.gz",
        vcf_snp_annotated_tbi = "results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only_annotated.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    params:
        ref_genome_id = config["reference_name"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        bgzip -c {input.vcf_snp_annotated_temp} > {output.vcf_snp_annotated}
        bcftools index -f -t {output.vcf_snp_annotated}
        """

rule annotate_indel_vcf:
    input:
        merged_indel_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only.vcf.gz",
        merged_indel_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only.vcf.gz.tbi",
        snpeff_config = config["snpeff_config"],
    output:
        vcf_indel_annotated_temp = temp("results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only_annotated.vcf"),
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
        snpEff ann -c {input.snpeff_config} -v -s {output.indel_annotation_stats} {params.ref_genome_id} {input.merged_indel_vcf} > {output.vcf_indel_annotated_temp}
        """

rule index_indel_vcf:
    input:
        vcf_indel_annotated_temp = "results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only_annotated.vcf",
    output:
        vcf_indel_annotated = "results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only_annotated.vcf.gz",
        vcf_indel_annotated_tbi = "results/{prefix}/merged_vcf/{prefix}_merged_pass_indel_only_annotated.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    params:
        ref_genome_id = config["reference_name"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        bgzip -c {input.vcf_indel_annotated_temp} > {output.vcf_indel_annotated}
        bcftools index -f -t {output.vcf_indel_annotated}
        """

