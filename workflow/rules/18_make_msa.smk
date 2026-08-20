
# merge the final consensus vcf files for each sample
rule merge_vcf:
    input:
        ref=REF_GENOME,
        vcfs=expand(
            "results/{prefix}/consensus/{sample}/{sample}_pass_snp_only.vcf.gz",
            sample=SAMPLES, prefix=PREFIX
        ),
        vcf_tbis=expand(
            "results/{prefix}/consensus/{sample}/{sample}_pass_snp_only.vcf.gz.tbi",
            sample=SAMPLES, prefix=PREFIX
        ),
    output:
        merged_vcf="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz",
        merged_vcf_tbi="results/{prefix}/merged_vcf/{prefix}_merged_pass_snp_only.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=180
    shell:
        """
        bcftools merge -Oz -o {output.merged_vcf} {input.vcfs}
        bcftools index -f -t {output.merged_vcf}
        """




# original method to make msa files - probably not applicable to a reference with multiple contigs
# rule make_msa:
#     input:
#         ref=REF_GENOME,
#         cons=expand(
#             "results/{prefix}/consensus/{sample}/{sample}_consensus.fa",
#             sample=SAMPLES, prefix=PREFIX
#         )
#     output:
#         msa="results/{prefix}/alignment/{prefix}_consensus_msa.fa"
#     threads: 1
#     resources:
#         mem_mb=1000,
#         runtime=180
#     shell:
#         """
#         cat {input.ref} {input.cons} > {output.msa}
#         """
