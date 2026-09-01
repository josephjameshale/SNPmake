
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


# concatenate the contigs in the reference genome
rule concatenate_reference_contigs:
    input:
        ref=REF_GENOME,
    output:
        ref_concat="results/{prefix}/consensus/{prefix}_ref_concat.fasta",
        contig_order="results/{prefix}/consensus/{prefix}_ref_contig_order.tsv",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    params:
        reference_name = config["reference_name"],
    conda:
        "envs/concatenate_contigs.yaml"
    shell:
        """
        python3.13 workflow/scripts/concatenate_contigs.py --input {input.ref} --output {output.ref_concat} --contig_order {output.contig_order} --sample_name {params.reference_name} --mode reference
        """

# concatenate the contigs in each sample consensus fasta file
rule concatenate_sample_contigs:
    input:
        consensus_fasta="results/{prefix}/consensus/{sample}/{sample}_consensus.fa",
        contig_order="results/{prefix}/consensus/{prefix}_ref_contig_order.tsv",
    output:
        sample_concat="results/{prefix}/consensus/{sample}/{sample}_consensus_concat.fa",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    conda:
        "envs/concatenate_contigs.yaml"
    shell:
        """
        python3.13 workflow/scripts/concatenate_contigs.py --input {input.consensus_fasta} --output {output.sample_concat} --contig_order {input.contig_order} --sample_name {wildcards.sample} --mode sample
        """

# combine the concatenated reference and sample fasta files into a single msa file
rule make_msa:
    input:
        ref_concat="results/{prefix}/consensus/{prefix}_ref_concat.fasta",
        sample_concat=expand(
            "results/{prefix}/consensus/{sample}/{sample}_consensus_concat.fa",
            sample=SAMPLES, prefix=PREFIX
        ),
    output:
        msa="results/{prefix}/alignment/{prefix}_consensus_msa.fa",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
    shell:
        """
        cat {input.ref_concat} {input.sample_concat} > {output.msa}
        """

