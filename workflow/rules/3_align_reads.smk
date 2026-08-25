
# note that -M was removed from the bwa-mem command
rule align_reads:
    input:
        r1=align_r1,
        r2=align_r2
    output:
        aligned_sam_out=temp("results/{prefix}/align_reads/{sample}/{sample}_aln.sam")
    params:
        ref_genome=config["reference_genome"],
        readgroup= lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA"
    log:
        "logs/{prefix}/align_reads/{sample}/{sample}.log"
    singularity:
        "docker://staphb/bwa:0.7.19"
    benchmark: 
        "benchmarks/{prefix}/align_reads/{sample}.benchmark.tsv"
    threads: 8
    resources:
        mem_mb=3000,
        runtime=15
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        bwa mem -R "{params.readgroup}" -t {threads} {params.ref_genome} {input.r1} {input.r2} > {output.aligned_sam_out} 2> {log}
        """
