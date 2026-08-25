
# Note that duplicates are marked but not entirely removed in this step
rule picard_remove_duplicates:
    input:
        sorted_bam_out = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam"
    output:
        picard_bam = temp("results/{prefix}/post_align/{sample}/picard_remove_duplicates/{sample}_aln_marked.bam"),
    params:
        picard_metrics_out = "results/{prefix}/post_align/{sample}/picard_remove_duplicates/{sample}_markduplicates_metrics"
    log:
        "logs/{prefix}/post_align/{sample}/{sample}_picard.log"
    benchmark: 
        "benchmarks/{prefix}/picard_remove_duplicates/{sample}_picard.benchmark.tsv"
    singularity:
        "docker://broadinstitute/picard:3.5.0"
    threads: 2
    resources:
        mem_mb=5000,
        runtime=15
    shell:
        """
        picard MarkDuplicates -REMOVE_DUPLICATES false -INPUT {input.sorted_bam_out} -OUTPUT {output.picard_bam} -METRICS_FILE {params.picard_metrics_out} -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT &> {log}
        """


rule samtools_sort_index:
    input:
        picard_bam = "results/{prefix}/post_align/{sample}/picard_remove_duplicates/{sample}_aln_marked.bam"
    output:
        sorted_picard_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
        sorted_picard_bam_bai = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam.bai",
    benchmark: 
        "benchmarks/{prefix}/picard_remove_duplicates/{sample}_samtools.benchmark.tsv"
    threads: 2
    resources:
        mem_mb=5000,
        runtime=15
    singularity:
        "docker://staphb/samtools:1.24"
    params:
        temp_dir = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_aln_sort_temp",
        sub_threads = lambda wildcards, threads: max(0, threads - 1),
    shell:
        """
        samtools sort {input.picard_bam} -m 500M -@ {params.sub_threads} -o {output.sorted_picard_bam} -T {params.temp_dir}
        samtools index {output.sorted_picard_bam}
        """

