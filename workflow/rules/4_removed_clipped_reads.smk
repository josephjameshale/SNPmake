
rule samtools_sam_to_bam:
    input:
        aligned_sam_out = "results/{prefix}/align_reads/{sample}/{sample}_aln.sam",
        ref_genome= config["reference_genome"],
    output:
        bam_out = temp("results/{prefix}/post_align/{sample}/aligned_bam/{sample}_aln.bam"),
        sorted_bam_out = temp("results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam")
    log:
        "logs/{prefix}/post_align/{sample}/{sample}_samtools.log",
    singularity:
        "docker://staphb/samtools:1.24"
    benchmark: 
        "benchmarks/{prefix}/remove_clipped_reads/{sample}_samtools.benchmark.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=90
    params:
        outdir_temp = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort_temp",
        sub_threads = lambda wildcards, threads: max(0, threads - 1),
    shell:
        """
        samtools sort {input.aligned_sam_out} -m 500M -@ {params.sub_threads} -o {output.sorted_bam_out} -T {params.outdir_temp} &> {log}
        """





# below is a version that preserves the samclip step
# rule samclip:
#     input:
#         aligned_sam_out = "results/{prefix}/align_reads/{sample}/{sample}_aln.sam"
#     output:
#         clipped_sam_out = temp("results/{prefix}/post_align/{sample}/samclip/{sample}_clipped.sam"),
#         # bam_out = temp("results/{prefix}/post_align/{sample}/aligned_bam/{sample}_aln.bam"),
#         # sorted_bam_out = temp("results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam")
#     params:
#         outdir_temp = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort_temp",
#         prefix = "{sample}",
#         ref_genome= config["reference_genome"]
#     singularity:
#         "docker://staphb/samclip:0.4.0"
#     benchmark: 
#         "benchmarks/{prefix}/remove_clipped_reads/{sample}_samclip.benchmark.tsv"
#     threads: 1
#     resources:
#         mem_mb=2000,
#         runtime=90
#     shell:
#         """
#         samclip --ref {params.ref_genome} --max 10 < {input.aligned_sam_out} > {output.clipped_sam_out}
#         """



# rule samtools_sam_to_bam:
#     input:
#         clipped_sam_out = "results/{prefix}/post_align/{sample}/samclip/{sample}_clipped.sam",
#     output:
#         bam_out = temp("results/{prefix}/post_align/{sample}/aligned_bam/{sample}_aln.bam"),
#         sorted_bam_out = temp("results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort.bam")
#     params:
#         outdir_temp = "results/{prefix}/post_align/{sample}/sorted_bam/{sample}_aln_sort_temp",
#         prefix = "{sample}",
#         ref_genome= config["reference_genome"]
#     log:
#         "logs/{prefix}/post_align/{sample}/{sample}_samtools.log",
#     singularity:
#         "docker://staphb/samtools:1.24"
#     benchmark: 
#         "benchmarks/{prefix}/remove_clipped_reads/{sample}_samtools.benchmark.tsv"
#     threads: 1
#     resources:
#         mem_mb=2000,
#         runtime=90
#     shell:
#         """
#         samtools view -Sb {input.clipped_sam_out} > {output.bam_out} && samtools sort {output.bam_out} -m 500M -@ 0 -o {output.sorted_bam_out} -T {params.outdir_temp} &> {log}
#         """
