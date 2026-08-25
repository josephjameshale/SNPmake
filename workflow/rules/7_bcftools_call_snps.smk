# changes made to this rule:
# specify ploidy of 1
# remove log
# replace the consensus caller (-c) with the multiallelic caller (-m) - note that this will require splitting multiallelic sites into biallelic sites in a later step
# remove the unnecessary bgzip intermediate step
# NOTE: despite the name, this output vcf contains both SNPs and indels!
rule bcftools_call_variants:
    input:
        index_sorted_dups_rmvd_bam = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam",
        dups_rmvd_sorted_bam_out_bai = "results/{prefix}/post_align/{sample}/sorted_bam_dups_removed/{sample}_final.bam.bai",
        ref_genome = config["reference_genome"],
    output:
        raw_vcf_gz= "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf.gz",
        raw_vcf_tbi= "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf.gz.tbi",
    params:
        annotation_str = "FORMAT/AD,FORMAT/DP,INFO/AD"
    singularity:
        "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/bcftools_call_snps/{sample}.benchmark.tsv"
    shell:
        """
        bcftools mpileup -Ou -a {params.annotation_str} -f {input.ref_genome} {input.index_sorted_dups_rmvd_bam} | bcftools call -Oz -v -m --ploidy 1 -o {output.raw_vcf_gz}
        bcftools index -t {output.raw_vcf_gz}
        """

# take the raw vcf file and normalize it to split multiallelic sites into biallelic sites
# also remove any any duplicate records - note that SNPs with alternate allele calls at the same position will still be retained here
# then subset to only SNPs, if desired
# this can be done by adding these steps:
    # bcftools view -v snps -Oz -o {output.snp_vcf_gz} {output.norm_vcf_gz}
    # bcftools index -t {output.snp_vcf_gz}
rule bcftools_normalize:
    input:
        raw_vcf_gz = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw.vcf.gz",
    output:
        norm_vcf_gz = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_norm.vcf.gz",
        norm_vcf_tbi = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_norm.vcf.gz.tbi",
        #snp_vcf_gz = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_snp.vcf.gz",
        #snp_vcf_tbi = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_snp.vcf.gz.tbi",
    params:
        ref_genome = config["reference_genome"]
    singularity:
        "docker://staphb/bcftools:1.23.1"
    shell:
        """
        bcftools norm --multiallelics -any --fasta-ref {params.ref_genome} --rm-dup exact -Oz -o {output.norm_vcf_gz} {input.raw_vcf_gz}
        bcftools index -t {output.norm_vcf_gz}
        """
    


