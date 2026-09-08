# NOTES
# as written, both SNPs and indels will be present in the 'snp_vcf' output files
# (both SNPs and indels are present in the 'final_raw_snp_vcf' input file)
# my goal for editing these scripts is mainly to write a SINGLE vcf file that contains PASS (if the snp passes the filtering methods) and FAIL otherwise
# the upstream steps will be edited to split the vcf files, such that final_raw_snp_vcf contains only SNPs
# I may also need to change some of the filtering parameters, since I used the multiallelic caller rather than the consensus caller

rule variant_hard_filter:
    input:
        norm_vcf_gz = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_norm.vcf.gz",
        norm_vcf_tbi = "results/{prefix}/bcftools_varcall/{sample}/{sample}_aln_mpileup_raw_norm.vcf.gz.tbi",
    output:
        vcf_filter_temp = temp("results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_temp.vcf.gz"),
        vcf_filter = "results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered.vcf.gz",
        vcf_filter_tbi = "results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered.vcf.gz.tbi",
        vcf_filter_pass_only = "results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_pass_only.vcf.gz",
        vcf_filter_pass_only_tbi = "results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_pass_only.vcf.gz.tbi",
    params:
        ref_genome = config["reference_genome"],
        qual_snp_filter = config["qual_snp_filter"],
        mq_snp_filter = config["mq_snp_filter"],
        dp_snp_filter = config["dp_snp_filter"],
        ad_snp_filter_alt_depth = config["ad_snp_filter_alt_depth"],
        ad_snp_filter_ratio = config["ad_snp_filter_ratio"],
        qual_indel_filter = config["qual_indel_filter"],
        mq_indel_filter = config["mq_indel_filter"],
        dp_indel_filter = config["dp_indel_filter"],
        ad_indel_filter_alt_depth = config["ad_indel_filter_alt_depth"],
        ad_indel_filter_ratio = config["ad_indel_filter_ratio"],
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60
    shell:
        r"""
        EVAL_STRING_SNP='TYPE == "snp" && (QUAL<{params.qual_snp_filter} || INFO/MQ<{params.mq_snp_filter} || FORMAT/DP[0]<{params.dp_snp_filter} || FORMAT/AD[0:1]<{params.ad_snp_filter_alt_depth} || FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])<{params.ad_snp_filter_ratio})'
        EVAL_STRING_INDEL='TYPE == "indel" && (QUAL<{params.qual_indel_filter} || INFO/MQ<{params.mq_indel_filter} || FORMAT/DP[0]<{params.dp_indel_filter} || FORMAT/AD[0:1]<{params.ad_indel_filter_alt_depth} || FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])<{params.ad_indel_filter_ratio})'

        # filter SNPs first
        bcftools filter -e "$EVAL_STRING_SNP" -s 'FAIL_SNP_QC' -Oz -o {output.vcf_filter_temp} {input.norm_vcf_gz}

        # then filter indels using this temp file
        bcftools filter -e "$EVAL_STRING_INDEL" -s 'FAIL_INDEL_QC' -Oz -o {output.vcf_filter} {output.vcf_filter_temp}
        bcftools index -t {output.vcf_filter}

        # make a vcf file that contains only PASS records
        bcftools view -f 'PASS' -Oz -o {output.vcf_filter_pass_only} {output.vcf_filter}
        bcftools index -t {output.vcf_filter_pass_only}
        """



