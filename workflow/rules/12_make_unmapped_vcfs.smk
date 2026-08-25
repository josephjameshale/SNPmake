# rather than making a separate vcf for low-coverage regions, modify the existing vcf to flag variants in low-coverage regions
# do this with bcftools annotate and the low-coverage bed file
rule annotate_low_coverage_variants:
    input:
        vcf_indelprox="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf.gz",
        vcf_indexprox_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox.vcf.gz.tbi",
        lowcoverage_bed = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz",
        lowcoverage_bed_tbi = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz.tbi",
        lowcoverage_header = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage_header.hdr",        
    output:
        vcf_indelprox_lowcov="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz",
        vcf_indexprox_lowcov_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz.tbi",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    shell:
        """
        bcftools annotate -a {input.lowcoverage_bed} -h {input.lowcoverage_header} -c CHROM,FROM,TO,+FILTER -Oz -o {output.vcf_indelprox_lowcov} {input.vcf_indelprox}
        bcftools index -t {output.vcf_indelprox_lowcov}
        """


