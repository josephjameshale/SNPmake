# use the filtered vcf file to build a bed file of all positions that failed any filter
rule build_all_fail_bed:
    input:
        vcf_indelprox_lowcov="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz",
        vcf_indexprox_lowcov_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz.tbi",
    output:
        all_fail_bed = "results/{prefix}/bedtools/{sample}/{sample}_all_fail.bed",
    singularity:
        "docker://staphb/bcftools:1.23.1"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        """
        ### Note that only SNPs are included here!!
        bcftools query -i 'FILTER!="PASS" && TYPE="snp"' -f '%CHROM\t%POS0\t%END\n' {input.vcf_indelprox_lowcov} > {output.all_fail_bed}
        """


rule combine_bed_masks:
    input:
        all_fail_bed = "results/{prefix}/bedtools/{sample}/{sample}_all_fail.bed",
        lowcoverage_bed = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz",
        lowcoverage_bed_tbi = "results/{prefix}/bedtools/{sample}/{sample}_lowcoverage.bed.gz.tbi",
    output:
        final_bed = "results/{prefix}/bedtools/{sample}/{sample}_final.bed",
    singularity:
        "docker://staphb/bedtools:2.31.1"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        """
        # combine the all_fail_bed with the first three columns of the lowcoverage_bed file, and sort and merge the intervals
        cat {input.all_fail_bed} <(zcat {input.lowcoverage_bed} | cut -f1-3) | bedtools sort -i - | bedtools merge -i - > {output.final_bed}
        """



rule consensus_fasta:
    input:
        vcf_indelprox_lowcov="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz",
        vcf_indexprox_lowcov_tbi="results/{prefix}/filtered_vcf/{sample}/{sample}_bcftools_filtered_indelprox_lowcov.vcf.gz.tbi",
        final_bed = "results/{prefix}/bedtools/{sample}/{sample}_final.bed",
        ref_genome = config["reference_genome"],
    output:
        vcf_consensus = "results/{prefix}/consensus/{sample}/{sample}_pass_snp_only.vcf.gz",
        vcf_consensus_tbi = "results/{prefix}/consensus/{sample}/{sample}_pass_snp_only.vcf.gz.tbi",
        fasta="results/{prefix}/consensus/{sample}/{sample}_consensus.fa",
        duplicated_positions="results/{prefix}/consensus/{sample}/{sample}_duplicated_positions.txt"
    singularity:
        "docker://staphb/bcftools:1.23.1"
    benchmark: 
        "benchmarks/{prefix}/consensus_ref_allele_unmapped_variants/{sample}.benchmark.tsv"
    threads: 1
    resources:
            mem_mb=1000,
            runtime=10
    shell:
        """

        ### Note that the command below extracts only SNPs!! Indels will not be included in the final vcf or fasta file.
        bcftools view -f PASS -v snps -i 'GT="alt"' -Oz -o {output.vcf_consensus} {input.vcf_indelprox_lowcov}
        bcftools index -f -t {output.vcf_consensus}

        # ensure that there are no duplicated positions in this vcf file
        bcftools query -f '%CHROM\t%POS\n' {output.vcf_consensus} | sort -k1,1 -k2,2n | uniq -d > {output.duplicated_positions}
        if [ -s {output.duplicated_positions} ]; then
            echo "Error: Duplicated positions found in consensus vcf file. See {output.duplicated_positions} for details."
            exit 1
        fi

        # make a consensus fasta file from the vcf file, masking low-coverage positions with N
        VCF_SAMPLE=$(bcftools query -l {output.vcf_consensus} | head -n 1)
        bcftools consensus -f {input.ref_genome} -s $VCF_SAMPLE -p {wildcards.sample}_ -m {input.final_bed} --mask-with N {output.vcf_consensus} > {output.fasta}
        """
























# old version
# def consensus_vcf_list(wc):
#     lst = [
#         f"results/{wc.prefix}/unmapped_vcf/{wc.sample}/{wc.sample}_unmapped.vcf.gz",
#         f"results/{wc.prefix}/filtered_mask_vcf/{wc.sample}/{wc.sample}_filtered_mask.vcf.gz",
#         f"results/{wc.prefix}/indel_prox_mask_vcf/{wc.sample}/{wc.sample}_indel_prox_mask.vcf.gz",
#     ]
#     if PHAGE_ENABLED:
#         lst.append(f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz")
#     lst.append(f"results/{wc.prefix}/filtered_vcf/{wc.sample}/{wc.sample}_5bp_indel_removed.vcf.gz")
#     return " ".join(lst)

# # Make reference allele consensus fasta for each sample using the unmapped variants, filtered mask variants, and pass variants
# rule consensus_ref_allele_unmapped_variants:
#     input:
#         ref=REF_GENOME,
#         pass_vcf="results/{prefix}/filtered_vcf/{sample}/{sample}_5bp_indel_removed.vcf.gz",
#         unmapped_vcf="results/{prefix}/unmapped_vcf/{sample}/{sample}_unmapped.vcf.gz",
#         filtered_vcf="results/{prefix}/filtered_mask_vcf/{sample}/{sample}_filtered_mask.vcf.gz",                                 
#         indel_prox_vcf="results/{prefix}/indel_prox_mask_vcf/{sample}/{sample}_indel_prox_mask.vcf.gz",
#         phage_vcf_gz=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz" if PHAGE_ENABLED else [],
#         phage_tbi=lambda wc: f"results/{wc.prefix}/phage_mask_vcf/{wc.sample}/{wc.sample}_phage_mask.vcf.gz.tbi" if PHAGE_ENABLED else [],
#     output:
#         fasta="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}/{sample}_ref_allele_unmapped_variants.fa"
#     params:
#         outdir="results/{prefix}/consensus_ref_allele_unmapped_variants/{sample}",
#         vcf_list=consensus_vcf_list
#     singularity:
#         "docker://staphb/bcftools:1.23.1"
#     benchmark: 
#         "benchmarks/{prefix}/consensus_ref_allele_unmapped_variants/{sample}.benchmark.tsv"
#     # envmodules:
#     #     "Bioinformatics",
#     #     "htslib",
#     #     "bcftools"
#     threads: 1
#     resources:
#             mem_mb=1000,
#             runtime=10
#     shell:
#         r"""
#         set -euo pipefail
#         mkdir -p {params.outdir}

#         VCF_SAMPLE=$(bcftools query -l {input.pass_vcf} | head -n 1)

#         # Build a single sample overlay VCF by concatenation (not merge)
#         bcftools concat -a -D {params.vcf_list} -O z -o {params.outdir}/{wildcards.sample}_tmp_concat.vcf.gz
#         bcftools sort -O z -o {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz {params.outdir}/{wildcards.sample}_tmp_concat.vcf.gz
#         tabix -f -p vcf {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz

#         # Drop het genotypes by setting them to missing
#         bcftools +setGT {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz -O z \
#         -o {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz \
#         -- -t q -n . -i 'GT="het"'
#         tabix -f -p vcf {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz

#         # Now missing GT (formerly het) becomes N in the consensus
#         bcftools consensus -s "$VCF_SAMPLE" -M N -f {input.ref} {params.outdir}/{wildcards.sample}_consensus_sources_homonly.vcf.gz > {output.fasta}

#         # bcftools consensus -s "$VCF_SAMPLE" -f {input.ref} {params.outdir}/{wildcards.sample}_consensus_sources.vcf.gz > {output.fasta} 
#         sed -i 's/>.*/>{wildcards.sample}/g' {output.fasta}
#         """