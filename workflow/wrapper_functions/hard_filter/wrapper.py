__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell

# params
ref_genome = snakemake.params.get("ref_genome", "")

# version 2

# snp filters
qual_snp_filter = snakemake.params.get("qual_snp_filter", "")
mq_snp_filter = snakemake.params.get("mq_snp_filter", "")
dp_snp_filter = snakemake.params.get("dp_snp_filter", "")
ad_snp_filter_total = snakemake.params.get("ad_snp_filter_total", "")
ad_snp_filter_ratio = snakemake.params.get("ad_snp_filter_ratio", "")
#fq_snp_filter = snakemake.params.get("fq_snp_filter", "")
#af_snp_filter = snakemake.params.get("af_snp_filter", "")

# indels filters
qual_indel_filter = snakemake.params.get("qual_indel_filter", "")
qd_indel_filter = snakemake.params.get("qd_indel_filter", "")
dp_indel_filter = snakemake.params.get("dp_indel_filter", "")
ad_indel_filter_total = snakemake.params.get("ad_indel_filter_total", "")
ad_indel_filter_ratio = snakemake.params.get("ad_indel_filter_ratio", "")
#mq_indel_filter = snakemake.params.get("mq_indel_filter", "")
#af_indel_filter = snakemake.params.get("af_indel_filter", "")


# first, filter snps using bcftools filter
eval_string = f'TYPE == "snp" && (QUAL<{qual_snp_filter} || INFO/MQ<{mq_snp_filter} || FORMAT/DP[0]<{dp_snp_filter} || FORMAT/AD[0:1]<{ad_snp_filter_total} || FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])<{ad_snp_filter_ratio})'
shell("bcftools filter -e \"{eval_string}\" -s 'FAIL_LOW_QUALITY' -Oz -o {snakemake.output.vcf_filter_temp} {snakemake.input.norm_vcf_gz}")

# next, filter indels using this temp file as input
eval_string = f'TYPE == "indel" && (QUAL<{qual_indel_filter} || INFO/QD<{qd_indel_filter} || FORMAT/DP[0]<{dp_indel_filter} || FORMAT/AD[0:1]<{ad_indel_filter_total} || FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1])<{ad_indel_filter_ratio})'
shell("bcftools filter -e \"{eval_string}\" -s 'FAIL_LOW_QUALITY' -Oz -o {snakemake.output.vcf_filter} {snakemake.output.vcf_filter_temp}")
shell("bcftools index -t {snakemake.output.vcf_filter}")

# if needed, filter the vcf file to only include variants that passed all filters
shell("bcftools view -f 'PASS' -Oz -o {snakemake.output.vcf_filter_pass_only} {snakemake.output.vcf_filter}")
shell("bcftools index -t {snakemake.output.snp_vcf_filter_pass_only}")






# # next, filter indels using GATK VariantFiltration
# eval_string2 = f' --filter-name FAIL_LOW_QUALITY --filter-expression "QUAL<{qual_indel_filter}"'
# eval_string2 += f' --filter-name FAIL_LOW_QD --filter-expression "QD<{qd_indel_filter}"'
# eval_string2 += f' --filter-name FAIL_LOW_DP --filter-expression "DP<{dp_indel_filter}"'
# eval_string2 += f' --genotype-filter-name FAIL_LOW_ALT_SUPPORT --genotype-filter-expression "AD[1] < {ad_indel_filter_total} || AD[1] < 0.90 * (AD[0] + AD[1])"'












# # snps
# dp_snp_filter = snakemake.params.get("dp_snp_filter", "")
# fq_snp_filter = snakemake.params.get("fq_snp_filter", "")
# mq_snp_filter = snakemake.params.get("mq_snp_filter", "")
# qual_snp_filter = snakemake.params.get("qual_snp_filter", "")
# af_snp_filter = snakemake.params.get("af_snp_filter", "")

# # indels
# dp_indel_filter = snakemake.params.get("dp_indel_filter", "")
# mq_indel_filter = snakemake.params.get("mq_indel_filter", "")
# qual_indel_filter = snakemake.params.get("qual_indel_filter", "")
# af_indel_filter = snakemake.params.get("af_indel_filter", "")


# # snps
# gatk_snp_filter_parameter_expression= "%s && %s && %s && %s && %s" % (dp_snp_filter, fq_snp_filter, mq_snp_filter, qual_snp_filter, af_snp_filter)
# shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter_snp_vcf} --variant {snakemake.input.final_raw_snp_vcf} --filter-expression \"{gatk_snp_filter_parameter_expression}\" --filter-name PASS_filter &> {snakemake.log.gatk_snp}") 
# shell("grep '#\|PASS_filter' {snakemake.output.filter_snp_vcf} > {snakemake.output.filter_snp_final}")
# shell("""
#        bgzip -c {snakemake.output.filter_snp_final} > {snakemake.output.zipped_filtered_snp_vcf} &&
#        tabix -p vcf -f {snakemake.output.zipped_filtered_snp_vcf}
#     """)

# shell("""
#        bgzip -c {snakemake.output.filter_snp_vcf} > {snakemake.output.zipped_filter_snp_vcf} &&
#        tabix -p vcf -f {snakemake.output.zipped_filter_snp_vcf}
#     """)

# # indels       
# gatk_indel_filter_parameter_expression="%s && %s && %s && %s" % (dp_indel_filter, mq_indel_filter, qual_indel_filter, af_indel_filter)
# shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter_indel_vcf} --variant {snakemake.input.final_raw_indel_vcf} --filter-expression \"{gatk_indel_filter_parameter_expression}\" --filter-name PASS_filter&> {snakemake.log.gatk_indels}")
# shell("grep '#\|PASS_filter' {snakemake.output.filter_indel_vcf} > {snakemake.output.filter_indel_final}")

# shell("""
#        bgzip -c {snakemake.output.filter_indel_final} > {snakemake.output.zipped_filtered_indel_vcf} &&
#        tabix -p vcf -f {snakemake.output.zipped_filtered_indel_vcf}
#     """)

# shell("""
#        bgzip -c {snakemake.output.filter_indel_vcf} > {snakemake.output.zipped_filter_indel_vcf} &&
#        tabix -p vcf -f {snakemake.output.zipped_filter_indel_vcf}
#     """)