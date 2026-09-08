__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import gzip
from snakemake.shell import shell
import pandas as pd

def add_filter(existing_filter, new_filter):
    """Append a FILTER label while preserving valid VCF semantics."""
    if existing_filter in (".", "PASS", ""):
        return new_filter
    labels = existing_filter.split(";")
    if new_filter not in labels:
        labels.append(new_filter)
    return ";".join(labels)

def vcf_is_snp_record(ref_record, alt_record):
    """Determine if a VCF record is a SNP based on the REF and ALT fields, assuming biallelic records only."""
    len_check = len(ref_record) == 1 and len(alt_record) == 1
    sym_check = alt_record not in ('.', '*')
    return len_check and sym_check

def vcf_gz_to_pandas_df(vcf_gz_file):
    # read in a vcf.gz file and return a pandas dataframe
    with gzip.open(vcf_gz_file, 'rt') as f:
        # skip header lines that start with ##
        for line in f:
            if line.startswith('#CHROM'):
                colnames = line.strip().split('\t')
                # remove the leading # from the first column name
                colnames[0] = colnames[0].lstrip('#')
                break
    # read the rest of the vcf.gz with pandas
    df = pd.read_csv(vcf_gz_file, sep='\t', comment='#', header=None, names=colnames)
    return(df)

def find_indel_prox_snps(vcf_df, prox_window = 5):
    # find SNPs that are within prox_window bp of an indel
    vcf_df_indels = vcf_df[vcf_df['INFO'].str.contains('INDEL')]
    vcf_df_snps = vcf_df[~vcf_df['INFO'].str.contains('INDEL')]
    # ensure everything in vcf_df_snps is a SNP record
    vcf_df_snps = vcf_df_snps[vcf_df_snps.apply(lambda row: vcf_is_snp_record(row['REF'], row['ALT']), axis=1)]
    # take both the start and end positions of each indel
    indel_positions = {}
    for index, row in vcf_df_indels.iterrows():
        indel_start = row['POS']
        indel_chrom = row['CHROM']
        indel_end = indel_start + len(row['REF']) - 1
        indel_positions[indel_chrom] = indel_positions.get(indel_chrom, []) + [indel_start, indel_end]
    # find SNPs that are within prox_window bp of an indel
    # if there are multiple records at the same position, they will all be masked
    indel_prox_snp_positions = set()
    for index, row in vcf_df_snps.iterrows():
        snp_chrom = row['CHROM']
        snp_pos = row['POS']
        for indel_pos in indel_positions.get(snp_chrom, []):
            if abs(snp_pos - indel_pos) <= prox_window:
                indel_prox_snp_positions.add((snp_chrom, snp_pos))
                break
    return(indel_prox_snp_positions)

def mask_indel_prox_snps(vcf_gz_file, output_vcf_file, excluded_positions_file, prox_window = 5):
    # convert vcf.gz to pandas df
    vcf_df = vcf_gz_to_pandas_df(vcf_gz_file)
    # find SNPs that are within prox_window bp of an indel
    indel_prox_snp_positions = find_indel_prox_snps(vcf_df, prox_window = prox_window)
    # write the positions to a file
    with open(excluded_positions_file, "w") as f:
        for chrom, pos in sorted(indel_prox_snp_positions):
            f.write(f"{chrom}\t{pos}\n")
    # mask the SNPs in the vcf_df
    # retain all comment lines and structure by reading line by line
    with gzip.open(vcf_gz_file, 'rt') as f_in, open(output_vcf_file, 'w') as f_out:
        filter_declared = False
        for line in f_in:
            if line.startswith('#'):
                if line.startswith('##FILTER=') and not filter_declared:
                    f_out.write('##FILTER=<ID=FAIL_INDEL_PROX,Description="SNP is within 5 bp of an indel">\n')
                    filter_declared = True
                f_out.write(line)
            else:
                fields = line.rstrip('\n').split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                info = fields[7]
                # skip any rows that are indels
                if 'INDEL' not in info and vcf_is_snp_record(fields[3], fields[4]):
                    if (chrom, pos) in indel_prox_snp_positions:
                        fields[6] = add_filter(fields[6], 'FAIL_INDEL_PROX')
                        f_out.write('\t'.join(fields) + '\n')
                    else:
                        f_out.write(line)
                else:
                    f_out.write(line)


input_vcf_gz_file = snakemake.input.vcf_filter

output_vcf_temp = snakemake.output.vcf_indelprox_temp
#output_vcf = snakemake.output.vcf_indelprox
excluded_positions_file = snakemake.output.masked_positions

indel_prox_window = snakemake.params.indel_prox_window

mask_indel_prox_snps(input_vcf_gz_file, output_vcf_temp, excluded_positions_file, prox_window = indel_prox_window)
