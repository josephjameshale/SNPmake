
import argparse
import pandas as pd

def read_snp_data(snp_data_file):
    # read in a snp_data file and return a pandas dataframe
    df = pd.read_csv(snp_data_file, sep='\t', dtype={'CHROM': str, 'POS': int})
    # check for any duplicated combinations of CHROM and POS
    if df.duplicated(subset=['CHROM', 'POS']).any():
        print(f'Error: Duplicate SNP positions found in {snp_data_file}.')
        quit(1)
    return(df)

def merge_snp_data(snp_data_files, output_file, fail_bed, fail_threshold = 80.0):
    # read in all snp_data files and concatenate them into a single dataframe, preserving all rows
    merged_df = read_snp_data(snp_data_files[0])
    for snp_data_file in snp_data_files[1:]:
        df = read_snp_data(snp_data_file)
        merged_df = pd.merge(merged_df, df, how='outer', on=['CHROM', 'POS'])
    # fill all missing values with 'ABSENT'
    merged_df.fillna('ABSENT', inplace=True)
    # use all columns that don't start with 'CHROM' or 'POS' as sample columns, and determine the percentage of 'PASS' and 'ABSENT' for each row
    # also determine the percentage of non-PASS and non-ABSENT values for each row
    sample_columns = [col for col in merged_df.columns if not col.startswith('CHROM') and not col.startswith('POS')]
    merged_df['PASS_PERCENT'] = merged_df[sample_columns].apply(lambda row: (row == 'PASS').sum() / len(row) * 100, axis=1)
    merged_df['ABSENT_PERCENT'] = merged_df[sample_columns].apply(lambda row: (row == 'ABSENT').sum() / len(row) * 100, axis=1)
    merged_df['FAIL_PERCENT'] = merged_df[sample_columns].apply(lambda row: ((row != 'PASS') & (row != 'ABSENT')).sum() / len(row) * 100, axis=1)
    merged_df['FAIL_PERCENT_PRESENT'] = merged_df[sample_columns].apply(lambda row: ((row != 'PASS') & (row != 'ABSENT')).sum() / (row != 'ABSENT').sum() * 100, axis=1)
    # write the merged dataframe to a file
    merged_df.to_csv(output_file, sep='\t', index=False)
    # determine which positions have FAIL_PERCENT_PRESENT >= to the fail threshold
    fail_positions = merged_df[merged_df['FAIL_PERCENT_PRESENT'] >= fail_threshold]
    fail_bed_df = pd.DataFrame()
    fail_bed_df['CHROM'] = fail_positions['CHROM']
    fail_bed_df['START'] = fail_positions['POS'] - 1
    fail_bed_df['END'] = fail_positions['POS']
    # write the fail bed file
    fail_bed_df.to_csv(fail_bed, sep='\t', index=False, header=False)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str, nargs='+',
        help='''Provide a list of paths to the input snp_data files.''',
        default=None,required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to the output file.''',
        default=None
        )
    parser.add_argument(
        '--fail_threshold','-f',type=float,
        help='''Provide a percentage threshold for failing positions. If a variant fails QC 
        in at least this percentage of isolates, it will be added to the fail bed file.''',
        default=None
        )
    parser.add_argument(
        '--fail_bed','-b',type=str,
        help='''Provide a path to the output fail bed file. This will be used to mask failing positions in the consensus fasta.''',
        default=None
        )
    args = parser.parse_args()
    merge_snp_data(args.input, args.output, args.fail_bed, args.fail_threshold)

if __name__ == "__main__":
    main()


