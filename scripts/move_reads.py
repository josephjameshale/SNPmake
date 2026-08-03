import os
import subprocess
import argparse
import pandas as pd

def move_reads(reads_dir, qc_file, target_clade, existing_output_dir, temp_dir, mode = 'move', prefixes = None, max_samples=200, suffix = '.fastq.gz'):
    m_df = pd.read_csv(qc_file, sep = '\t')
    # subset to only rows with 'PASS' in the QC_EVALUATION column
    m_pass = m_df[m_df['QC_EVALUATION'] == 'PASS']
    target_samples = m_pass[m_pass['auriclass_clade'] == target_clade]['Sample'].tolist()
    if mode == 'move':
        scmd = 'mv'
    elif mode == 'copy':
        scmd = 'cp'
    sample_count = 0
    for sample in target_samples:
        if prefixes is not None:
            if not any(sample.startswith(prefix) for prefix in prefixes):
                continue
        if existing_output_dir is not None:
            if os.path.isdir(os.path.join(existing_output_dir, sample)):
                # skip the sample if it's already present in the snpkit output
                continue
        for filename in os.listdir(reads_dir):
            if filename.startswith(sample):
                source_filename = None
                if filename.endswith('_R1.fastq.gz') or filename.endswith('_R1_trim_paired.fastq.gz'):
                    if filename.split('_R1')[0] == sample:
                        source_filename = filename
                        file_suffix = '_R1' + suffix
                elif filename.endswith('_R2.fastq.gz') or filename.endswith('_R2_trim_paired.fastq.gz'):
                    if filename.split('_R2')[0] == sample:
                        source_filename = filename
                        file_suffix = '_R2' + suffix
                if source_filename is None:
                    continue
                source_path = os.path.join(reads_dir, source_filename)
                target_path = os.path.join(temp_dir, f'{sample}{file_suffix}')
                # ensure the target is not being overwritten
                if os.path.exists(target_path):
                    print(f'Warning: target file {target_path} already exists! This file was not overwritten with {source_path}.')
                    continue
                subprocess.call([scmd, source_path, target_path])
                sample_count += 1
        if sample_count >= max_samples:
            break


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reads_dir','-rd',type=str,
        help='''Provide a path to the raw reads. Files to process must end in _R1.fastq.gz or _R2.fastq.gz''',
        default=None,required=True
        )
    parser.add_argument(
        '--temp_dir','-td',type=str,
        help='''Provide a path to move the reads to.''',
        default=None,required=True
        )
    parser.add_argument(
        '--qc_file','-q',type=str,
        help='''Provide a path to the QC file.''',
        default=None,required=True
        )
    parser.add_argument(
        '--clade','-c',type=str,
        help='''Provide the clade. This must exactly match the clade name in the master QC file.''',
        default=None,required=True
        )
    parser.add_argument(
        '--existing_output','-eo',type=str,
        help='''Provide a path to an output directory. If there are directories with the same names as samples in reads_dir, these samples will be skipped.''',
        default=None
        )
    parser.add_argument(
        '--mode','-m',type=str,choices=['move','copy'],
        help='''Specify if reads should be moved or copied to the temp directory.''',
        default='move'
        )
    parser.add_argument(
        '--prefixes','-p',type=str,nargs='+',
        help='''Provide a list of file name prefixes. If given, only files with these prefixes will be moved or copied.''',
        default=None
        )
    parser.add_argument(
        '--max_samples','-ms',type=int,
        help='''Provide the maximum number of samples to process.''',
        default=200
        )
    parser.add_argument(
        '--suffix','-s',type=str,
        help='''Provide the desired suffix for the files.''',
        default='.fastq.gz'
        )
    args = parser.parse_args()
    move_reads(args.reads_dir, args.qc_file, args.clade, args.existing_output, args.temp_dir, args.mode, args.prefixes, args.max_samples, args.suffix)

if __name__ == "__main__":
    main()


