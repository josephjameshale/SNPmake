
import os
import subprocess
import argparse
import pandas as pd

def copy_gm_key():
    if not os.path.exists('./.gm_key'):
        homedirkey = os.path.expanduser('~') + '/.gm_key'
        if not os.path.exists(homedirkey):
            print('GeneMark license not found in home directory. Please download .gm_key from https://genemark.bme.gatech.edu/license_download.cgi')
        else:
            subprocess.call(['cp',homedirkey,'./.gm_key'])
            print('Copied GeneMark license to working directory')
    else:
        print('GeneMark license found in funQCD directory')

def make_samples_csv(path):
    flist = os.listdir(path)
    sample_id_dict = {}
    if os.path.exists('config/samples.csv'):
        print('Overwriting config/samples.csv with new version')
    for f in flist:
        if f.endswith('_R1.fastq.gz') or f.endswith('_R2.fastq.gz'):
            sample_id = '_R'.join(f.split('_R')[:-1])
            sample_id_dict[sample_id] = 1
            # this should always return the full text to the left of '_R1.fastq.gz' or '_R2.fastq.gz', even if it contains '_R' somewhere other than the end
        elif f.endswith('.fasta'):
            sample_id = f.split('.fasta')[0]
            sample_id_dict[sample_id] = 'assembly'
    with open('config/samples.csv','w') as fhout:
        _ = fhout.write('sample_id\n')
        for sid in sorted(sample_id_dict.keys()):
            _ = fhout.write(f'{sid}\n')

def add_path_to_config(path,prefix,config_file = 'config/config.yaml'):
    lines = []
    with open(config_file,'r') as fh:
        for line in fh:
            if line.startswith('short_reads:'):
                line = f'short_reads: {path}\n'
            if line.startswith('prefix:'):
                line = f'prefix: {prefix}\n'
            lines.append(line)
    with open(config_file,'w') as fh_out:
        for line in lines:
            _ = fh_out.write(line)

def modify_fastq_names(path, trim_string):
    # ensure all files start with the string preceding trim_string, and end in either _R1.fastq.gz or _R2.fastq.gz as appropriate
    flist = os.listdir(path)
    for fname in flist:
        if not fname.endswith('.fastq.gz') or fname.startswith('.'):
            continue
        fname_prefix = fname.split(trim_string)[0]
        # extract the character following the last instance of '_R' in the filename, and ensure it is either 1 or 2
        if fname.count('_R') == 0:
            print(f'Error: cannot determine read direction for {fname}')
            quit(1)
        if fname.count('_R') > 1:
            print(f'Warning: ambiguous read direction for {fname}. The last occurrence will be used to deteremine read direction.')
        fname_suffix = fname.split('_R')[-1]
        if fname_suffix[0] not in ['1','2']:
            print(f'Error: cannot determine read direction for {fname}')
            quit(1)
        new_fname = f'{fname_prefix}_R{fname_suffix[0]}.fastq.gz'
        if os.path.exists(os.path.join(path,new_fname)):
            print(f'Error: {fname} shares the same name as another file after modification')
            quit(1)
        os.rename(os.path.join(path,fname),os.path.join(path,new_fname))
        print(f'Renamed {fname} to {new_fname}')


def variant_calling_setup(funQCD_dir, temp_reads_dir):
    # search through the funQCD output for all auriclass output directories
    # use these to determine which isolates to put in each clade directory in temp_reads_dir
    # assume that trimmed reads are in a subdirectory of funQCD_dir called 'trimmomatic'
    for sample_name in os.listdir(os.path.join(funQCD_dir,'auriclass')):
        auriclass_report_file = os.path.join(funQCD_dir,'auriclass',sample_name,f'{sample_name}_report.tsv')
        auriclass_df = pd.read_csv(auriclass_report_file,sep='\t')
        sample_clade = auriclass_df['Clade'][0]
        sample_clade = sample_clade.replace(' ','_')
        # copy this sample's reads to a temp_reads_dir subdirectory for this clade
        clade_dir = os.path.join(temp_reads_dir,sample_clade)
        if not os.path.exists(clade_dir):
            os.makedirs(clade_dir)
        reads_r1 = os.path.join(funQCD_dir,'trimmomatic',sample_name,f'{sample_name}_R1_trim_paired.fastq.gz')
        reads_r2 = os.path.join(funQCD_dir,'trimmomatic',sample_name,f'{sample_name}_R2_trim_paired.fastq.gz')
        if not os.path.exists(reads_r1) or not os.path.exists(reads_r2):
            print(f'Warning: could not locate trimmed reads for {sample_name} in {funQCD_dir}')
            continue
        subprocess.call(['cp',reads_r1,clade_dir])
        subprocess.call(['cp',reads_r2,clade_dir])
        # add this sample to a samples.csv file in temp_reads_dir
        clade_samples_csv = os.path.join(temp_reads_dir,f'samples_{sample_clade}.csv')
        if not os.path.isfile(clade_samples_csv):
            with open(clade_samples_csv,'w') as fh:
                _ = fh.write('sample_id\n')
        with open(clade_samples_csv,'a') as fh:
            _ = fh.write(f'{sample_name}\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--funQCD','-f',type=str,
        help='''Provide a path to funQCD output directory with auriclass results.''',
        default=None,required=True
        )
    parser.add_argument(
        '--temp_reads_dir','-t',type=str,
        help='''Provide a path to a temporary directory of reads to be processed with variant calling.''',
        default=None
        )
    args = parser.parse_args()
    if not os.path.isdir(args.funQCD):
        print(f'Could not locate directory at {args.funQCD}')
        quit(1)
    if not os.path.isdir(args.temp_reads_dir):
        os.makedirs(args.temp_reads_dir)
    variant_calling_setup(args.funQCD, args.temp_reads_dir)

if __name__ == "__main__":
    main()


