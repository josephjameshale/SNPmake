import os
import argparse

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


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to the directory of raw reads.''',
        required=True
        )
    parser.add_argument(
        '--trim_string','-t',type=str,
        help='''(Optional) Alternate string to trim from the end of each filename. The default is _R.''',
        default='_R'
        )
    args = parser.parse_args()
    modify_fastq_names(args.input, args.trim_string)

if __name__ == "__main__":
    main()


