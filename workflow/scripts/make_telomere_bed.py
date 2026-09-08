
import os
import argparse
import gffutils as gff
from Bio import SeqIO


def find_files(input_dir):
    fasta_fname, gff_fname = None, None
    # input_dir may actually be provided as a path to a fasta file
    if os.path.isfile(input_dir) and (input_dir.endswith('.fasta') or input_dir.endswith('.fa')):
        input_dir = os.path.dirname(input_dir)
    for fname in os.listdir(input_dir):
        if fname.endswith('.fasta') or fname.endswith('.fa'):
            if fasta_fname is not None:
                print(f'Error: found multiple fasta files in {input_dir}: {fasta_fname} and {fname}')
                quit(1)
            fasta_fname = os.path.join(input_dir, fname)
        elif fname.endswith('.gff') or fname.endswith('.gff3'):
            if gff_fname is not None:
                print(f'Error: found multiple gff files in {input_dir}: {gff_fname} and {fname}')
                quit(1)
            gff_fname = os.path.join(input_dir, fname)
    if fasta_fname is None:
        print(f'Error: could not find a fasta file in {input_dir}')
        quit(1)
    if gff_fname is None:
        print(f'Warning: could not find a gff file in {input_dir}. Genes near telomeres may be masked.')
    if fasta_fname is not None and gff_fname is not None:
        if os.path.basename(fasta_fname).split('.')[0] != os.path.basename(gff_fname).split('.')[0]:
            print(f'Error: fasta and gff files do not have the same base name: {fasta_fname} and {gff_fname}')
            quit(1)
    return fasta_fname, gff_fname


def mask_telomeres(input_fasta, input_gff, output_bed, bases_masked=5000, feature_buffer=2000):
    # take a gff file and a fasta file
    # create a bed file of regions to mask at the beginning and end of each scaffold
    # stop the interval early if it is within [feature_buffer] of any annotated feature in the gff file
    db = None
    if input_gff is not None:
        db = gff.create_db(input_gff, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    with open(output_bed, 'w') as fh_out:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            record_length = len(record.seq)
            mask_start = bases_masked
            mask_end = record_length - bases_masked
            if db is not None:
                for feature in db.region(seqid=record.id, featuretype=('gene','CDS','mRNA','exon')):
                    # note that gff files are 1-based
                    feature_start = feature.start - feature_buffer - 1
                    feature_end = feature.end + feature_buffer
                    if feature_start < mask_start:
                        print(f'Found {feature.id} at {feature_start}')
                        mask_start = feature_start
                    if feature_end > mask_end:
                        print(f'Found {feature.id} at {feature_end}')
                        mask_end = feature_end
            if mask_end <= mask_start:
                print(f'Warning: scaffold {record.id} will be fully masked.')
                _ = fh_out.write(f'{record.id}\t0\t{record_length}\n')
            else:
                if mask_start > 0:
                    print(f'{record.id}: masking from 0 to {mask_start}')
                    _ = fh_out.write(f'{record.id}\t0\t{mask_start}\n')
                if mask_end < record_length:
                    print(f'{record.id}: masking from {mask_end} to {record_length}')
                    _ = fh_out.write(f'{record.id}\t{mask_end}\t{record_length}\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to an input directory. This should contain exactly one .fasta file and one .gff file for the reference genome.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to the output BED file.''',
        required=True
        )
    parser.add_argument(
        '--scaffold_mask_size','-s',type=int,
        help='''Specify how many bases to mask at each scaffold.''',
        default=5000
        )
    parser.add_argument(
        '--feature_buffer_size','-f',type=int,
        help='''Specify the buffer size around each feature to include in the BED file.''',
        default=2000
        )
    args = parser.parse_args()
    fasta_fname, gff_fname = find_files(args.input)
    mask_telomeres(fasta_fname, gff_fname, args.output, bases_masked=args.scaffold_mask_size, feature_buffer=args.feature_buffer_size)


if __name__ == "__main__":
    main()


