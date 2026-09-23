
import argparse
import os
import shutil
from Bio import SeqIO

def find_files(input_dir, prefix, mode):
    if mode == 'msa':
        msa_name = None
        # search through the alignment subdirectory
        msa_dir = os.path.join(input_dir,'alignment')
        if os.path.isdir(msa_dir):
            for f in os.listdir(msa_dir):
                if f == f'{prefix}_consensus_msa.fa':
                    msa_name = os.path.join(msa_dir, f)
        # if the msa file is not found, search through the input directory
        if msa_name is None:
            for f in os.listdir(input_dir):
                if f == f'{prefix}_consensus_msa.fa':
                    msa_name = os.path.join(input_dir, f)
        if msa_name is None:
            print(f'Error: could not locate {prefix}_consensus_msa.fa in {msa_dir}')
            quit(1)
        return msa_name
    if mode == 'contig_order':
        contig_order_name = None
        # search through the consensus subdirectory
        con_dir = os.path.join(input_dir,'consensus')
        if os.path.isdir(con_dir):
            for f in os.listdir(con_dir):
                if f == f'{prefix}_ref_contig_order.tsv':
                    contig_order_name = os.path.join(con_dir, f)
        # if the contig order file is not found, search through the input directory
        if contig_order_name is None:
            for f in os.listdir(input_dir):
                if f == f'{prefix}_ref_contig_order.tsv':
                    contig_order_name = os.path.join(input_dir, f)
        if contig_order_name is None:
            print(f'Error: could not locate {prefix}_ref_contig_order.tsv in {con_dir}')
            quit(1)
        return contig_order_name

def read_contig_order(contig_order_file):
    contig_order = []
    with open(contig_order_file, 'r') as fh:
        for line in fh:
            contig, length = line.rstrip().split('\t')
            contig_order.append((contig, int(length)))
    if len(contig_order) < 1:
        print(f'Error: unable to read contig order file {contig_order_file}')
        quit(1)
    return contig_order

def check_contig_order(contig_order1, contig_order2):
    order1 = read_contig_order(contig_order1)
    order2 = read_contig_order(contig_order2)
    if order1 != order2:
        print('Error: the two input SNPmake runs do not appear to have the same reference genome or contig order.')
        print(f'Please ensure that {contig_order1} and {contig_order2} match.')
        quit(1)

# def check_contig_order(contig_order1, contig_order2):
#     with open(contig_order1, 'r') as fh1, open(contig_order2, 'r') as fh2:
#         for line1, line2 in zip(fh1, fh2):
#             if line1 != line2:
#                 print('Error: the two input SNPmake runs do not appear to have the same reference genome or contig order.')
#                 print(f'Please ensure that {contig_order1} and {contig_order2} match.')
#                 quit(1)

def combine_msa(msa1, msa2, output_msa):
    with open(output_msa, 'w') as fh_out:
        ref_record = False
        ref_record2 = False
        ids_seen = set()
        # assume that the reference sequence is the first record in each msa file
        for record in SeqIO.parse(msa1, 'fasta'):
            if not ref_record:
                ref_name = record.id
                ref_seq = record.seq
                ref_len = len(record.seq)
                print(f'Found reference sequence {ref_name} with length {ref_len} in {msa1}')
                ref_record = True
                ids_seen.add(record.id)
                SeqIO.write(record, fh_out, 'fasta')
                continue
            if len(record.seq) != ref_len:
                print(f'Error: incorrect sequence length for {record.id} in {msa1}')
                os.remove(output_msa)
                quit(1)
            if record.id in ids_seen:
                print(f'Skipping duplicate sequence {record.id} in {msa1}')
                continue
            ids_seen.add(record.id)
            SeqIO.write(record, fh_out, 'fasta')
        for record in SeqIO.parse(msa2, 'fasta'):
            if not ref_record2:
                ref_record2 = True
                if record.id != ref_name or record.seq != ref_seq:
                    print(f'Error: the reference sequence {record.id} in {msa2} does not match the reference sequence {ref_name}in {msa1}')
                    os.remove(output_msa)
                    quit(1)
                print(f'Found reference sequence {ref_name} with length {len(record.seq)} in {msa2}')
                continue
            if len(record.seq) != ref_len:
                print(f'Error: incorrect sequence length for {record.id} in {msa2}')
                os.remove(output_msa)
                quit(1)
            if record.id in ids_seen:
                print(f'Skipping duplicate sequence {record.id} in {msa2}')
                continue
            ids_seen.add(record.id)
            SeqIO.write(record, fh_out, 'fasta')

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input1','-i1',type=str,
        help='''Provide one SNPmake output directory, or any directory containing MSA and contig order files. Directory name should be unique.''',
        required=True
        )
    parser.add_argument(
        '--input2','-i2',type=str,
        help='''Provide a second SNPmake output directory, or any directory containing MSA and contig order files. Directory name should be unique.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide an output directory for the combined multiple sequence alignment.''',
        required=True
        )
    args = parser.parse_args()
    prefix1 = os.path.basename(os.path.normpath(args.input1))
    prefix2 = os.path.basename(os.path.normpath(args.input2))
    msa1 = find_files(args.input1, prefix1, 'msa')
    msa2 = find_files(args.input2, prefix2, 'msa')
    contig_order1 = find_files(args.input1, prefix1, 'contig_order')
    contig_order2 = find_files(args.input2, prefix2, 'contig_order')
    check_contig_order(contig_order1, contig_order2)
    # if all checks were passed, make the output directories
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    prefix_output_dir = os.path.join(args.output, f'{prefix1}_{prefix2}')
    if not os.path.exists(prefix_output_dir):
        os.makedirs(prefix_output_dir)
    # make the output MSA and contig files in the same format as the input files
    output_msa = os.path.join(prefix_output_dir, f'{prefix1}_{prefix2}_consensus_msa.fa')
    print(f'Combining {msa1} and {msa2} into {output_msa}')
    combine_msa(msa1, msa2, output_msa)
    # also copy one of the contig order files to the output directory
    output_contig_order = os.path.join(prefix_output_dir, f'{prefix1}_{prefix2}_ref_contig_order.tsv')
    print(f'Copying {contig_order1} to {output_contig_order}')
    shutil.copyfile(contig_order1, output_contig_order)

if __name__ == "__main__":
    main()


