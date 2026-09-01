
import argparse
import os
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

def make_contig_order(seq_dict, contig_order_file):
    with open(contig_order_file,'w') as fh:
        for contig, record in seq_dict.items():
            _ = fh.write(f'{contig}\t{len(record.seq)}\n')

# def make_contig_order(seq_dict, contig_order_file):
#     # determine length of each contig, keeping the same order as the input fasta file
#     contig_lengths = [(x, len(seq_dict[x])) for x in seq_dict]
#     # if needed, contig can be sorted from longest to shortest
#     # contig_lengths = sorted(contig_lengths, key = lambda x:-x[1])
#     with open(contig_order_file,'w') as fh:
#         for contig, length in contig_lengths:
#             _ = fh.write(f'{contig}\t{length}\n')

def concatenate_contigs(input_fasta, output_fasta, contig_order_file, sample_name, mode):
    # read in the fasta file and store the sequences in a dictionary
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    # create the contig order file if in reference mode
    if mode == 'reference':
        # insertion order should be preserved in fasta_dict, so the order of the contigs should not change
        make_contig_order(fasta_dict, contig_order_file)
    # read in the contig order file
    if not os.path.isfile(contig_order_file):
        print('Failed to locate or generate contig order file')
        quit(1)
    contig_order = []
    with open(contig_order_file, 'r') as fh:
        for line in fh:
            contig, length = line.rstrip().split('\t')
            contig_order.append((contig, int(length)))
    if len(contig_order) < 1:
        print('Unable to read contig order file')
        quit(1)
    # concatenate the fasta records based on contig_order
    concat_record = Seq.Seq('')
    for contig_name, contig_length in contig_order:
        # each contig in the input fasta should match exactly one record in fasta_dict
        matching_fasta_record = [x for x in fasta_dict if x == contig_name]
        if len(matching_fasta_record) == 1 and len(fasta_dict[matching_fasta_record[0]].seq) == contig_length:
            concat_record += fasta_dict[matching_fasta_record[0]].seq
        else:
            if len(matching_fasta_record) == 0:
                matching_fasta_record = ['not found']
            print('Unable to match input fasta record name to reference contig name:')
            print(f'reference contig: {contig_name}')
            print(f'input fasta records: {",".join(matching_fasta_record)}')
            quit(1)
    concat_out = SeqRecord(concat_record, id=sample_name, description='')
    # make sure the output record is the expected length
    expected_length = sum([x[1] for x in contig_order])
    if len(concat_out.seq) != expected_length:
        print(f'Error: concatenated fasta record length {len(concat_out.seq)} does not match expected length {expected_length}')
        quit(1)
    with open(output_fasta, 'w') as fh:
        SeqIO.write(concat_out, fh, 'fasta')
    

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide an input fasta file.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide an output fasta file. This will have the contigs concatenated together.''',
        required=True
        )
    parser.add_argument(
        '--contig_order','-co',type=str,
        help='''Provide a two-column tsv file with the order of the contigs in the reference genome. If run in reference mode, a file will be created at the path provided.''',
        required=True
        )
    parser.add_argument(
        '--sample_name','-sn',type=str,
        help='''Provide a sample name for the output fasta file. This will be used as the fasta header in the output file.''',
        required=True
        )
    parser.add_argument(
        '--mode','-m',type=str,choices=['reference','sample'],
        help='''Specify whether to run in reference mode or sample mode. A contig order file is required for sample mode, but will be created in reference mode.''',
        required=True
        )
    args = parser.parse_args()
    concatenate_contigs(args.input, args.output, args.contig_order, args.sample_name, args.mode)


if __name__ == "__main__":
    main()


