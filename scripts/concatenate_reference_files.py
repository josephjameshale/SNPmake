import os
import subprocess
import argparse
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy
import pandas as pd


def concatenate_reference_files(fasta_file, gff_file, prefix, outdir, buffer_size=1000):
    # take a reference fasta file containing multiple contigs and a gff file
    # concatenate the contigs into a single entry in a new fasta file, adding a buffer of Ns between each contig
    # also make a dataframe that contains the start and end coordinates of each contig before and after concatenation
    # also make a new gff file with the coordinates of each feature adjusted to match the new concatenated fasta file
    gff_db = None
    if gff_file is not None:
        gff_db = gff.create_db(gff_file,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    # fasta_records = {record.id:record for record in SeqIO.parse(fasta_file,'fasta')}
    fasta_outfile = os.path.join(outdir, f'{prefix}.fasta')
    gff_outfile = os.path.join(outdir, f'{prefix}.gff')
    fasta_offsets = {}
    fasta_chunks = []
    contig_df = pd.DataFrame(columns=['contig','original_start','original_end','new_start','new_end'])
    # record all lengths for all contigs
    current_gp_start = 1
    current_gp_end = 1
    current_offset = 0
    for record in SeqIO.parse(fasta_file,'fasta'):
        fasta_offsets[record.id] = current_offset
        current_offset += len(record.seq) + buffer_size
        fasta_chunks.append(str(record.seq))
        current_gp_end = current_gp_start + (len(record.seq)-1)
        contig_df.loc[len(contig_df)] = {'contig':record.id, 'original_start':1, 'original_end':len(record.seq), 'new_start':current_gp_start, 'new_end':current_gp_end}
        current_gp_start = current_gp_end + 1 + buffer_size
    # add buffer and write to fasta
    fasta_concat_seq = Seq(('N' * buffer_size).join(fasta_chunks))
    concat_record = SeqRecord(fasta_concat_seq,id=prefix,name=prefix,description="")
    SeqIO.write(concat_record, fasta_outfile, "fasta")
    # also write the contig dataframe to a csv file
    contig_df.to_csv(os.path.join(outdir, f'{prefix}_contig_coordinates.csv'), index=False)
    # use fasta_offsets to update the gff coordinates, if provided (dictionary keys should be in the same order as the appended sequences)
    if gff_db is not None:
        with open(gff_outfile, 'w') as gff_out:
            _ = gff_out.write('##gff version 3\n')
            for oldfeature in gff_db.all_features():
                feature = copy.copy(oldfeature)
                contig = oldfeature.seqid
                feature.start = oldfeature.start + fasta_offsets[contig]
                feature.end = oldfeature.end + fasta_offsets[contig]
                feature.seqid = prefix
                _ = gff_out.write(str(feature) + "\n")

    
def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta','-f',type=str,
        help='''Provide a path to the FASTA file.''',
        required=True
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide a path to the GFF file.''',
        required=True
        )
    parser.add_argument(
        '--prefix','-p',type=str,
        help='''Provide a prefix for the output files.''',
        required=True
        )
    parser.add_argument(
        '--buffer','-b',type=int,
        help='''Provide the buffer size for concatenating reference files. The specified number of Ns will be added between each contig in the reference fasta file.''',
        default=1000
        )
    parser.add_argument(
        '--outdir','-o',type=str,
        help='''Provide a path to the output directory.''',
        required=True
        )
    args = parser.parse_args()
    conatatenate_reference_files(args.fasta, args.gff, args.prefix, args.outdir, args.buffer)

if __name__ == "__main__":
    main()


