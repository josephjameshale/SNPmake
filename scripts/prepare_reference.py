import os
import argparse
from Bio import SeqIO
import mask_telomeres
import concatenate_reference_files
import subprocess

# dependencies: gffutils, biopython, bwa, samtools, repeatmasker (singularity container), repeatmasker library

def find_files(input_dir):
    fasta_fname, gff_fname = None, None
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

def process_reference(input_dir, fasta_fname, gff_fname, singularity_repeatmasker, repeatmasker_lib):
    # move the original fasta and gff files to a new directory called "original_reference"
    original_reference_dir = os.path.join(input_dir, 'original_reference_files')
    fname = os.path.basename(fasta_fname).split('.')[0]
    os.makedirs(original_reference_dir, exist_ok=True)
    original_fasta_fname = os.path.join(original_reference_dir, os.path.basename(fasta_fname))
    os.rename(fasta_fname, original_fasta_fname)
    if gff_fname is not None:
        original_gff_fname = os.path.join(original_reference_dir, os.path.basename(gff_fname))
        os.rename(gff_fname, original_gff_fname)
    else:
        original_gff_fname = None
    # mask telomeres with default settings of 5000 bases masked and 2000 base buffer around features
    telomere_masked_fasta = os.path.join(original_reference_dir, f'{fname}_telomere_masked.fasta')
    mask_telomeres.mask_telomeres(original_fasta_fname, original_gff_fname, telomere_masked_fasta)
    # run repeatmasker on the telomere masked fasta file
    cmd1 = ['singularity', 'exec', singularity_repeatmasker, 'RepeatMasker', '-dir', original_reference_dir, '-lib', repeatmasker_lib, '-pa', '1', telomere_masked_fasta]
    print(' '.join(cmd1))
    subprocess.run(cmd1)
    # concatenate the contigs in the masked fasta file together, using a default buffer size of 1000 N's
    concatenate_reference_files.concatenate_reference_files(telomere_masked_fasta, original_gff_fname, fname, input_dir)
    # the concatenated and masked fasta should now be at fasta_fname
    # index the new fasta file with bwa and samtools
    cmd2 = ['bwa', 'index', fasta_fname]
    print(' '.join(cmd2))
    subprocess.run(cmd2)
    cmd3 = ['samtools', 'faidx', fasta_fname]
    print(' '.join(cmd3))
    subprocess.run(cmd3)
    cmd4 = ['samtools', 'dict', fasta_fname, '-o', os.path.join(input_dir, f'{fname}.dict')]
    print(' '.join(cmd4))
    subprocess.run(cmd4)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to an input directory. This should contain exactly one .fasta file and one .gff file for the reference genome.''',
        required=True
        )
    parser.add_argument(
        '--singularity_repeatmasker','-s',type=str,
        help='''Provide a path to the singluarity container with RepeatMasker.''',
        default='/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2025_funQCD_database/tetools_2.0.0_07-15-26.sif'
        )
    parser.add_argument(
        '--repeatmasker_lib','-r',type=str,
        help='''Provide a path to the RepeatMasker library.''',
        default='/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2025_funQCD_database/lib/repeat_libraries/fungi_b8441/b8441_fungi_repeatlib.fa'
        )
    args = parser.parse_args()
    fasta_fname, gff_fname = find_files(args.input)
    # check number of sequences in fasta file
    num_seqs = sum(1 for record in SeqIO.parse(fasta_fname, 'fasta'))
    if num_seqs == 1:
        print(f'Reference genome appears to be already masked and concatenated (found 1 contig in {fasta_fname}).')
    else:
        process_reference(args.input, fasta_fname, gff_fname, args.singularity_repeatmasker, args.repeatmasker_lib)

if __name__ == "__main__":
    main()


