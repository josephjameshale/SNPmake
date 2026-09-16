import os
import argparse
import subprocess

def prepare_snpeff_database(input_dir):
    # take an input directory containing a reference genome fasta file and a gff3 file
    # create the appropriate snpeff database directory structure and create symbolic links for sequences.fa and genes.gff
    # assume the directory name is the name of the reference genome
    ref_genome_name = os.path.basename(os.path.normpath(input_dir))
    if not os.path.isfile(os.path.join(input_dir, f"{ref_genome_name}.fasta")):
        raise ValueError(f"Reference genome fasta file {ref_genome_name}.fasta not found in {input_dir}")
    else:
        fasta_path = os.path.join(os.path.abspath(input_dir), f"{ref_genome_name}.fasta")
    if not os.path.isfile(os.path.join(input_dir, f"{ref_genome_name}.gff3")):
        raise ValueError(f"Reference genome gff3 file {ref_genome_name}.gff3 not found in {input_dir}")
    else:
        gff3_path = os.path.join(os.path.abspath(input_dir), f"{ref_genome_name}.gff3")
    # Create the snpeff database directory structure
    os.makedirs(os.path.join(input_dir, "snpeff"), exist_ok=True)
    os.makedirs(os.path.join(input_dir, "snpeff", "data"), exist_ok=True)
    os.makedirs(os.path.join(input_dir, "snpeff", "data", ref_genome_name), exist_ok=True)
    # Create symbolic links for sequences.fa and genes.gff
    os.symlink(fasta_path, os.path.join(input_dir, "snpeff", "data", ref_genome_name, "sequences.fa"))
    os.symlink(gff3_path, os.path.join(input_dir, "snpeff", "data", ref_genome_name, "genes.gff"))
    # if .proteins.fa or .cds-transcripts.fa files exist, create symbolic links for them as well
    proteins_path = os.path.join(os.path.abspath(input_dir), f"{ref_genome_name}.proteins.fa")
    cds_path = os.path.join(os.path.abspath(input_dir), f"{ref_genome_name}.cds-transcripts.fa")
    if os.path.isfile(proteins_path):
        os.symlink(proteins_path, os.path.join(input_dir, "snpeff", "data", ref_genome_name, "protein.fa"))
    if os.path.isfile(cds_path):
        os.symlink(cds_path, os.path.join(input_dir, "snpeff", "data", ref_genome_name, "cds.fa"))
    # Create the snpeff.config file
    with open(os.path.join(input_dir, "snpeff", "snpEff.config"), "w") as f:
        f.write(f"data.dir = {os.path.join(os.path.abspath(input_dir), 'snpeff', 'data')}\n")
        f.write(f"{ref_genome_name}.genome : {ref_genome_name}\n")


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to the input directory.''',
        required=True
        )
    parser.add_argument(
        '--mode','-m',choices=['prepare', 'build'],type=str,
        help='''Specify if you want to prepare the directories or attempt to build the database.''',
        required=True
        )
    parser.add_argument(
        '--container','-c',type=str,
        help='''Provide a path to a SNPeff singularity container. Only needed for build mode.''',
        default='/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2025_funQCD_database/singularity_containers/snpeff_5-4c.sif'
        )    
    parser.add_argument(
        '--skip-cds',action='store_true',
        help='''Skip the CDS check when building the database.''',
        )
    parser.add_argument(
        '--skip-protein',action='store_true',
        help='''Skip the protein check when building the database.''',
        )
    args = parser.parse_args()
    if args.mode == 'prepare':
        prepare_snpeff_database(args.input)
    elif args.mode == 'build':
        if not os.path.isdir(os.path.join(args.input, "snpeff")):
            raise ValueError(f"SNPeff database directory not found in {args.input}. Please run in prepare mode first.")
        skip_line = ''
        if args.skip_cds:
            skip_line += '-noCheckCDS '
        if args.skip_protein:
            skip_line += '-noCheckProtein '
        fasta_name = os.path.basename(os.path.normpath(args.input))
        config_name = os.path.join(os.path.abspath(args.input), "snpeff", "snpEff.config")
        bind_string = f"{os.path.abspath(args.input)}:{os.path.abspath(args.input)}"
        os.chdir(os.path.join(args.input, "snpeff"))
        subprocess.run(['singularity', 'exec', '--bind', bind_string, args.container, 'snpEff', 'build', '-c', config_name, '-gff3', '-v', skip_line, os.path.basename(os.path.normpath(args.input))])

if __name__ == "__main__":
    main()


