import os
import argparse
import subprocess

def make_dist_matrix(input_dir):
    # find the .msa file in the alignment subdirectory
    alignment_dir = os.path.join(input_dir, 'alignment')
    msa_fname = None
    for fname in os.listdir(alignment_dir):
        if fname.endswith('.fa'):
            msa_fname = os.path.join(alignment_dir, fname)
            break
    if msa_fname is None:
        print(f'Error: could not find the MSA file in {alignment_dir}')
        quit(1)
    # run snp-sites to reduce the msa to only variant sites
    fname = os.path.basename(msa_fname).split('_consensus_msa')[0]
    msa_var_fname = os.path.join(alignment_dir, f'{fname}_variant_only_alignment.fa')
    cmd1 = ['snp-sites', '-o', msa_var_fname, msa_fname]
    subprocess.run(cmd1)
    # also make a vcf of all variant positions
    vcf_var_fname = os.path.join(alignment_dir, f'{fname}_variant_positions.vcf')
    cmd1_5 = ['snp-sites', '-v', '-o', vcf_var_fname, msa_fname]
    subprocess.run(cmd1_5)
    # run snp-dists to generate the distance matrix
    dist_matrix_fname = os.path.join(alignment_dir, f'{fname}_distance_matrix.tsv')
    with open(dist_matrix_fname, 'w') as dist_matrix_file:
        cmd2 = ['snp-dists', msa_var_fname]
        subprocess.run(cmd2, stdout=dist_matrix_file)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a results directory from SNPmake. A distance matrix will be generated in the alignment subdirectory.''',
        required=True
        )
    args = parser.parse_args()
    make_dist_matrix(args.input)

if __name__ == "__main__":
    main()


