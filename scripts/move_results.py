import os
import subprocess
import argparse


def move_variant_calling(source_dir,batch_name,target_dir):
    src_dir = os.path.join(source_dir,batch_name)
    dest_dir = os.path.join(target_dir,batch_name)
    if os.path.isdir(dest_dir):
        print(f'Warning: {dest_dir} already exists!')
    os.makedirs(dest_dir,exist_ok=True)

    alignment_dir_src = os.path.join(src_dir,'alignment')
    subprocess.run(['cp','-r',alignment_dir_src,dest_dir])

    con_dir_src = os.path.join(src_dir,'consensus')
    con_dir_dst = os.path.join(dest_dir,'consensus')
    os.makedirs(con_dir_dst,exist_ok=True)
    for fname in os.listdir(con_dir_src):
        con_dir_src_file = os.path.join(con_dir_src,fname)
        if not os.path.isdir(con_dir_src_file):
            subprocess.run(['cp',con_dir_src_file,con_dir_dst])

    merged_dir_src = os.path.join(src_dir,'merged_vcf')
    merged_dir_dst = os.path.join(dest_dir,'merged_vcf')
    os.makedirs(merged_dir_dst,exist_ok=True)
    for fname in os.listdir(merged_dir_src):
        if fname.endswith(('_annotated.vcf.gz','_annotated.vcf.gz.tbi','_only.vcf.gz','_only.vcf.gz.tbi')):
            subprocess.run(['cp',os.path.join(merged_dir_src,fname),merged_dir_dst])

    filtered_dir_src = os.path.join(src_dir,'filtered_vcf')
    filtered_dir_dst = os.path.join(dest_dir,'filtered_vcf')
    os.makedirs(filtered_dir_dst,exist_ok=True)
    for fname in os.listdir(filtered_dir_src):
        filtered_sample_src_dir = os.path.join(filtered_dir_src,fname)
        if os.path.isdir(filtered_sample_src_dir):
            filtered_sample_dst_dir = os.path.join(filtered_dir_dst,fname)
            os.makedirs(filtered_sample_dst_dir,exist_ok=True)
            for sub_fname in os.listdir(filtered_sample_src_dir):
                if sub_fname.endswith(('_lowcov.vcf.gz','_indelprox_masked_positions.txt')):
                    subprocess.run(['cp',os.path.join(filtered_sample_src_dir,sub_fname),filtered_sample_dst_dir])

    masked_dir_src = os.path.join(src_dir,'bedtools')
    for fname in os.listdir(masked_dir_src):
        masked_sample_src_dir = os.path.join(masked_dir_src,fname)
        if os.path.isdir(masked_sample_src_dir):
            filtered_sample_dst_dir = os.path.join(filtered_dir_dst,fname)
            os.makedirs(filtered_sample_dst_dir,exist_ok=True)
            for sub_fname in os.listdir(masked_sample_src_dir):
                if sub_fname.endswith(('_final_mask.bed.gz','_GenomeCoverage.bedgraph')):
                    subprocess.run(['cp',os.path.join(masked_sample_src_dir,sub_fname),filtered_sample_dst_dir])


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--source_dir','-s',type=str,
        help='''Provide a source directory. This should be the results/ directory from SNPmake.''',
        required=True
        )
    parser.add_argument(
        '--dest_dir','-d',type=str,
        help='''Provide a destination directory.''',
        required=True
        )
    parser.add_argument(
        '--batch_name','-b',type=str,
        help='''Provide the batch name. This should be the name of directory inside results/ that you want to move.''',
        required=True
        )
    args = parser.parse_args()
    move_variant_calling(args.source_dir,args.batch_name,args.dest_dir)

if __name__ == "__main__":
    main()
