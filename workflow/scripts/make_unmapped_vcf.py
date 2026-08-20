
import argparse
import gzip


def read_bed_intervals(bed_file):
    """Read zero-based, half-open BED intervals by contig."""
    intervals = {}
    with open(bed_file) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith(("#", "track", "browser")):
                continue
            line = line.strip()
            fields = line.split('\t')
            if len(fields) < 3:
                print(f'Error: incorrect format for BED file low-coverage intervals')
                quit(1)
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end <= start:
                print(f'Error: invalid BED interval in {bed_file}:{line_number}: {chrom}:{start}-{end}')
                quit(1)
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def overlaps_low_coverage(chrom, pos, ref, intervals):
    """
    Determine whether the reference span of a VCF record overlaps a BED
    interval.
    VCF POS is one-based. BED coordinates are zero-based and end-exclusive.
    """
    variant_start = pos - 1
    variant_end = variant_start + len(ref)

    return any(
        variant_start < bed_end and variant_end > bed_start
        for bed_start, bed_end in intervals.get(chrom, [])
    )


def add_filter(existing_filter, new_filter):
    """Append a FILTER label while preserving valid VCF semantics."""
    if existing_filter in (".", "PASS", ""):
        return new_filter
    
    labels = existing_filter.split(";")
    if new_filter not in labels:
        labels.append(new_filter)
    return ";".join(labels)


def flag_low_coverage_variants(input_vcf_gz,low_coverage_bed,output_vcf):
    intervals = read_bed_intervals(low_coverage_bed)
    filter_declared = False
    filter_index = None

    with gzip.open(input_vcf_gz, "rt") as input_handle, open(output_vcf, "w") as output_handle:
        for line in input_handle:
            if line.startswith("##"):
                if line.startswith(f"##FILTER=<ID=") and filter_declared is False:
                    output_handle.write(f"##FILTER=<ID=FAIL_LOW_COVERAGE,Description=\"Variant overlaps a reference interval below the minimum coverage threshold\">\n")
                    output_handle.write(line)
                    filter_declared = True
                else:
                    output_handle.write(line)

            elif line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                try:
                    filter_index = header_fields.index("FILTER")
                except ValueError as exc:
                    raise ValueError(
                        "VCF header does not contain a FILTER column"
                    ) from exc
                output_handle.write(line)

            elif line.startswith("#"):
                output_handle.write(line)

            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < len(header_fields):
                    raise ValueError(
                        "Incorrect number of columns in VCF record"
                    )
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                if overlaps_low_coverage(chrom, pos, ref, intervals):
                    fields[filter_index] = add_filter(existing_filter = fields[filter_index], new_filter = "FAIL_LOW_COVERAGE")

                output_handle.write("\t".join(fields) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add a VCF FILTER label to variants overlapping low-coverage BED intervals."
        )
    )
    parser.add_argument(
        "--vcf-file",
        required=True,
        help="Input bgzip- or gzip-compressed VCF"
    )
    parser.add_argument(
        "--low-cov-intervals",
        required=True,
        help="Three-column BED file of low-coverage intervals"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output uncompressed VCF"
    )

    args = parser.parse_args()
    
    flag_low_coverage_variants(args.vcf_file,args.low_cov_intervals,args.out)


if __name__ == "__main__":
    main()






















import gzip
import subprocess
import sys
import argparse

def mask_indel_prox_snps(vcf_gz_file, low_coverage_intervals_file, output_vcf_file):
    # read a vcf line-by-line
    # if the variant overlaps a low-coverage interval, append ",FAIL_LOW_COVERAGE" to the FILTER column
    # otherwise, write the line as-is to the output vcf file
    low_cov_intervals = {}
    with open(low_coverage_intervals_file) as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            low_cov_intervals[chrom] = low_cov_intervals.get(chrom, []) + [(int(start), int(end))]
    # retain all comment lines and structure by reading line by line
    with gzip.open(vcf_gz_file, 'rt') as f_in, open(output_vcf_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('##'):
                f_out.write(line)
            elif line.startswith('#CHROM'):
                # determine the index of the FILTER column
                header_fields = line.strip().split('\t')
                filter_index = header_fields.index('FILTER')
                info_index = header_fields.index('INFO')
                f_out.write(line)
            else:
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                filt = fields[filter_index]
                info = fields[info_index]
                # in order to account for indels, check both the start and end positions of the variant
                overlap_found = False
                if any(pos >= start and pos <= end for start,end in low_cov_intervals.get(chrom, [])):
                    overlap_found = True
                pos2 = pos + len(ref) - 1
                if pos2 != pos:
                    if any(pos2 >= start and pos2 <= end for start,end in low_cov_intervals.get(chrom, [])):
                        overlap_found = True
                # if an overlap is found, append ",FAIL_LOW_COVERAGE" to the original filter column
                if overlap_found:
                    line = line.replace(filt, filt + ',FAIL_LOW_COVERAGE')
                f_out.write(line)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf_file", required=True, help="VCF file to perform masking on. This should be a compressed file that ends in .vcf.gz")
    p.add_argument("--low_cov_intervals", required=True, help="File with low-coverage intervals. This should be a .bed file with 3 columns: contig_name, window_start, window_end")
    p.add_argument("--out", required=True, help="Output VCF filepath. This will be compressed and indexed after masking.")
    args = p.parse_args()
    mask_indel_prox_snps(args.vcf_file, args.low_cov_intervals, args.out)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr)
        raise
















# old version
# import argparse
# import subprocess
# import sys
# import tempfile

# def run(cmd):
#     return subprocess.run(cmd, check=True, text=True, capture_output=True)

# def fetch_ref_bases(ref_fasta, chrom, positions):
#     """
#     Fetch reference base for each position using samtools faidx in batch via -r regions_file.
#     Returns dict[pos] = base.
#     """
#     regions = [f"{chrom}:{p}-{p}" for p in positions]
#     with tempfile.NamedTemporaryFile(mode="w", delete=False) as tf:
#         for r in regions:
#             tf.write(r + "\n")
#         region_file = tf.name

#     try:
#         res = run(["samtools", "faidx", "-r", region_file, ref_fasta])
#     finally:
#         try:
#             import os
#             os.unlink(region_file)
#         except Exception:
#             pass

#     bases = {}
#     current_pos = None
#     for line in res.stdout.splitlines():
#         if not line:
#             continue
#         if line.startswith(">"):
#             hdr = line[1:].strip()              
#             chrom_part, rng = hdr.split(":")
#             start, end = rng.split("-")
#             current_pos = int(start)
#         else:
#             base = line.strip().upper()
#             if current_pos is not None:
#                 bases[current_pos] = base[0] if base else "N"
#                 current_pos = None
#     return bases

# def main():
#     p = argparse.ArgumentParser()
#     p.add_argument("--ref", required=True)
#     p.add_argument("--unmapped", required=True)
#     p.add_argument("--chrom", required=True, help="Reference contig name (e.g., contig_1)")
#     p.add_argument("--vcf_sample", required=True, help="Sample name to use in VCF header (must match PASS VCF sample)")
#     p.add_argument("--out", required=True, help="Output VCF (uncompressed .vcf)")
#     args = p.parse_args()

#     # read positions
#     positions = []
#     with open(args.unmapped) as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             positions.append(int(line))

#     positions = sorted(set(positions))
#     bases = fetch_ref_bases(args.ref, args.chrom, positions) if positions else {}

#     with open(args.out, "w") as out:
#         out.write("##fileformat=VCFv4.2\n")
#         out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
#         out.write('##FILTER=<ID=UNMAPPED,Description="Position flagged as unmapped/low coverage; masked in consensus">\n')
#         out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.vcf_sample}\n")

#         for pos in positions:
#             ref_base = bases.get(pos, "N")
#             # Mask as N and set GT to 1/1 (homozygous alt) to indicate masked position in consensus
#             out.write(f"{args.chrom}\t{pos}\t.\t{ref_base}\tN\t.\tUNMAPPED\t.\tGT\t1/1\n")

# if __name__ == "__main__":
#     try:
#         main()
#     except subprocess.CalledProcessError as e:
#         sys.stderr.write(e.stderr)
#         raise