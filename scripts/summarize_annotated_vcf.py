
import argparse
import gzip

IMPACT_RANK = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0,
}

def open_vcf(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "r")

def is_missing_vcf_format(input_string):
    # determine if the input string consists only of ':' and '.', the characters used for empty entries in a vcf
    return set(input_string).issubset(set(':.'))

def summarize_vcf(vcf_file,output_file):
    # read in a vcf line-by-line
    # return a tsv that shows which isolates have which variant, and if the variant has a notable consequence from annotation
    with open_vcf(vcf_file) as fh, open(output_file,'w') as fh_out:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                # these are the column names
                # assume that the first nine columns are CHROM, POS, etc through FORMAT
                line = line.rstrip().split('\t')
                if line[0] != '#CHROM' or line[1] != 'POS' or line[8] != 'FORMAT':
                    raise ValueError('Unexpected VCF column headers')
                isolate_names = line[9:]
                _ = fh_out.write('CHROM\tPOS\tREF\tALT\tANN_GENE\tANN_DESC\tANN_CATEGORY\tANN_NOTE\t' + '\t'.join(isolate_names) + '\n')
            else:
                line = line.rstrip().split('\t')
                CHROM = line[0]
                POS = line[1]
                REF = line[3]
                ALT = line[4]
                INFO = line[7]
                ISOL_DATA = line[9:]
                ISOL_DATA = ['0' if is_missing_vcf_format(x) else '1' for x in ISOL_DATA]
                # take only the ANN field from format
                ANN = INFO.split('ANN=')[-1]
                ANN = ANN.split(',')
                ANN_DATA = ANN[0].split('|')
                for new_ann_data in ANN:
                    new_ann_data_list = new_ann_data.split('|')
                    if IMPACT_RANK[new_ann_data_list[2]] > IMPACT_RANK[ANN_DATA[2]]:
                        ANN_DATA = new_ann_data_list
                ANN_GENE = ANN_DATA[3]
                ANN_DESC = ANN_DATA[1]
                ANN_CATEGORY = ANN_DATA[2]
                ANN_NOTE = ANN_DATA[-1]
                if ANN_NOTE == '':
                    ANN_NOTE = 'NA'
                towrite = '\t'.join([CHROM,POS,REF,ALT,ANN_GENE,ANN_DESC,ANN_CATEGORY,ANN_NOTE] + ISOL_DATA) + '\n'
                _ = fh_out.write(towrite)
                
                

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to the input vcf.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output summary .tsv file.''',
        required=True
        )
    args = parser.parse_args()
    summarize_vcf(args.input,args.output)
    
    

if __name__ == "__main__":
    main()




