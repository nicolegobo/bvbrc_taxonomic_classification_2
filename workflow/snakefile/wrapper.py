#!/opt/patric-common/runtime/bin/python3
import csv
import glob 
import json 
import os
import sys


def check_input_fastqs():
    input_paths =[]
    #### code block for checking the input fastqs ####
    ### relative path running from service-script ###
    input_paths = list(glob.glob("staging/*fastq*"))

    zipped_fastqs = []
    # iterate over every file 
    for input_path in input_paths:
        # find the extension
        if input_path.endswith('fastq.gz'):
            zipped_fastqs.append(input_path)
        elif input_path.endswith('.fastq'):
            zipped_path = input_path + '.gz'
            # check if the zipped file already exists in the list of input files 
            if os.path.exists(zipped_path) == True:
                msg = "ignoring {} duplicte file \n.".format(input_path)
                sys.stderr.write(msg)
                print(msg)
                pass
            else:
                msg = '{} identified as unzipped, zipping for analysis \n'.format(input_path)
                sys.stderr.write(msg)
                gzip_command = 'gzip {}'.format(input_path)
                os.system(gzip_command)
                # append to the list
                zipped_fastqs.append(zipped_path)
        else:
            msg = "Error {} not end in 'fastq' or 'fastq.gz'. \n {} ignored".format(input_path, input_path)
            sys.stderr.write(msg)
            pass
    return

def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict

# run the snakefile command
def run_snakefile(input_dict):
    # if there are paired end reads process them
    SNAKEMAKE_PATH = '/opt/patric-common/runtime/artic-ncov2019/bin/snakemake'
    ### relative path running from service-script ###
    SNAKEFILE_DIR = '../workflow/snakefile/'
    if os.path.exists('staging/pe_reads'):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'pe_fastq_processing')
        cmd = "{} \
            --snakefile {} \
            --use-singularity \
            --verbose \
            --printshellcmds \
            --debug \
            ".format(SNAKEMAKE_PATH, SNAKEFILE)
        os.system(cmd)

    # if there are single end reads process them
    if os.path.exists('staging/se_reads'):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'se_fastq_processing')
        cmd = "{} \
            --snakefile {} \
            --use-singularity \
            --verbose \
            --printshellcmds \
            --debug \
            ".format(SNAKEMAKE_PATH, SNAKEFILE)
        os.system(cmd)


    # paired and single end reads will be processed together
    if input_dict["analysis_type"] == "pathogen":
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'pathogen_analysis')

    if input_dict["analysis_type"] == 'microbiome':
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'microbiome_analysis')

    cmd = "{} \
        --snakefile {} \
        --use-singularity \
        --verbose \
        --printshellcmds \
        --debug \
        ".format(SNAKEMAKE_PATH, SNAKEFILE)
    os.system(cmd)
    return


def set_up_sample_dictionary(input_dict):
    # set up the sample dictionary
    #### paired reads ####
    paired_sample_dict = {}
    if len(input_dict['paired_end_libs']) != 0:
        ws_paired_reads = []
        ws_paired_reads = input_dict['paired_end_libs']
        ### relative path running from service-script ###
        pe_dir_command = 'mkdir -p staging/pe_reads'
        os.system(pe_dir_command)
        for i in range(len(ws_paired_reads)):
            read1_filename = ws_paired_reads[i]['read1'].split('/')[-1]
            read2_filename = ws_paired_reads[i]['read2'].split('/')[-1]
            pe_r1_samplename = "paired_sample_{}_R1.fastq.gz".format(i+1)
            pe_r2_samplename = "paired_sample_{}_R2.fastq.gz".format(i+1)

            paired_sample_dict[read1_filename] = pe_r1_samplename
            paired_sample_dict[read2_filename] = pe_r2_samplename
            ### commands for running from service-script ###
            cp_pe_r1_command = ('cp staging/{} staging/pe_reads/{}').format(read1_filename, paired_sample_dict[read1_filename])
            os.system(cp_pe_r1_command)
            cp_pe_r2_command = ('cp staging/{} staging/pe_reads/{}').format(read2_filename, paired_sample_dict[read2_filename])
            os.system(cp_pe_r2_command)
    
    #### single reads ####
    single_end_sample_dict = {}
    if len(input_dict['single_end_libs']) != 0:
        ### commands for running from service-script ###
        se_dir_command = 'mkdir -p staging/se_reads'
        os.system(se_dir_command)
        ws_single_end_reads = []
        ws_single_end_reads = input_dict['single_end_libs']
        for i in range(len(ws_single_end_reads)):
            se_filename = ws_single_end_reads[i]['read'].split('/')[-1]
            se_samplename = "single_sample_{}.fastq.gz".format(i+1)
            single_end_sample_dict[se_filename] = se_samplename
            ### for running from the service-script ###
            cp_se_command = ('cp staging/{} staging/se_reads/{}').format(se_filename, single_end_sample_dict[se_filename])
            os.system(cp_se_command)
    # export the sample dictionary to .CSV
    with open('output/sample_key.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['User sample name', 'Analysis sample name'])
        for key, value in paired_sample_dict.items():
            writer.writerow([key, value])
        for key, value in single_end_sample_dict.items():
            writer.writerow([key, value])
    return 

def load_hisat_indicies(input_dict):
    # get the refrence genome from parameters
    host_dict = {"homo_sapiens": "/9606/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.ht2.tar",
        "mus_musculus": "/10090/GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.ht2.tar",
        "rattus_norvegicus" : "/10116/Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.ht2.tar",
        "caenorhabditis_elegans" : "/6239/WBcel235/GCF_000002985.6_WBcel235_genomic.ht2.tar",
        "drosophila_melanogaster_strain" : "/7227/Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.ht2.tar",
        "danio_rerio_strain_tuebingen": "/7955/GRCz11/GCF_000002035.6_GRCz11_genomic.ht2.tar",
        "gallus_gallus" : "/9031/GRCg6a/GCF_000002315.6_GRCg6a_genomic.ht2.tar",
        "macaca_mulatta" : "/9544/Mmul_10/GCF_003339765.1_Mmul_10_genomic.ht2.tar",
        "mustela_putorius_furo" : "/9669/MusPutFur1.0/GCF_000215625.1_MusPutFur1.0_genomic.ht2.tar",
        "sus_scrofa" : "/9823/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.ht2.tar"}
    print(host_dict[input_dict['host_genome']])

    msg = "using ftp to get hisat indicies for host {}".format(input_dict['host_genome'])
    sys.stderr.write(msg)
    # get hisat indicies
    curl_cmd = 'curl -o hisat_indices_tar_file ftp.bvbrc.org/host_genomes{}'.format(host_dict[input_dict['host_genome']])
    print(curl_cmd)
    os.system(curl_cmd)
    # cmd 2
    tar_cmd = 'tar -xvf hisat_indices_tar_file'
    os.system(tar_cmd)
    return

# run the script from service-script/app_taxonomic_classification perl script
def main(argv):
    inputfile = argv
    raw_params = inputfile[0]
    input_dict = format_inputs(raw_params)
    check_input_fastqs()
    set_up_sample_dictionary(input_dict)
    load_hisat_indicies(input_dict)
    run_snakefile(input_dict)


if __name__ == "__main__":
    main(sys.argv[1:])