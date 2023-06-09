#!/opt/patric-common/runtime/bin/python3
import csv
import glob 
import json 
import os
import sys
import shutil
import subprocess

def check_input_fastqs(input_dir, filename):
    input_path = f"{input_dir}/{filename}"
    if input_path.endswith('fastq.gz'):
        return input_path
    
    elif input_path.endswith('.fastq'):
        zipped_path = input_path + '.gz'
        # check if the zipped file already exists in the list of input files 
        if os.path.exists(zipped_path) == True:
            msg = "zipped file exists using {}.gz instead \n".format(input_path)
            sys.stderr.write(msg)
            return zipped_path

        else:
            msg = '{} identified as unzipped, zipping for analysis \n'.format(input_path)
            sys.stderr.write(msg)
            cmd = ["gzip", input_path]
            subprocess.run(cmd)
            return zipped_path

    else:
        msg = "Error {} not end in 'fastq' or 'fastq.gz'. \n {} ignored".format(input_path, input_path)
        sys.stderr.write(msg)
        pass

def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict

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

    msg = "using ftp to get hisat indicies for host {}".format(input_dict['host_genome'])
    sys.stderr.write(msg)
    # get hisat indicies
    curl_cmd = ['curl', '-o', 'hisat_indices_tar_file', 'ftp.bvbrc.org/host_genomes{}'.format(host_dict[input_dict['host_genome']])];
    subprocess.run(curl_cmd);
    # cmd 2
    tar_cmd = ['tar', '-xvf', 'hisat_indices_tar_file']
    subprocess.run(tar_cmd)
    return

# run the snakefile command
def run_snakefile(input_dict, config):
    # if there are paired end reads process them
    input_dir = config['input_data_dir']
    SNAKEMAKE_PATH = config['snakemake']
    ### relative path running from service-script ###
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'pe_fastq_processing')
        cmd = [SNAKEMAKE_PATH,
        "--snakefile",  SNAKEFILE,
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--debug"]
        subprocess.run(cmd)

    # if there are single end reads process them
    if os.path.exists('staging/se_reads'):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'se_fastq_processing')
        cmd = [SNAKEMAKE_PATH,
        "--snakefile",  SNAKEFILE,
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--debug"]
        subprocess.run(cmd)

    # paired and single end reads will be processed together
    if input_dict["analysis_type"] == "pathogen":
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'pathogen_analysis')
        cmd = [SNAKEMAKE_PATH,
        "--snakefile",  SNAKEFILE,
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--debug"]
        subprocess.run(cmd)

    if input_dict["analysis_type"] == 'microbiome':
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, 'microbiome_analysis')
        cmd = [SNAKEMAKE_PATH,
        "--snakefile",  SNAKEFILE,
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--debug"]
        subprocess.run(cmd)
    return


def set_up_sample_dictionary(input_dir, input_dict, output_dir):
    # set up the sample dictionary
    #### paired reads ####
    if len(input_dict['paired_end_libs']) != 0:
        paired_sample_dict = {}
        ws_paired_reads = []
        ws_paired_reads = input_dict['paired_end_libs']
        ### relative path running from service-script ###
        pe_path= f"{input_dir}/pe_reads"
        os.makedirs(pe_path, exist_ok = True)

        for i in range(len(ws_paired_reads)):
            read1_filename = ws_paired_reads[i]['read1'].split('/')[-1]
            read1_filepath = check_input_fastqs(input_dir, read1_filename)
            shutil.copy(read1_filepath, '/home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_taxonomic_classification_2/service-scripts/test')

            read2_filename = ws_paired_reads[i]['read2'].split('/')[-1]
            read2_filepath = check_input_fastqs(input_dir, read2_filename)
            sample_id = ws_paired_reads[i]['sample_id'].split('/')[-1]
            pe_r1_samplename = f"{sample_id}_R1.fastq.gz"
            pe_r2_samplename = f"{sample_id}_R2.fastq.gz"

            paired_sample_dict[read1_filename] = pe_r1_samplename
            paired_sample_dict[read2_filename] = pe_r2_samplename

            shutil.copy(read1_filepath, f"{input_dir}/pe_reads/{pe_r1_samplename}")
            shutil.copy(read2_filepath, f"{input_dir}/pe_reads/{pe_r2_samplename}")

    
    #### single reads ####
    if len(input_dict['single_end_libs']) != 0:
        single_end_sample_dict = {}
        ws_single_end_reads = []
        ws_single_end_reads = input_dict['single_end_libs']
        ### relative path running from service-script ###
        se_path= f"{input_dir}/se_reads"
        os.makedirs(se_path, exist_ok = True)

        for i in range(len(ws_single_end_reads)):
            se_filename = ws_single_end_reads[i]['read'].split('/')[-1]
            se_filepath = check_input_fastqs(input_dir, se_filename)
            sample_id = ws_single_end_reads[i]['sample_id']
            se_samplename = f"{sample_id}.fastq.gz"
            single_end_sample_dict[se_filename] = se_samplename

            shutil.copy(se_filepath, f"{input_dir}/se_reads/{se_samplename}")

    # export the sample dictionary to .CSV
    with open(f"{output_dir}/sample_key.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['User sample name', 'Analysis sample name'])
        if len(input_dict['paired_end_libs']) != 0:
            for key, value in paired_sample_dict.items():
                writer.writerow([key, value])
        if len(input_dict['single_end_libs']) != 0:
            for key, value in single_end_sample_dict.items():
                writer.writerow([key, value])
    return 


# run the script from service-script/app_taxonomic_classification perl script
# It takes a single argument which is the pathname of a config.json file
# This contains the app parameters in the params slot.
def main(argv):
    params_file = argv[0]
    try:
        fh = open(params_file)
        config = json.load(fh)
        fh.close()
    except IOError:
        print(f"Could not read params file {params_file}")
        sys.exit(1)
    except json.JasonDecodeError as e:
        print(f"JSON parse error in pyfile {params_file}: {e}")
        sys.exit(1)

    input_dict = config['params']
    input_dir = config['input_data_dir']
    output_dir = config['output_data_dir']
    set_up_sample_dictionary(input_dir, input_dict, output_dir)
    
    try:
        load_hisat_indicies(input_dict)
    except:
        pass
    run_snakefile(input_dict, config)


if __name__ == "__main__":
    main(sys.argv[1:])
