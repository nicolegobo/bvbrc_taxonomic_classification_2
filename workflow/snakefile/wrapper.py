#!/opt/patric-common/runtime/bin/python3
import csv
import glob 
import json 
import os
import re
import sys
import shutil
import subprocess

def check_input_fastqs(input_dir, filename, files_to_zip):
    input_path = f"{input_dir}/{filename}"
    if input_path.endswith("fastq.gz") or input_path.endswith("fq.gz"):
        return input_path
    
    elif input_path.endswith(".fastq") or input_path.endswith(".fq"):
        zipped_path = input_path + ".gz"
        # check if the zipped file already exists in the list of input files 
        if os.path.exists(zipped_path) == True:
            msg =f"zipped file exists using {input_path}.gz instead \n"
            sys.stderr.write(msg)
            return zipped_path

        else:
            msg = f"{input_path} identified as unzipped, zipping for analysis \n"
            sys.stderr.write(msg)
            files_to_zip.append(input_path)
            return zipped_path
    else:
        msg = f"Error {input_path} not end in 'fastq' or 'fastq.gz'. \n {input_path} ignored \n"
        sys.stderr.write(msg)
        sys.exit(1)
        pass


def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict


def load_hisat_indicies(input_dict):
    # get the refrence genome from parameters
    if input_dict["host_genome"] == "no_host":
        msg = "No host genome selected. Not getting hisat indices \n"
        sys.stderr.write(msg)
    else:
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

        msg = "using ftp to get hisat indicies for host {} \n".format(input_dict["host_genome"])
        sys.stderr.write(msg)
        # get hisat indicies
        curl_cmd = ["curl", "-o", "hisat_indices_tar_file", "ftp.bvbrc.org/host_genomes{}".format(host_dict[input_dict["host_genome"]])];
        subprocess.run(curl_cmd);
        # cmd 2
        tar_cmd = ["tar", "-xvf", "hisat_indices_tar_file"]
        subprocess.run(tar_cmd)
    return


def post_processing_check(all_sample_ids, output_dir):
    dict_samples = {}
    complete = []
    incomplete = []
    for sample_name in all_sample_ids:
        # check for problematic Kraken results
        krona_path = f"{output_dir}/{sample_name}_sankey.html"
        if os.path.isfile(krona_path) == True:
            msg = f"{sample_name}"
            complete.append(msg)
        else:
            dict_samples[sample_name] = False
            msg = f"{sample_name}"
            incomplete.append(msg)
    ## write the complete and incomplete samples to the error report for user to see ##
    if len(incomplete) != 0:
        msg = f"Reliable kraken results produced for the following samples: {complete}. \n \
                Review the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
        sys.exit(1)
    else:
        msg = f"Reliable kraken results produced for the following samples: {complete}. \n \
        Post Kraken processing complete"
        sys.stderr.write(msg)
        return True


def preprocessing_check(input_dir, output_dir, input_dict):
    ## Using the input dictionary instead of path from staging directory due to file name confusion
    # Parse "sample_id" from "paired_end_libs"
    paired_sample_ids = [item["sample_id"] for item in input_dict.get("paired_end_libs", [])]

    # Parse "sample_id" from "single_end_libs"
    single_sample_ids = [item["sample_id"] for item in input_dict.get("single_end_libs", [])]

    # Parse "sample_id" from "srr_libs"
    srr_sample_ids = [item["sample_id"] for item in input_dict.get("srr_libs", [])]
    # Merge all sample IDs into one list
    all_sample_ids = paired_sample_ids + single_sample_ids + srr_sample_ids

    # Edit the sample ids to match the sample ids defined in set-up-sample-dictionary(input_dir input_dict output_dir cores)
    clean_sample_ids = []
    for sample_id in all_sample_ids:
        # Define a regular expression pattern to match all special characters except underscore
        pattern = r"[^a-zA-Z0-9_]"
        # Use the re.sub() function to replace all matches of the pattern with an empty string
        sample_id = re.sub(pattern, "", sample_id)
        clean_sample_ids.append(sample_id)
    all_sample_ids = clean_sample_ids
        
    # file check for preprocessing/kraken run
    dict_samples = {}
    complete = []
    incomplete = []

    # check if Kraken report exists for each sample ID
    for sample_name in all_sample_ids:
        report_path = f"{output_dir}/{sample_name}/kraken_output/{sample_name}_k2_output.txt"
        if os.path.isfile(report_path) == True:
            complete.append(sample_name)
        else:
            dict_samples[sample_name] = False
            incomplete.append(sample_name)
    if len(incomplete) != 0:
        msg = f"Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing \n \
                FASTQ Processing is complete for the following samples: {complete}. \n \
                FASTQ Processing is INCOMPLETE for the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
        sys.exit(1)
        ## return False
    else:
        msg = "preprocessing complete \n"
        sys.stderr.write(msg)
        return True, all_sample_ids

def run_16s_snakefile(input_dict, input_dir, output_dir,  config):
    # process any paired reads
    input_dir = config["input_data_dir"]
    SNAKEMAKE_PATH = config["snakemake"]
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"

    common_params = [
        SNAKEMAKE_PATH,
        "--cores", str(config["cores"]),
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going"
        ]
    if config["cores"] == 1:
        common_params.append("--debug")
        
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_fastq_processing_16s")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)

    # process any single end reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_fastq_processing_16s")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
    kraken_check_result = preprocessing_check(input_dir, output_dir, input_dict)
    kraken_check = kraken_check_result[0]
    all_sample_ids = kraken_check_result[1]
    
    if  kraken_check == True:
        # if all of the kraken report files exist then the final step will trigger.
        # paired and single end reads will be processed together
        if input_dict["analysis_type"] == "16S":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "16s_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            subprocess.run(cmd)
    # file check if anything is wrong with the kraken outputs it should fail to 
    # produce krona plot
    krona_check = post_processing_check(all_sample_ids, output_dir)
    if krona_check == True:
        msg = "16s analysis is complete \n"
        sys.stderr.write(msg)
    return

# run the snakefile command for WGS
def run_wgs_snakefile(input_dict, input_dir, output_dir,  config):

    input_dir = config["input_data_dir"]
    SNAKEMAKE_PATH = config["snakemake"]
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"

    common_params = [
        SNAKEMAKE_PATH,
        "--cores", str(config["cores"]),
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going"
        ]

    if config["cores"] == 1:
        common_params.append("--debug")
    
    # process any paired reads
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_fastq_processing")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)

    # process any single reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_fastq_processing")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        subprocess.run(cmd)
    kraken_check_result = preprocessing_check(input_dir, output_dir, input_dict)
    ### check!
    kraken_check = kraken_check_result[0]
    all_sample_ids = kraken_check_result[1]

    if  kraken_check == True:
        # if all of the kraken report files exist then the final step will trigger.
        # paired and single end reads will be processed together
        if input_dict["analysis_type"] == "pathogen":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pathogen_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            subprocess.run(cmd)

        if input_dict["analysis_type"] == "microbiome":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "microbiome_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            subprocess.run(cmd)
    # final file check
    krona_check = post_processing_check(all_sample_ids, output_dir)
    if krona_check == True:
        msg = " WGS analysis is complete \n"
        sys.stderr.write(msg)
    return


def set_up_sample_dictionary(input_dir, input_dict, output_dir, cores):
    # set up the sample dictionary
    #### paired reads ####
    files_to_zip = []
    to_copy = []
    #sample_ids = []

    if len(input_dict["paired_end_libs"]) != 0:
        paired_sample_dict = {}
        ws_paired_reads = []
        ws_paired_reads = input_dict["paired_end_libs"]

        pe_path= f"{input_dir}/pe_reads"
        os.makedirs(pe_path, exist_ok = True)

        for i in range(len(ws_paired_reads)):
            read1_filename = ws_paired_reads[i]["read1"].split("/")[-1]
            read1_filepath = check_input_fastqs(input_dir, read1_filename, files_to_zip)

            read2_filename = ws_paired_reads[i]["read2"].split("/")[-1]
            read2_filepath = check_input_fastqs(input_dir, read2_filename, files_to_zip)
            sample_id = ws_paired_reads[i]["sample_id"].split("/")[-1]
            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
            pe_r1_samplename = f"{sample_id}_R1.fastq.gz"
            pe_r2_samplename = f"{sample_id}_R2.fastq.gz"
            #sample_ids.append(sample_id)
            paired_sample_dict[read1_filename] = pe_r1_samplename
            paired_sample_dict[read2_filename] = pe_r2_samplename

            to_copy.append([read1_filepath, f"{input_dir}/pe_reads/{pe_r1_samplename}"])
            to_copy.append([read2_filepath, f"{input_dir}/pe_reads/{pe_r2_samplename}"])


    #### single reads ####
    if len(input_dict["single_end_libs"]) != 0:
        single_end_sample_dict = {}
        ws_single_end_reads = []
        ws_single_end_reads = input_dict["single_end_libs"]
        ### relative path running from service-script ###
        se_path= f"{input_dir}/se_reads"
        os.makedirs(se_path, exist_ok = True)

        for i in range(len(ws_single_end_reads)):
            se_filename = ws_single_end_reads[i]["read"].split("/")[-1]
            se_filepath = check_input_fastqs(input_dir, se_filename, files_to_zip)
            sample_id = ws_single_end_reads[i]["sample_id"]
            #sample_ids.append(sample_id)

            # Define a regular expression pattern to match all special characters except underscore
            pattern = r"[^a-zA-Z0-9_]"
            # Use the re.sub() function to replace all matches of the pattern with an empty string
            sample_id = re.sub(pattern, "", sample_id)
            se_samplename = f"{sample_id}.fastq.gz"
            single_end_sample_dict[se_filename] = se_samplename
            to_copy.append([se_filepath, f"{input_dir}/se_reads/{se_samplename}"])

    #
    # in parallel, gzip anything that needs to be gzipped
    # We deferred the copy until after the zip
    #
    par_inp = "\n".join(files_to_zip + [""])
    subprocess.run(["parallel", "-j", str(cores), "gzip", "{}"], input=par_inp,text=True)

    for (src, dest) in to_copy:
        shutil.copy(src, dest)

    # export the sample dictionary to .CSV
    with open(f"{output_dir}/sample_key.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["User sample name", "Analysis sample name"])
        if len(input_dict["paired_end_libs"]) != 0:
            for key, value in paired_sample_dict.items():
                writer.writerow([key, value])
        if len(input_dict["single_end_libs"]) != 0:
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

    input_dict = config["params"]
    input_dir = config["input_data_dir"]
    output_dir = config["output_data_dir"]

    set_up_sample_dictionary(input_dir, input_dict, output_dir, min(8, int(config["cores"])))
    try:
        load_hisat_indicies(input_dict)
    except:
        pass
    if input_dict["analysis_type"] == "pathogen" or input_dict["analysis_type"] == "microbiome":
        run_wgs_snakefile(input_dict, input_dir, output_dir, config)
    elif input_dict["analysis_type"] == "16S":
        run_16s_snakefile(input_dict, input_dir, output_dir,  config)
    else:
        msg = f"unsuported analysis type {input_dict['analysis_type']} \n"
        sys.stderr.write(msg)
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
