snakefile_dir = os.path.join(shapemake_dir, 'snakemake', 'snakefiles')

config = {
    'bioinformatics_dir': '/Users/andylavender/dev/shapemake_dev/bioinformatics',
    'rna_biology_dir': '/Users/andylavender/dev/shapemake_dev/rna_biology',

    'shapemake_dir': shapemake_dir,
    'env_dir': env_dir,
    'script_dir': script_dir,
    'snakefile_dir': snakefile_dir,

    'bcl_convert_path': 'bcl-convert',
    'bowtie2_build_path': 'bowtie2-build',
    'bowtie2_path': 'bowtie2',
    'fastqc_path': 'fastqc',
    'multiqc_path': 'multiqc',
    'partition_path': 'partition-smp',
    'probability_plot_path': 'ProbabilityPlot',
    'scale_nuc_path': 'scale_nuc.py',
    'shapely_cli_path': 'shapely/cli.py',
    'shapemapper_path': 'shapemapper',
    'snakemake_path': 'snakemake',
}