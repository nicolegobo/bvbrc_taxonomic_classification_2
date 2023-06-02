# Taxonomic Classification

The BV-BRC taxonomic classification application.


The perl wrapper creates a config.json that drives the python snakemake-invocation code.

That config.json looks like this:

```
{
   "input_data_dir" : "/homes/olson/P3/dev-slurm/dev_container/modules/bvbrc_taxonomic_classification_2/staging",
   "output_data_dir" : "/homes/olson/P3/dev-slurm/dev_container/modules/bvbrc_taxonomic_classification_2/output",
   "params" : {
      "algorithm" : "Kraken2",
      "analysis_type" : "pathogen",
      "database" : "standard",
      "input_type" : "reads",
      "output_file" : "species_identification_example",
      "output_path" : "/olson@patricbrc.org/home/test",
      "paired_end_libs" : [
         {
            "read1" : "/olson@patricbrc.org/PATRIC-QA/applications/App-TaxonomicClassification/data/SRR12486981_1.fastq.gz",
            "read2" : "/olson@patricbrc.org/PATRIC-QA/applications/App-TaxonomicClassification/data/SRR12486981_2.fastq.gz"
         },
         {
            "read1" : "/olson@patricbrc.org/PATRIC-QA/applications/App-TaxonomicClassification/data/SRR12486983_1.fastq.gz",
            "read2" : "/olson@patricbrc.org/PATRIC-QA/applications/App-TaxonomicClassification/data/SRR12486983_2.fastq.gz"
         }
      ],
      "single_end_libs" : [],
      "srr_ids" : [],
      "srr_libs" : []
   },
   "snakemake" : "/disks/patric-common/runtime/artic-ncov2019/bin/snakemake",
   "workflow_dir" : "/home/olson/P3/dev-slurm/dev_container/modules/bvbrc_taxonomic_classification_2/workflow"
}
```