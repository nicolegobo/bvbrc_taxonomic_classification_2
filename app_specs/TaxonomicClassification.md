
# Application specification: TaxonomicClassification

This is the application specification for service with identifier TaxonomicClassification.

The backend script implementing the application is [App-TaxonomicClassification.pl](../service-scripts/App-TaxonomicClassification.pl).

The raw JSON file for this specification is [TaxonomicClassification.json](TaxonomicClassification.json).

This service performs the following task:   Compute taxonomic classification for read data

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| host_genome | Host Genome | enum  | :heavy_check_mark: | no_host |
| analysis_type | Analysis Type | enum  | :heavy_check_mark: | 16S |
| paired_end_libs |  | group  |  |  |
| single_end_libs |  | group  |  |  |
| srr_libs |  | group  |  |  |
| database | Database | enum  | :heavy_check_mark: | SILVA |
| save_classified_sequences | Save the classified sequences | bool  |  | 0 |
| save_unclassified_sequences | Save the unclassified sequences | bool  |  | 0 |
| confidence_interval |  |   |  | 0.1 |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

