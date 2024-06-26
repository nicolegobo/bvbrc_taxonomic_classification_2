# Taxonomic Classification

## Overview

The Taxonomic Classification Service accepts reads or SRR values from sequencing of a metagenomic sample and uses [Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) to assign the reads to taxonomic bins, providing an initial profile of the possible constituent organisms present in the sample. We support taxonomic classification for whole genome sequencing data (WGS) and for 16s rRNA sequencing. It is important that you select the sequence type. Then the analysis options and database options will change support your sequence type.

## About this Module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There is one application service specification defined here:
1.  [Taxonomic Classification](app_specs/TaxonomicClassification.md): Service that that provides the backend for the BV-BRC web inerface; it takes reads as input.

The code in this module provides the BV-BRC application service wrapper scripts for the genome annotation service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-TaxonomicClassification.pl](service-scripts/App-TaxonomicClassification.pl) | App script for the [taxonomic classification service](https://www.bv-brc.org/docs/quick_references/services/taxonomic_classification_service.html) |

## See also 
* [Taxonomic Classification Service](https://www.bv-brc.org/app/TaxonomicClassification)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/taxonomic_classification_service.html)
* [Taxonomic Classification Service Tutorial](https://www.bv-brc.org/docs/tutorial/taxonomic_classification/taxonomic_classification.html)

## References
Andrews, S. FastQC: a quality control tool for high throughput sequence data. (2010). http://www.bioinformatics.babraham.ac.uk/projects/fastqc

Lu, J. et al. Metagenome analysis using the Kraken software suite. Nature protocols 17, 2815-2839 (2022).

Kim, D., Paggi, J. M., Park, C., Bennett, C. & Salzberg, S. L. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature biotechnology 37, 907-915 (2019).

Wood, D. E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome biology 20, 257 (2019).

Ondov, B. D., Bergman, N. H. & Phillippy, A. M. Interactive metagenomic visualization in a Web browser. BMC bioinformatics 12, 385 (2011).

Breitwieser, F. P. & Salzberg, S. L. Pavian: Interactive analysis of metagenomics data for microbiomics and pathogen identification. bioRxiv, 084715, doi:10.1101/084715 (2016).

Lu, J., Breitwieser, F. P., Thielen, P. & Salzberg, S. L. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3, e104 (2017).

Ewels, P., Magnusson, M., Lundin, S. & Käller, M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047-3048 (2016).

Krueger, F. Trim Galore: a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries. URL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. (Date of access: 28/04/2016) (2012).

O’Leary, N. A. et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic acids research 44, D733-D745 (2016).

Quast, C. et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic acids research 41, D590-D596 (2012).

DeSantis, T. Z. et al. Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. Applied and environmental microbiology 72, 5069-5072 (2006).

Anderson, M. J., Ellingsen, K. E. & McArdle, B. H. Multivariate dispersion as a measure of beta diversity. Ecology letters 9, 683-693 (2006).

Bray, J. R. & Curtis, J. T. An ordination of the upland forest communities of southern Wisconsin. Ecological monographs 27, 326-349 (1957).
