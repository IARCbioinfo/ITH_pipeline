# ITH_pipeline
## Nextflow pipeline to study intra-tumoral heterogeneity

[comment]: <> [![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/template-nf)
[comment]: <> [![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)
[comment]: <> [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1404)
[comment]: <> [![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

[comment]: <> ![Workflow representation](template-nf.png)

## Description
Nextflow pipeline to study intra-tumoral heterogeneity through subclonality reconstruction, using HATCHet, DeCiFer and clonEvol.
PIPELINE IN DEVELOPMENT.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [HATCHet](https://github.com/raphael-group/hatchet)
- [DeCiFer](https://github.com/raphael-group/decifer)
- [ClonEvol](https://github.com/hdng/clonevol)

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input
  | Type      | Description     |
  |-----------|---------------|
  | --bam_folder    | Folder containing BAM files to pass to HATCHet |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --ref    |            /whole/path/to/genome.fa | Reference genome in fasta format |
| --correspondance    |            TN_pairs.txt | Tabulated file containing two columns: tumor and normal bam names |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --cpu   |            1 | Number of cpu to use |
| --config    |            null | Use custom configuration file |
| --mem   |            4 | Size of memory used in GB |
| --output_folder   |            . | Path to output folder  |
| --samtools_folder   |            /usr/bin/ | samtools installation dir |
| --bcftools_folder   |            /usr/bin/ | bcftools installation dir |
| --hatchet_folder   |            /usr/bin/ | hatchet installation dir |
| --bnpy_folder   |            /usr/bin/ | bnpy installation dir |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run iarcbioinfo/ITH_pipeline --bam_folder bam_folder/ --ref genome.fa --correspondance pairs.txt --cpu 24 --samtools_folder ~/bin/samtools-1.7 --bcftools_folder ~/bin/bcftools-1.7 --hatchet_folder ~/bin/hatchet --bnpy_folder ~/bin/bnpy-dev
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | ......    | ...... |


## Detailed description
...

## Directed Acyclic Graph
[comment]: <> [![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/template-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Tiffany Delhomme   |            delhommet@students.iarc | Developer to contact for support |

## References (optional)

## FAQ (optional)
