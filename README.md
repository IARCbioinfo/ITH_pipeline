# ITH_pipeline
## Nextflow pipeline to study intra-tumoral heterogeneity

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

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Tiffany Delhomme   |            delhommet@students.iarc | Developer to contact for support |

## References (optional)

## FAQ (optional)
