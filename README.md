<div id="top"></div>

<div align="center">
<h1 align="center"> VirusWarn-Flu </h1>
<h3 align="center"> Mutation-Based Early Warning System to Prioritize Concerning Influenza Variants from Sequencing Data </h3>
</div>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

The goal of VirusWarn-Flu is to detect concerning Influenza variants from sequencing data.
It does so by parsing Influenza genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. 
The tool is based on <a href="https://github.com/rki-mf1/VirusWarn-SC2"><strong>VirusWarn-SC2</strong></a>, which was invented for SARS-CoV-2.


# Documentation

VirusWarn-Flu is part of *VirusWarn*

<a href="https://rki-mf1.github.io/viruswarn-doc/"><strong>Explore »</strong></a>


# Getting Started

⚠️ **Note**: 🔌 Right now, VirusWarn-Flu is tested on Linux and Mac system only 💻 

## Quick Installation

To run the pipeline, you need to have `Nextflow` and either `conda`, `Docker` or `Singularity`.

<details><summary><strong>Click!</strong> If you want to install <code>Nextflow</code> directly, you can use the following one-liner. </summary>

```bash
wget -qO- https://get.nextflow.io | bash
```
</details>

<details><summary><strong>Click!</strong> If you want to set up <code>conda</code> to run the pipeline and install all other dependencies through it, you can use the following steps. </summary>

Use the following bash commands if you are working on **Linux**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Use the following bash commands if you are working on **Mac**:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

Then, `Nextflow` an be installed over `conda`:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```
</details>

### Get / Update VirusWarn-Flu

```bash
nextflow pull rki-mf1/VirusWarn-Flu
```

### Call help

```bash
nextflow run rki-mf1/VirusWarn-Flu -r <version> --help
```

## Running VirusWarn-Flu

Once nextflow, is installed you are good to go! VirusWarn-Flu is run with only one command, using either conda (or mamba), Docker or Singularity. 

With a `conda`, please run:

```bash
nextflow run rki-mf1/VirusWarn-Flu -r <version> \
     -profile conda,local \
     --fasta 'test/openflu_h1n1.fasta' \
     --metadata 'test/metadata_h1n1.xlsx'
```

With a `Docker`, please run:

```bash
nextflow run rki-mf1/VirusWarn-Flu -r <version> \
     -profile docker,local \
     --fasta 'test/openflu_h1n1.fasta' \
     --metadata 'test/metadata_h1n1.xlsx'
```

With a `Singularity`, please run:

```bash
nextflow run rki-mf1/VirusWarn-Flu -r <version> \
     -profile singularity,local \
     --fasta 'test/openflu_h1n1.fasta' \
     --metadata 'test/metadata_h1n1.xlsx'
```

### Running VirusWarn-Flu with splitting (Input from OpenFlu)

```bash
nextflow run rki-mf1/VirusWarn-Flu -r <version> \
     -profile conda,local \
     --fasta 'test/openflu_combi.fasta' \
     --metadata 'test/metadata_combi.xlsx' \
     --split 'OpenFlu'
```

## Parameter list

```
fasta                    REQUIRED! The path to the fasta file with the sequences for VirusWarn-Flu.
                         [ default: '' ]
ref                      If you want to use the recent references from Nextclade, choose ''.
                         H1N1: A/Wisconsin/588/2019 (MW626065)
                         H3N2: A/Darwin/6/2021 (EPI1857216)
                         If you want to use the older references for H1N1 and H3N2, choose 'old'.
                         H1N1: A/California/7/2009 (CY121680)
                         H3N2: A/Wisconsin/67/2005 (CY163680)
                         For Victoria, only B/Brisbane/60/2008 (KX058884) is available.
                         [ default: '' ]
metadata                 The path to a metadate file for the sequences with collection dates.
                         Required to generate a heatmap in the report.
                         [ default: '' ]
subtype                  If the input fasta file only contains sequences of one subtype, 
                         define the subtype to choose the right references and tables.
                         The options are H1N1 and H3N2 for Influenza A,
                         Victoria for Influenza B.
                         [ default: 'h1n1' ]
split                    If the input fasta file contains sequences of more than one subtype, 
                         enable the split parameter to write them into one file per subtype and 
                         ensure the use of the right references and tables.
                         The options are FluPipe, GISAID and OpenFlu.
                         [ default: '' ]
qc                       If set to true, a QC report will be generated from the Nextclade output.
                         [ default: true ]
strict                   Run process with strict alert levels (without orange). Choose 'y'.
                         [ default: 'n' ]
season                   The Influenza season from which the input sequences are.
                         Important for checking on substitutions that are fixed in the population.
                         [ default: '23/24' ]
```


# Data

For further information on the tables that are used for the ranking, please take a look at the subfolders [`A(H1N1)pdm09`](data/A(H1N1)pdm09/), [`A(H3N2)`](data/A(H3N2)/) and [`B(Victoria)`](data/B(Victoria)/) of the folder [`data`](data/) depending on the subtype you are interested in and the READMEs that are provided there.

For further instructions for test runs and information on the used data, please take a look at the folder [`test`](test/) and the README that is provided there.

An example of the HTML report can be found in the folder [`example`](example/).


# How to interprete result.

VirusWarn-Flu output an alert level in four different colours which can be classified into 3 ratings.

| Alert color | Description |      Impact | 
| ----------- | ----------- | ----------- |
| Pink | Variants with fixed MOCs and ROIs from the previous season. | HIGH |
| Red | Variants with a high number of MOCs and ROIs that can be dangerous.     | HIGH |
| Orange | Variants in the orange level have less MOC and ROI than the ones in the red level and are therefore considered less dangerous but still concerning.   | MODERATE |
| Yellow | Variants that accumulate a high number of ROIs or PMs are sorted in the pink level for further inspection.   | MODERATE |
| Grey | The remaining variants are assigned to the black category.             | LOW |


# Contact

Did you find a bug? 🐛 Suggestion/Feedback/Feature request? 👨‍💻
Please visit [GitHub Issues](https://github.com/rki-mf1/VirusWarn-Flu/issues)

For business inquiries or professional support requests 🍺
Please feel free to contact us!


# Acknowledgments

* Original tool: [VirusWarn-SC2](https://github.com/rki-mf1/VirusWarn-SC2) (former VOCAL - Variant Of Concern ALert and prioritization)

    * Original Idea: SC2 Evolution Working group at the Robert Koch Institute in Berlin

    * Funding: Supported by the European Centers for Disease Control [grant number ECDC GRANT/2021/008 ECD.12222].


# Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
