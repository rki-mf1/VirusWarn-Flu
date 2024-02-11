<div id="top"></div>

<div align="center">
<h1 align="center"> FluWarnSystem: Alert system and prioritization for concerning Influenza variants </h1>
</div>
The goal of VOCAL-Flu is to detect concerning Influenza variants from sequencing data.
It does so by parsing Influenza genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. 
The tool is based on <a href="https://github.com/rki-mf1/vocal"><strong>VOCAL</strong></a>: Variant Of Concern ALert and prioritization, which was invented for SARS-CoV-2.

# Getting Started

⚠️**Note**: 🔌 Right now, FluWarnSystem is tested on Linux and Mac system only 💻 

## Quick Installation

To run the pipeline, you need to have `conda` and `Nextflow` installed and set up.
All other dependencies will be installed over `conda` in the pipeline.

To install `conda`, use the following bash commands:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Then, `Nextflow` an be installed over `conda`:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```

The usher repository can be cloned from Git:
```bash
git clone https://github.com/rki-mf1/vocal-flu.git
```

### Call help

```bash
nextflow run main.nf --help
```

## Running FluWarnSystem

```bash
cd vocal-flu
conda activate nextflow
nextflow run main.nf \
     -profile conda,local \
     --fasta 'path/to/input/fasta' 
```

### Running FluWarnSystem with splitting (Input from FluPipe)

```bash
cd vocal-flu
conda activate nextflow
nextflow run main.nf \
     -profile conda,local \
     --fasta 'path/to/input/fasta' \
     --split 'FluPipe'
```

### Running FluWarnSystem with PSL alignment

```bash
cd vocal-flu
conda activate nextflow
nextflow run main.nf \
     -profile conda,local \
     --fasta 'path/to/input/fasta' \
     --psl true
```

## Parameter list

```
fasta                    REQUIRED! The path to the fasta file with the sequences for FluWarnSystem.
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
```


# Contact

Did you find a bug?🐛 Suggestion/Feedback/Feature request?👨‍💻
Please visit [GitHub Issues](https://github.com/rki-mf1/vocal-flu/issues)

For business inquiries or professional support requests 🍺
Please contact Dr. Hölzer, Martin(<HoelzerM@rki.de>) or Dr. Richard, Hugues (<RichardH@rki.de>)


# Acknowledgments

* Original tool: [VOCAL](https://github.com/rki-mf1/vocal): Variant Of Concern ALert and prioritization 

    * Original Idea: SC2 Evolution Working group at the Robert Koch Institute in Berlin

    * Funding: Supported by the European Centers for Disease Control [grant number ECDC GRANT/2021/008 ECD.12222].


# Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.