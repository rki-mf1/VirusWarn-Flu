<div id="top"></div>

<div align="center">
<h1 align="center"> VOCAL-Flu: Alert system and prioritization for concerning Influenza variants </h1>
</div>
The goal of VOCAL-Flu is to detect concerning Influenza variants from sequencing data.
It does so by parsing Influenza genomes and detecting amino acids mutations in the spike proteins that can be associated with a phenotypic change. The phenotypic changes are annotated according to the knowledge accumulated on previous variants. Owing to the limited size of the genome, convergent evolution is expected to take place. 
The tool is based on [VOCAL](https://github.com/rki-mf1/vocal): Variant Of Concern ALert and prioritization, which was invented for SARS-CoV-2.


# Documentation VOCAL

<a href="https://rki-mf1.github.io/vocal-doc/"><strong>Explore the docs »</strong></a>


# Getting Started

⚠️**Note**: 🔌 Right now, VOCAL tested on Linux and Mac system only 💻 

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

## Running VOCAL-Flu

```bash
cd vocal-flu
conda activate nextflow
nextflow run main.nf \
     -profile conda,local \
     --fasta 'path/to/input_fasta' \
```

### Running VOCAL-Flu with PSL alignment

```bash
cd vocal-flu
conda activate nextflow
nextflow run main.nf \
     -profile conda,local \
     --fasta 'path/to/input_fasta' \
     --psl 'y'
```

## Parameter list

```
fasta                    The path to the fasta file with the sequences for VOCAL-Flu.
                         [ default: '' ]
subtype                  The Influenza subtype the input sequences are from.
                         The options are H1N1 and H3N2 for Influenza A,
                         Victoria and Yamagata for Influenza B.
                         [ default: 'H1N1' ]
psl                      The parameter enables to generate a PSL file with an alignment.
                         The options are yes/y and no/n.
                         [ default: 'n' ]
```


# Contact

Did you find a bug?🐛 Suggestion/Feedback/Feature request?👨‍💻 please visit [GitHub Issues](https://github.com/rki-mf1/vocal-flu/issues)

For business inquiries or professional support requests 🍺 please contact 
Dr. Hölzer, Martin(<HoelzerM@rki.de>) or Dr. Richard, Hugues (<RichardH@rki.de>)


# Acknowledgments

* Original tool: [VOCAL](https://github.com/rki-mf1/vocal): Variant Of Concern ALert and prioritization 

    * Original Idea: SC2 Evolution Working group at the Robert Koch Institute in Berlin

    * Funding: Supported by the European Centers for Disease Control [grant number ECDC GRANT/2021/008 ECD.12222].