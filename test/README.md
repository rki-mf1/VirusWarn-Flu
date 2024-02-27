# Test Datasets

In this folder, there are test datasets for FluWarnSystem which were downloaded from [OpenFlu](https://openflu.vital-it.ch).

The corresponding metadata files couldn't be downloaded due to an error, so alternative metadata files with only Isolate Name and Collection Date were prepared.

## Influenza A H1N1

With the following commands and data set, the mode subtype can be tested for Influenza A H1N1pdm.

### Running test set

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_h1n1.fasta'
```

### Running test set with metadata

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_h1n1.fasta' \
     --metadata 'test/metadata_h1n1.xlsx'
```

### Dataset

The dataset includes 56 sequences and can be downloaded with the following filters:

| Filter                  | Value           |
|:------------------------|----------------:|
| **Virus type**          | A               |
| **H subtype**           | 1               |
| **N subtype**           | 1               |
| **Country**             | Germany         |
| **Collection Date From**| 2019-01-01      |
| **Collection Date To**  | 2019-12-31      |


## Influenza A H3N2

With the following commands and data set, the mode subtype can be tested for Influenza A H3N2.

### Running test set

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_h3n2.fasta' \
     --subtype 'h3n2'
```

### Running test set with metadata

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_h3n2.fasta' \
     --metadata 'test/metadata_h3n2.xlsx' \
     --subtype 'h3n2'
```

### Dataset

The dataset includes 79 sequences and can be downloaded with the following filters:

| Filter                  | Value           |
|:------------------------|----------------:|
| **Virus type**          | A               |
| **H subtype**           | 3               |
| **N subtype**           | 2               |
| **Country**             | Germany         |
| **Collection Date From**| 2019-01-01      |
| **Collection Date To**  | 2019-12-31      |


## Influenza B Victoria

With the following commands and data set, the mode subtype can be tested for Influenza B Victoria.

### Running test set

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_vic.fasta' \
     --subtype 'vic'
```

### Running test set with metadata

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_vic.fasta' \
     --metadata 'test/metadata_vic.xlsx' \
     --subtype 'vic'
```

### Dataset

The dataset includes 2 sequences and can be downloaded with the following filters:

| Filter                  | Value           |
|:------------------------|----------------:|
| **Virus type**          | B               |
| **Lineage**             | Victoria        |
| **Country**             | Germany         |
| **Collection Date From**| 2019-01-01      |
| **Collection Date To**  | 2019-12-31      |


## Combi Set

With the following commands and data set, the mode split can be tested for OpenFlu.

### Running test set

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_combi.fasta' \
     --split 'OpenFlu'
```

### Running test set with metadata

```bash
nextflow run main.nf \
     -profile conda,local \
     --fasta 'test/openflu_combi.fasta' \
     --metadata 'test/metadata_combi.xlsx' \
     --split 'OpenFlu'
```

### Dataset

The combi dataset includes all sequences from the three datasets above.