/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// help message
if (params.help) { exit 0, helpMSG() }

// Parameters sanity checking
Set valid_params = ['cores', 'max_cores', 'memory', 'help',
                    'fasta', 'subtype', 'psl',
                    'output', 'vocal_dir', 'annot_dir', 'report_dir', 'runinfo_dir',
                    'publish_dir_mode', 'conda_cache_dir',
                    'cloudProcess', 'cloud-process']

def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FLUWARNSYSTEM_SUB } from '../subworkflows/local/FluWarnSystem_sub'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLUWARNSYSTEM {

//
// LOAD REFERENCES AND TABLES for choosen Influenza subtype
//
if (params.subtype == 'H1N1' || 'h1n1') {
    log.info"INFO: VOCAL-Flu is running for Influenza A(H1N1)pdm09"
    ref_nt = Channel.fromPath( file("data/A(H1N1)pdm09/A-Brisbane-2-2018_nucleotide.fa", checkIfExists: true) )
    ref_aa = Channel.fromPath( file("data/A(H1N1)pdm09/A-Brisbane-2-2018_protein.fa", checkIfExists: true) )
    mutation_table = Channel.fromPath( file("data/A(H1N1)pdm09/table_h1n1_mutations_annotation.tsv", checkIfExists: true) )
    roi_table = Channel.fromPath( file("data/A(H1N1)pdm09/table_h1n1_roi.csv", checkIfExists: true) )
} else if (params.subtype == 'H3N2' || 'h3n2') {
    log.info"INFO: VOCAL-Flu is running for Influenza A(H3N2)"
    ref_nt = Channel.fromPath( file("data/A(H3N2)/A-Kansas-14-2017_nucleotide.fa", checkIfExists: true) )
    ref_aa = Channel.fromPath( file("data/A(H3N2)/A-Kansas-14-2017_protein.fa", checkIfExists: true) )
    mutation_table = Channel.fromPath( file("data/A(H3N2)/table_h3n2_mutations_annotation.tsv", checkIfExists: true) )
    roi_table = Channel.fromPath( file("data/A(H3N2)/table_h3n2_roi.csv", checkIfExists: true) )
} else if (params.subtype == 'Victoria' || 'victoria') {
    log.info"INFO: VOCAL-Flu is running for Influenza B(Victoria)"
    ref_nt = Channel.fromPath( file("data/B(Victoria)/B-Brisbane-60-2008_nucleotide.fa", checkIfExists: true) )
    ref_aa = Channel.fromPath( file("data/B(Victoria)/B-Brisbane-60-2008_protein.fa", checkIfExists: true) )
    mutation_table = Channel.fromPath( file("data/B(Victoria)/table_victoria_mutations_annotation.tsv", checkIfExists: true) )
    roi_table = Channel.fromPath( file("data/B(Victoria)/table_victoria_roi.csv", checkIfExists: true) )
} else if (params.subtype == 'Yamagata' || 'yamagata') {
    log.info"INFO: VOCAL-Flu is running for Influenza B(Yamagata)"
    ref_nt = Channel.fromPath( file("data/B(Yamagata)/B-Florida-4-2006_nucleotide.fa", checkIfExists: true) )
    ref_aa = Channel.fromPath( file("data/B(Yamagata)/B-Florida-4-2006_protein.fa", checkIfExists: true) )
    mutation_table = Channel.fromPath( file("data/B(Yamagata)/table_yamagata_mutations_annotation.tsv", checkIfExists: true) )
    roi_table = Channel.fromPath( file("data/B(Yamagata)/table_yamagata_roi.csv", checkIfExists: true) )
} else {
    exit 1, 
    "ERROR: $params.subtype is an invalid input for the parameter subtype!\n Please choose between H1N1, H3N2, Victoria and Yamagata!\n"
}

//
// RUN VOCAL
//
input_fasta = Channel.fromPath( file("${params.fasta}", checkIfExists: true) )
rmd = Channel.fromPath( file("bin/report.Rmd", checkIfExists: true) )

FLUWARNSYSTEM_SUB ( ref_nt, ref_aa, input_fasta, mutation_table, roi_table, rmd )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\u001B[31m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

    Workflow: VOCAL-Flu

    ${c_yellow}Usage examples:${c_reset}

    ${c_yellow}Input options:${c_reset}
    ${c_green} --subtype ${c_reset}         Define the subtype of the input sequences to choose the right references and tables.
                                            Options for Influenza A: 'H1N1' and 'H3N2'
                                            Options for Influenza B: 'Victoria' and 'Yamagata'
                                            [ default: $params.subtype ]
    ${c_green} --psl ${c_reset}             Run process with ('y') or without ('n') psl format.
                                            [ default: $params.psl ]

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]

    ${c_yellow}Output options:${c_reset}
    --output                 Name of the result folder [default: $params.output]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: $params.publish_dir_mode]
                                ${c_dim}With 'symlink' results are lost when removing the work directory.${c_reset}

    ${c_yellow}Caching:${c_reset}
    --conda_cache_dir        Location for storing the conda environments [default: $params.conda_cache_dir]
    
    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executors${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executor${c_reset} (choose one):
        local
        slurm
    
    ${c_blue}Engines${c_reset} (choose one):
        conda
        mamba
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/