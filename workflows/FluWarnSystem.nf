/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// help message
if (params.help) { exit 0, helpMSG() }

// Parameters sanity checking
Set valid_params = ['cores', 'max_cores', 'memory', 'help',
                    'fasta', 'refh1n1', 'refh3n2', 'refvic',
                    'metadata', 'subtype', 
                    'psl', 'split', 'complete', 'n',
                    'output', 'split_dir', 'nextclade_dir', 
                    'annot_dir', 'report_dir', 'runinfo_dir',
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

include { FLUWARNSYSTEM_SUB }   from '../subworkflows/local/FluWarnSystem_sub'
include { FLUWARNSYSTEM_SPLIT } from '../subworkflows/local/FluWarnSystem_split'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLUWARNSYSTEM {

    input_fasta = Channel.fromPath( file("${params.fasta}", checkIfExists: true) )

    if (params.refh1n1 != '') {
        ref_h1n1 = Channel.fromPath( file("${params.refh1n1}", checkIfExists: true) )
    } else {
        ref_h1n1 = params.refh1n1
    }
    control_h1n1 = Channel.fromPath( file("data/A(H1N1)pdm09/HA_moc_roi_seq.fasta", checkIfExists: true) )
    moc_table_h1n1 = Channel.fromPath( file("data/A(H1N1)pdm09/table_h1n1_moc_HA.tsv", checkIfExists: true) )
    roi_table_h1n1 = Channel.fromPath( file("data/A(H1N1)pdm09/table_h1n1_roi.csv", checkIfExists: true) )     

    if (params.refh3n2 != '') {
        ref_h3n2 = Channel.fromPath( file("${params.refh3n2}", checkIfExists: true) )
    } else {
        ref_h3n2 = params.refh3n2
    }
    control_h3n2 = Channel.fromPath( file("data/A(H3N2)/HA_moc_roi_seq.fasta", checkIfExists: true) )
    moc_table_h3n2 = Channel.fromPath( file("data/A(H3N2)/table_h3n2_moc_HA.tsv", checkIfExists: true) )
    roi_table_h3n2 = Channel.fromPath( file("data/A(H3N2)/table_h3n2_roi.csv", checkIfExists: true) )

    if (params.refvic != '') {
        ref_vic = Channel.fromPath( file("${params.refvic}", checkIfExists: true) )
    } else {
        ref_vic = params.refvic
    }
    control_vic = Channel.fromPath( file("data/B(Victoria)/HA_moc_roi_seq.fasta", checkIfExists: true) )
    moc_table_vic = Channel.fromPath( file("data/B(Victoria)/table_vic_moc_HA.tsv", checkIfExists: true) )
    roi_table_vic = Channel.fromPath( file("data/B(Victoria)/table_vic_roi.csv", checkIfExists: true) )

    rmd = Channel.fromPath( file("bin/report.Rmd", checkIfExists: true) )

    if (params.metadata != '') {
        metadata = Channel.fromPath( file("${params.metadata}", checkIfExists: true) )
    } else {
        metadata = params.metadata
    }

    if (params.split == ''){

        if (params.subtype == 'h1n1') {
            log.info"INFO: FluWarnSystem is running for Influenza A(H1N1)pdm09"
            ref_nt = ref_h1n1
            control = control_h1n1
            mutation_table = moc_table_h1n1
            roi_table = roi_table_h1n1
        } else if (params.subtype == 'h3n2') {
            log.info"INFO: FluWarnSystem is running for Influenza A(H3N2)"
            ref_nt = ref_h3n2
            control = control_h3n2
            mutation_table = moc_table_h3n2
            roi_table = roi_table_h3n2
        } else if (params.subtype == 'vic') {
            log.info"INFO: FluWarnSystem is running for Influenza B(Victoria)"
            ref_nt = ref_vic
            control = control_vic
            mutation_table = moc_table_vic
            roi_table = roi_table_vic
        } else {
            exit 1, 
            "ERROR: $params.subtype is an invalid input for the parameter subtype!\n Please choose between H1N1, H3N2, Victoria and Yamagata!\n"
        }

        FLUWARNSYSTEM_SUB ( ref_nt, input_fasta, control, mutation_table, roi_table, rmd, metadata )

    } else if (params.split == 'FluPipe' || 'flupipe' || 'GISAID' || 'gisaid' || 'OpenFlu' || 'openflu') {

        if (params.psl == true) {
            exit 1, 
            "ERROR: $params.psl is an invalid input for the parameter split!\n Please choose between FluPipe, GISAID and OpenFlu!\n"
        }

        FLUWARNSYSTEM_SPLIT ( 
            input_fasta, 
            ref_h1n1, control_h1n1, moc_table_h1n1, roi_table_h1n1,
            ref_h3n2, control_h3n2, moc_table_h3n2, roi_table_h3n2, 
            ref_vic, control_vic, moc_table_vic, roi_table_vic, 
            rmd, metadata 
        )

    } else {
        exit 1, 
            "ERROR: $params.split is an invalid input for the parameter split!\n Please choose between FluPipe, GISAID and OpenFlu!\n"
    }

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

    Workflow: FluWarnSystem

    ${c_yellow}Usage examples:${c_reset}

    ${c_yellow}Input options:${c_reset}
    ${c_green} --fasta ${c_reset}           REQUIRED! Path to the input fasta file.
                                            [ default: $params.fasta ]
    ${c_green} --refh1n1 ${c_reset}         Path to the reference sequence for H1N1.
                                            Otherwise A/Wisconsin/588/2019 (MW626062.1) is used.
                                            [ default: $params.refh1n1 ]
    ${c_green} --refh3n2 ${c_reset}         Path to the reference sequence for H3N2.
                                            Otherwise A/Darwin/6/2021 (EPI_ISL_1563628) is used.
                                            [ default: $params.refh3n2 ]
    ${c_green} --refvic ${c_reset}          Path to the reference sequence for Victoria.
                                            Otherwise B/Brisbane/60/2008 (KX058884.1) is used.
                                            [ default: $params.refvic ]                                        
    ${c_green} --metadata ${c_reset}        The path to a metadata file in GISAID form for the sequences with 
                                            collection dates.
                                            Required to generate a heatmap in the report.
                                            [ default: $params.metadata ]
    ${c_green} --subtype ${c_reset}         If the input fasta file only contains sequences of one subtype, 
                                            define the subtype to choose the right references and tables.
                                            Options for Influenza A: 'H1N1' and 'H3N2'
                                            Options for Influenza B: 'Victoria'
                                            [ default: $params.subtype ]
    ${c_green} --psl ${c_reset}             Run process with ('y') or without ('n') psl format.
                                            [ default: $params.psl ]
    ${c_green} --split ${c_reset}           If the input fasta file contains sequences of more than one subtype, 
                                            enable the split parameter to write them into one file per subtype and 
                                            ensure the use of the right references and tables.
                                            Options: 'FluPipe', 'GISAID' and 'OpenFlu'
                                            [ default: $params.split ]
    ${c_green} --complete ${c_reset}        FluWarnSystem only considers sequences within a defined range of length
                                            and writes the rest into incomplete_seq.fasta if set to 'y'.
                                            If set to 'n', all sequences are considered.
                                            [ default: $params.complete ]
    ${c_green} --n ${c_reset}               Number of nucleotides a sequence can differ from the length of 
                                            the reference sequence to be considered as complete.
                                            [ default: $params.n ]

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