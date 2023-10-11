include { VOCAL }       from '../../modules/local/vocal'
include { PSL }         from '../../modules/local/psl'
include { ANNOTATION }  from '../../modules/local/annotation'
include { REPORT }      from '../../modules/local/report'

workflow VOCAL_FLOW {
    take:
        ref_nt
        input_fasta
        mutation_table

    main:
        if (params.psl == 'yes' || 'y') {
            PSL ( ref_nt, input_fasta )
            VOCAL ( 'PSL', input_fasta, PSL.out.output_psl )
        } else if (params.psl == 'no' || 'n') {
            VOCAL ( '', input_fasta )
        } else {
            exit 1,
            "ERROR: $params.psl is an invalid input for the parameter psl!\n Please choose between yes/y and no/n!\n"
        }

        ANNOTATION ( VOCAL.out.variant_table, mutation_table)

        REPORT ( ANNOTATION.out.variants_with_phenotypes )

    emit:
        report = REPORT.out.report

}