include { VOCAL }       from '../../modules/local/vocal'
include { PSL }         from '../../modules/local/psl'
include { ANNOTATION }  from '../../modules/local/annotation'
include { REPORT }      from '../../modules/local/report'

workflow FLUWARNSYSTEM_SUB {
    take:
        ref_nt
        ref_aa
        input_fasta
        mutation_table
        roi_table
        rmd
        metadata

    main:
        if (params.psl == 'y') {
            PSL ( ref_nt, ref_aa, input_fasta )
            ANNOTATION ( PSL.out.variant_table, mutation_table, roi_table)
        } else if (params.psl == 'n') {
            VOCAL ( input_fasta, ref_nt, ref_aa )
            ANNOTATION ( VOCAL.out.variant_table, mutation_table, roi_table)
        } else {
            exit 1,
            "ERROR: $params.psl is an invalid input for the parameter psl!\n Please choose between yes/y and no/n!\n"
        }

        REPORT ( ANNOTATION.out.variants_with_phenotypes, rmd, input_fasta, mutation_table, roi_table, metadata )

    emit:
        report = REPORT.out.report

}