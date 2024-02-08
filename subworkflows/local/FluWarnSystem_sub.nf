include { NEXTCLADE }   from '../../modules/local/nextclade'
include { ANNOTATION }  from '../../modules/local/annotation'
include { REPORT }      from '../../modules/local/report'

workflow FLUWARNSYSTEM_SUB {
    take:
        input_fasta
        moc_table
        roi_table
        rmd
        metadata

    main:
        NEXTCLADE ( input_fasta, params.subtype )

        ANNOTATION ( NEXTCLADE.out.variant_table, moc_table, roi_table )

        REPORT ( 
            NEXTCLADE.out.variant_table,
            ANNOTATION.out.variants_phenotypes, rmd, 
            input_fasta, moc_table, roi_table, 
            metadata, params.subtype, 'report.html' 
        )

    emit:
        report = REPORT.out.report

}