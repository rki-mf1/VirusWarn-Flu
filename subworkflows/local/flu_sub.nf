include { NEXTCLADE }   from '../../modules/local/nextclade'
include { QC }          from '../../modules/local/qc'
include { ANNOTATION }  from '../../modules/local/annotation'
include { REPORT }      from '../../modules/local/report'

workflow FLU_SUB {
    take:
        input_fasta
        moc_table
        roi_table
        fixed_table
        rmd
        qc_rmd
        metadata

    main:
        NEXTCLADE ( input_fasta, params.subtype )

        if (params.qc) {
            QC ( input_fasta, NEXTCLADE.out.nextclade, qc_rmd, params.subtype )
        }

        ANNOTATION ( NEXTCLADE.out.variant_table, moc_table, roi_table )

        REPORT ( 
            NEXTCLADE.out.variant_table,
            ANNOTATION.out.variants_phenotypes, rmd, 
            input_fasta, moc_table, roi_table, 
            metadata, fixed_table, 'report.html' 
        )

    emit:
        report = REPORT.out.report

}