include { SPLIT }                                                                                       from '../../modules/local/split'
include { QC as QC_H1N1; QC as QC_H3N2; QC as QC_VIC }                                                  from '../../modules/local/qc'
include { NEXTCLADE as NEXTCLADE_H1N1; NEXTCLADE as NEXTCLADE_H3N2; NEXTCLADE as NEXTCLADE_VIC }        from '../../modules/local/nextclade'
include { ANNOTATION as ANNOTATION_H1N1; ANNOTATION as ANNOTATION_H3N2; ANNOTATION as ANNOTATION_VIC }  from '../../modules/local/annotation'
include { REPORT as REPORT_H1N1; REPORT as REPORT_H3N2; REPORT as REPORT_VIC }                          from '../../modules/local/report'

workflow FLU_SPLIT {
    take:
        input_fasta
        moc_table_h1n1
        roi_table_h1n1
        fixed_table_h1n1
        moc_table_h3n2
        roi_table_h3n2
        fixed_table_h3n2
        moc_table_vic
        roi_table_vic
        fixed_table_vic
        rmd
        qc_rmd
        metadata

    main:
        SPLIT ( input_fasta )

        if (SPLIT.out.h1n1_ha) {
            NEXTCLADE_H1N1 ( SPLIT.out.h1n1_ha, 'h1n1' )

            if (params.qc) {
                QC_H1N1 ( input_fasta, NEXTCLADE_H1N1.out.nextclade, qc_rmd, 'h1n1' )
            }

            ANNOTATION_H1N1 ( NEXTCLADE_H1N1.out.variant_table, moc_table_h1n1, roi_table_h1n1 )

            REPORT_H1N1 ( 
                NEXTCLADE_H1N1.out.variant_table,
                ANNOTATION_H1N1.out.variants_phenotypes, rmd, 
                SPLIT.out.h1n1_ha, moc_table_h1n1, roi_table_h1n1, 
                metadata, fixed_table_h1n1, 'report-h1n1.html'
            )
        }

        if (SPLIT.out.h3n2_ha) {
            NEXTCLADE_H3N2 ( SPLIT.out.h3n2_ha, 'h3n2' )

            if (params.qc) {
                QC_H3N2 ( input_fasta, NEXTCLADE_H3N2.out.nextclade, qc_rmd, 'h3n2' )
            }

            ANNOTATION_H3N2 ( NEXTCLADE_H3N2.out.variant_table, moc_table_h3n2, roi_table_h3n2 )

            REPORT_H3N2 ( 
                NEXTCLADE_H3N2.out.variant_table,
                ANNOTATION_H3N2.out.variants_phenotypes, rmd, 
                SPLIT.out.h3n2_ha, moc_table_h3n2, roi_table_h3n2, 
                metadata, fixed_table_h3n2, 'report-h3n2.html'
            )
        }

        if (SPLIT.out.vic_ha) {
            NEXTCLADE_VIC ( SPLIT.out.vic_ha, 'vic' )

            if (params.qc) {
                QC_VIC ( input_fasta, NEXTCLADE_VIC.out.nextclade, qc_rmd, 'vic' )
            }

            ANNOTATION_VIC ( NEXTCLADE_VIC.out.variant_table, moc_table_vic, roi_table_vic )

            REPORT_VIC ( 
                NEXTCLADE_VIC.out.variant_table,
                ANNOTATION_VIC.out.variants_phenotypes, rmd, 
                SPLIT.out.vic_ha, moc_table_vic, roi_table_vic, 
                metadata, fixed_table_vic, 'report-vic.html'
            )
        }

    emit:
        report_h1n1 = REPORT_H1N1.out.report
        report_h3n2 = REPORT_H3N2.out.report
        report_vic = REPORT_VIC.out.report

}