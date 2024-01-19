include { SPLIT }                                                       from '../../modules/local/split'
include { VOCAL as VOCAL_H1N1; VOCAL as VOCAL_H3N2 }                    from '../../modules/local/vocal'
include { ANNOTATION as ANNOTATION_H1N1; ANNOTATION as ANNOTATION_H3N2} from '../../modules/local/annotation'
include { REPORT as REPORT_H1N1; REPORT as REPORT_H3N2 }                from '../../modules/local/report'

workflow FLUWARNSYSTEM_SPLIT {
    take:
        input_fasta
        ref_h1n1
        control_h1n1
        moc_table_h1n1
        roi_table_h1n1
        ref_h3n2
        control_h3n2
        moc_table_h3n2
        roi_table_h3n2
        rmd
        metadata

    main:
        SPLIT ( input_fasta )

        if (SPLIT.out.h1n1_ha) {
            VOCAL_H1N1 ( SPLIT.out.h1n1_ha, ref_h1n1, control_h1n1 )

            ANNOTATION_H1N1 ( VOCAL_H1N1.out.variant_table, moc_table_h1n1, roi_table_h1n1 )

            REPORT_H1N1 ( 
                ANNOTATION_H1N1.out.variants_phenotypes, rmd, 
                SPLIT.out.h1n1_ha, moc_table_h1n1, roi_table_h1n1, 
                metadata, 'report-h1n1.html'
            )
        }

        if (SPLIT.out.h3n2_ha) {
            VOCAL_H3N2 ( SPLIT.out.h3n2_ha, ref_h3n2, control_h3n2 )

            ANNOTATION_H3N2 ( VOCAL_H3N2.out.variant_table, moc_table_h3n2, roi_table_h3n2 )

            REPORT_H3N2 ( 
                ANNOTATION_H3N2.out.variants_phenotypes, rmd, 
                SPLIT.out.h3n2_ha, moc_table_h3n2, roi_table_h3n2, 
                metadata, 'report-h3n2.html'
            )
        }

    emit:
        report_h1n1 = REPORT_H1N1.out.report
        report_h3n2 = REPORT_H3N2.out.report

}