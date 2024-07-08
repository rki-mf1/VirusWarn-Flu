process ANNOTATION {
    label 'annotation'

    publishDir "${params.output}/${params.annot_dir}", mode: params.publish_dir_mode

    input:
        path variant_table
        path moc_table
        path roi_table
        path fixed_table

    output:
        path "variants_with_phenotypes.tsv",                    emit: variants_phenotypes
        path "variants_with_phenotypes_without_fixed.tsv",      emit: variants_phenotypes_filtered

    script:
    """
    echo "Step 2: Annotate mutation phenotypes"

    Mutations2Function.py \
        -i ${variant_table} \
        -a ${moc_table} \
        --roi_table ${roi_table} \
        -o "variants_with_phenotypes.tsv"

    echo "Filtering substitutions which are fixed in the population..."

    remove_fixed.py \
        -i "variants_with_phenotypes.tsv" \
        -f ${fixed_table} \
        -s ${params.season} \
        -o "variants_with_phenotypes_without_fixed.tsv"

    """

    stub:
    """
    touch variants_with_phenotypes.tsv 
    touch variants_with_phenotypes_without_fixed.tsv
    """
}