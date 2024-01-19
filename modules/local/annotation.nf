process ANNOTATION {
    label 'annotation'

    publishDir "${params.output}/${params.annot_dir}", mode: params.publish_dir_mode

    input:
        path variant_table
        path mutation_table
        path roi_table

    output:
        path "variants_with_phenotypes.tsv",   emit: variants_phenotypes

    script:
    """
    echo "Step 2: Annotate mutation phenotypes"

    Mutations2Function.py \
        -i ${variant_table} \
        -a ${mutation_table} \
        --roi_table ${roi_table} \
        -o "variants_with_phenotypes.tsv"
    """

    stub:
    """
    touch variants_with_phenotypes.tsv 
    """
}