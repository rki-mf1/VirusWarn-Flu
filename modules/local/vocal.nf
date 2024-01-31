process VOCAL {
    label 'vocal'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta
        path ref
        path control

    output:
        path "variant_table.tsv",   emit: variant_table

    script:
    """
    echo "Step 1: Annotate mutations in the proteins"

    vocal.py \
        -i ${input_fasta} \
        -r ${ref} \
        --control ${control} \
        --complete ${params.complete} \
        -n ${params.n} \
        -o "variant_table.tsv"
    """

    stub:
    """
    touch variant_table.tsv
    """
}