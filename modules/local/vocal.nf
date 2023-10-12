process VOCAL {
    label 'vocal'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta
        path ref_nt
        path ref_aa

    output:
        path "variant_table.tsv",   emit: variant_table

    script:
    """
    echo "Step 1: Annotate mutations in the proteins"

    vocal.py \
        -i ${input_fasta} \
        --ref_nt ${ref_nt} \
        --ref_aa ${ref_aa} \
        -o "variant_table.tsv"
    """

    stub:
    """
    touch variant_table.tsv
    """
}