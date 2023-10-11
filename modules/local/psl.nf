process PSL {
    label 'psl'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        path ref_nt
        path input_fasta

    output:
        path "output.psl",          emit: output_psl

    script:
    """
    echo "Step 1: Generate a PSL file with alingments"

    pblat ${ref_nt} \
        ${input_fasta} \
        -threads=4 \
        "output.psl"
    """

    stub:
    """
    touch output.psl
    touch variant_table.tsv
    """
}