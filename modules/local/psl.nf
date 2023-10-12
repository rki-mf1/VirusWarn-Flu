process PSL {
    label 'psl'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        path ref_nt
        path ref_aa
        path input_fasta

    output:
        path "variant_table.tsv",   emit: variant_table

    script:
    """
    echo "Step 1: Generate a PSL file with alingments"

    pblat ${ref_nt} \
        ${input_fasta} \
        -threads=4 \
        "output.psl"

    vocal.py \
        -i ${input_fasta} \
        --ref_nt ${ref_nt} \
        --ref_aa ${ref_aa} \
        --PSL "output_psl" \
        -o "variant_table.tsv"
    """

    stub:
    """
    touch output.psl
    touch variant_table.tsv
    """
}