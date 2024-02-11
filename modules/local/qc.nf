process QC {
    label 'qc'

    publishDir "${params.output}/${params.qc_dir}", mode: params.publish_dir_mode

    input:
        path input_fasta
        path nextclade
        path qc_rmd
        val subtype

    output:
        path "qc_table.tsv",   emit: qc_table
        path "qc-report.html", emit: qc_report

    script:
    """
    echo "Quality Control Step: Summarize QC results from Nextclade"

    qc.py \
        -i ${nextclade} \
        -o "qc_table.tsv"

    Rscript --vanilla -e \
        "rmarkdown::render(input = \'${qc_rmd}\', \\
        output_file = \'qc-report.html\', \\
        params = list(
            fasta = \'${input_fasta}\', \\
            qc_table = \'qc_table.tsv\', \\
            subtype = \'${subtype}\')
        )"
    """

    stub:
    """
    touch qc_table.tsv
    touch qc-report.html
    """
}