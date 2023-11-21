process REPORT {
    label 'report'

    publishDir "${params.output}/${params.report_dir}", mode: params.publish_dir_mode

    input:
        path variants_with_phenotypes
        path rmd

    output:
        path "fluwarnsystem-alerts-samples-all.csv",                emit: alerts_samples
        path "fluwarnsystem-alerts-clusters-summaries-all.csv",     emit: alerts_clusters
        path "fluwarnsystem-report.html",                           emit: report

    script:
    """
    echo "Step 3: Detect and alert emerging variants"

    prepare_report.R \
        -f ${variants_with_phenotypes} \
        -s "fluwarnsystem-alerts-samples-all.csv" \
        -c "fluwarnsystem-alerts-clusters-summaries-all.csv"

    Rscript --vanilla -e \
        "rmarkdown::render(input = \'${rmd}\', output_file = \'fluwarnsystem-report.html\', params = list(alert_samples = \'fluwarnsystem-alerts-samples-all.csv\', alert_clusters = \'fluwarnsystem-alerts-clusters-summaries-all.csv\'))"
    """

    stub:
    """
    touch fluwarnsystem-alerts-samples-all.csv 
    touch fluwarnsystem-alerts-clusters-summaries-all.csv 
    touch fluwarnsystem-report.html 
    """
}