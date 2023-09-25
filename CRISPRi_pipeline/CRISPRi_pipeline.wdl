version 1.0

workflow CRISPRi_pipeline {
    input {
        String directory_10x
        File guide_list
        File NC_list
        Int genes_per_cell_threshold=200
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/crispri_pipeline:0.4'
    }

    call read_counts_and_filter {
        input:
            directory_10x = directory_10x,
            guide_list = guide_list,
            NC_list = NC_list,
            genes_per_cell_threshold = genes_per_cell_threshold,
            docker_image = docker_image
    }

    output {
        File counts_h5ad=read_counts_and_filter.counts_h5ad
        File background_counts_h5ad=read_counts_and_filter.background_counts_h5ad
    }
}

task read_counts_and_filter {
    input {
        String directory_10x
        File guide_list
        File NC_list
        Int genes_per_cell_threshold=200
        String docker_image
    }

    String crispr_analysis_path = directory_10x + "/crispr_analysis.tar.gz"
    String filtered_feature_bc_matrix_path = directory_10x + "/filtered_feature_bc_matrix.h5"

    command {
        set -ex
        (git clone https://github.com/broadinstitute/CRISPRi-pipeline.git /app ; cd /app)
        gsutil cp ${crispr_analysis_path} crispr_analysis.tar.gz
        gsutil cp ${filtered_feature_bc_matrix_path} filtered_feature_bc_matrix.h5
        tar -xf crispr_analysis.tar.gz protospacer_calls_per_cell.csv
        python3 -u /app/CRISPRi_pipeline/filter_counts.py ${guide_list} ${NC_list} filtered_feature_bc_matrix.h5 protospacer_calls_per_cell.csv ${genes_per_cell_threshold}
    }
    output {
        File counts_h5ad = "filtered_counts.h5ad"
        File background_counts_h5ad = "filtered_background_counts.h5ad"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}