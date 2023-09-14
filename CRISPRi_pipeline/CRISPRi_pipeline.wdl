version 1.0

workflow CRISPRi_pipeline {
    input {
        String directory_10x
        File guide_list
        File NC_list
    }

    call read_counts_and_filter {
        input:
            directory_10x=directory_10x,
            guide_list=guide_list
    }
}

task read_counts_and_filter {
    input {
        String directory_10x
        File guide_list
        File NC_list
    }

    String crispr_analysis_path = directory_10x + "/crispr_analysis.tar.gz"
    String filtered_feature_bc_matrix_path = directory_10x + "/filtered_feature_bc_matrix.h5"
    command {
        (git clone https://github.com/broadinstitute/CRISPRi-pipeline.git /app ; cd /app)
        gsutil cp ${crispr_analysis_path} crispr_analysis.tar.gz
        gsutil cp ${filtered_feature_bc_matrix_path} filtered_feature_bc_matrix.h5
        tar -xf crispr_analysis.tar.gz protospacer_calls_per_cell.csv
        python3 -u /app/CRISPRi_pipeline/filter_counts.py ${guide_list} ${NC_list} filtered_feature_bc_matrix.h5 protospacer_calls_per_cell.csv
    }
}