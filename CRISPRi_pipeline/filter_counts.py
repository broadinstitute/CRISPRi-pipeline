import pandas as pd
import numpy as np
import scanpy as sc
import argparse


def get_cell_cycles(adata):
    cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    if not set(cell_cycle_genes) & set(s_genes) or not set(cell_cycle_genes) & set(
        g2m_genes
    ):
        return False
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    return True


def filter_try_all_counts(day_counts, day_CRISPR_counts, work_list, min_genes):
    counts = day_counts.copy()
    CRISPR_counts = day_CRISPR_counts.copy()

    counts.var_names_make_unique()
    # add feature call data to anndata object
    guide_data = counts.obs.merge(
        CRISPR_counts[["cell_barcode", "feature_call"]].set_index("cell_barcode"),
        left_index=True,
        right_index=True,
    ).copy()
    # make guide call have universal name standard '_'
    guide_data.feature_call = guide_data.feature_call.str.replace("-", "_")

    ## TODO: change to be receptive to a guide list input
    ## This is all unfinished work in progress
    working_cells = guide_data.feature_call.isin(work_list.split(","))
    guide_data["working_features"] = np.nan
    guide_data.working_features[working_cells] = guide_data.feature_call[
        working_cells
    ].str.split("_", expand=True)[0]
    guide_data.working_features[~working_cells] = "No_working_guide"
    counts = counts[guide_data.index, :]
    counts.obs = guide_data

    sc.pp.filter_cells(counts, min_genes=min_genes)
    print(counts)
    cc_counts = counts.copy()

    cc = get_cell_cycles(cc_counts)
    if cc:
        counts.obs[["S_score", "G2M_score", "phase"]] = cc_counts.obs[
            ["S_score", "G2M_score", "phase"]
        ]

    if CRISPR_counts.feature_call.str.contains("NC").any():
        NC_barcodes = list(
            set(
                CRISPR_counts.cell_barcode[
                    CRISPR_counts.feature_call.str.contains("NC")
                ]
            )
            & set(counts.obs_names)
        )
        counts.obs.replace(to_replace=r"NC.*", value="NC", regex=True, inplace=True)
        counts_NC = counts[NC_barcodes].copy()
        sc.pp.filter_genes(counts_NC, min_cells=counts_NC.shape[0] // 20)
    else:
        NC_barcodes = list(
            set(counts.obs[counts.obs.working_features == "No_working_guide"].index)
        )
        counts_NC = counts[NC_barcodes].copy()
        sc.pp.filter_genes(counts_NC, min_cells=counts_NC.shape[0] // 20)

    print(counts_NC)
    genes = list(set(counts_NC.var_names))
    counts = counts[:, genes].copy()

    return counts, counts_NC


def main():
    # Read in args
    parser = argparse.ArgumentParser(
        description="Read in 10x h5, filter based on features and target guides"
    )
    parser.add_argument(
        "guide_list",
        type=str,
        help="Tab separated file with guide name, target of guide. Columns = [guide] [target]",
    )
    parser.add_argument(
        "NC_list", type=str, help="Tab separated table with NC names, no header."
    )
    parser.add_argument("counts", type=str, help="h5ad file of counts")
    parser.add_argument(
        "feature_calls",
        type=str,
        help="Feature call csv of some kind. Must contain cell_barcode and feature_call column. For now, should also contain num_features column.",
    )
    args = parser.parse_args()

    # read in feature call csv and h5ad
    feature_calls = pd.read_csv(args.feature_call)
    counts = sc.read_10x_h5(args.counts)

    # for now, filter to cells with only one feature call
    feature_calls = feature_calls.query("num_features == 1")
    # filter to cells present in h5ad
    feature_calls = feature_calls[feature_calls.cell_barcode.isin(counts.obs_names)]

    # Read in the names of all guides and their targets
    guide_list = pd.read_table(args.guide_list, header=True)
    # Read in names/categories for nontargeting/negative control guides
    NC_df = pd.read_table(args.NC_list, header=None)
    NC_names = NC_df.values.squeeze()

    working_guides = guide_list[~guide_list.target.isin(NC_names)].guide.values


if __name__ == "__main__":
    main()
