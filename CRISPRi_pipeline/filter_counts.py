import pandas as pd
import numpy as np
import scanpy as sc
import argparse


def get_cell_cycles(adata):
    cell_cycle_genes = [x.strip() for x in open("../CRISPRi_pipeline/data/regev_lab_cell_cycle_genes.txt")]
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


def filter_try_all_counts(counts, feature_calls, guide_list, NC_names, min_genes):
    print("Making variable names unique")
    counts.var_names_make_unique()

    # add feature call data to anndata object
    guide_data = counts.obs.merge(
        feature_calls[["cell_barcode", "feature_call"]].set_index("cell_barcode"),
        left_index=True,
        right_index=True,
    )
    print("Adding target info")
    # add target info for each guide
    guide_data['feature_target'] = guide_data.feature_call.map(guide_list.set_index('feature_call').feature_target)

    # Mask for guides we want to treat as non-working/nontargeting/NC
    # all will be labeled as "No_working_guide" and treated as such in analysis
        # if a 'nontargeting' guide is not passed in the NC file, for example, it will be treated as its own 'working' targeting guide
    NC_cells_mask = guide_data.feature_call.isin(NC_names)
    guide_data["working_features"] = np.nan
    guide_data.loc[:, "working_features"][~NC_cells_mask] = guide_data.feature_target
    guide_data.loc[:, "working_features"][NC_cells_mask] = "No_working_guide"
    # make guide call have universal name standard '_'
    guide_data.feature_call = guide_data.feature_call.str.replace("-", "_")


    # add guide info to counts object
    counts = counts[guide_data.index, :]
    counts.obs = guide_data

    # filter cells
    print("Filtering cells based on number of genes")
    sc.pp.filter_cells(counts, min_genes=min_genes)
    cc_counts = counts.copy()

    print("Getting cell cycle data")
    cc = get_cell_cycles(cc_counts)
    if cc:
        counts.obs[["S_score", "G2M_score", "phase"]] = cc_counts.obs[
            ["S_score", "G2M_score", "phase"]
        ]

    # All "no working guides" are now considered NC
    NC_barcodes = list(
        set(counts.obs[counts.obs.working_features == "No_working_guide"].index)
    )
    # get a counts_NC anndata object with non-working guide cells only, for background features and other QC
    counts_NC = counts[NC_barcodes].copy()
    # filter genes in the NC anddata
    sc.pp.filter_genes(counts_NC, min_cells=counts_NC.shape[0] // 20)

    genes = list(set(counts_NC.var_names))
    # only include genes that have good background coverage
    counts = counts[:, genes].copy()

    return counts, counts_NC


def main():
    # Read in args
    parser = argparse.ArgumentParser(
        description="Read in 10x h5, filter based on features and target guides"
    )
    parser.add_argument(
        "guide_list", type=str,
        help="Tab & line separated file with guide name, target of guide. Columns = [guide] [target]",
    )
    parser.add_argument(
        "NC_list", type=str,
        help="Tab & line separated table with NC names, no header."
    )
    parser.add_argument("counts", type=str,
                        help="h5ad file of counts")
    parser.add_argument(
        "feature_calls", type=str,
        help="Feature call csv of some kind. Must contain cell_barcode and feature_call column. For now, should also contain num_features column.",
    )
    parser.add_argument("min_genes", type=int,
                        help="minimum # genes to keep a cell")
    args = parser.parse_args()

    # read in feature call csv and h5ad
    feature_calls = pd.read_csv(args.feature_calls)
    counts = sc.read_10x_h5(args.counts)

    # for now, filter to cells with only one feature call
    feature_calls = feature_calls.query("num_features == 1")
    # filter to cells present in h5ad
    feature_calls = feature_calls[feature_calls.cell_barcode.isin(counts.obs_names)]

    # Read in the names of all guides and their targets
    guide_list = pd.read_table(args.guide_list, header=0, index_col=0)
    # set universal column names
    assert len(guide_list.columns) == 2
    guide_list.columns = ["feature_call", "feature_target"]

    # Read in names/categories for nontargeting/negative control guides
    NC_df = pd.read_table(args.NC_list, header=None, index_col=0)
    NC_names = NC_df.values.squeeze()

    final_counts, NC_counts = filter_try_all_counts(counts, feature_calls, guide_list, NC_names, args.min_genes)
    final_counts.write_h5ad('filtered_counts.h5ad')
    NC_counts.write_h5ad('filtered_background_counts.h5ad')
    print("Done.")




if __name__ == "__main__":
    main()
