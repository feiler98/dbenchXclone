# to do list:
"""docarize a function of xclone where it outputs h5ad or csv of the result of the analysis
    must be in python 3.9
    two differend human genomes available
    hg38 must be implemented in var
    obs only needs sample names and externally derived cell classification (N or T for normal and tumor samples)

    Documentation:
    ==============
    https://xclone-cnv.readthedocs.io/en/latest/preprocessing.html
"""

# imports
# ----------------------------------------------------------------------------------------------------------------------
import xclone
from pathlib import Path
from utility import benchmark_method, genomic_position_from_gtf
import pandas as pd
import scanpy as sc
from scipy import sparse
import itertools
# ----------------------------------------------------------------------------------------------------------------------


def grid_by_dict(pars_dict: dict) -> list:
    keys=pars_dict.keys()
    combinations = itertools.product(*pars_dict.values())
    list_of_kwargs = [dict(zip(keys, cc)) for cc in combinations]
    return list_of_kwargs


def val_build_project() -> (Path, Path):
    cwd_path = Path.cwd()
    print(f"Current working directory of running script {Path(__file__).name}: {cwd_path}")
    path_out = cwd_path / "app" / "out"
    path_in = cwd_path / "data_input"

    if not path_in.exists():
        raise ValueError(f"Data dir '{str(path_in)}' does not exist!")

    if not path_out.exists():
        path_out.mkdir(parents=True, exist_ok=True)
        print(f"Data out-dir '{str(path_out)}' has been created...")

    return path_in, path_out


def get_hg_38_desc_paths(target_path: Path) -> dict:
    """
    These fetched .txt files correlate to .csv RCM files --> describe normal cells within the datasets.
    """
    return {p.stem: p for p in target_path.rglob("*__hg_38.txt")}


def csvs_to_adatas(target_path: Path) -> dict:
    """
    Generates a dictionary with adata and their respective reference catalogue of normal cells (cell_names).
    """
    dict_hg38_desc = get_hg_38_desc_paths(target_path)
    dict_accepted_files = {}
    for k, path_txt in dict_hg38_desc.items():
        path_rcm = Path(target_path) / f"{k}__RCM.csv"
        if path_rcm.exists():
            adata = sc.read_csv(path_rcm).T
            adata.obs["cell_names"] = adata.obs.index
            with open(path_txt, "r") as f:
                list_norm_cells = list(map(lambda x: x.replace("\n", ""), f.readlines()))
            adata.obs["ref_cells"] = ["N" if idx in list_norm_cells else "unknown" for idx in adata.obs.index]
            adata.X = sparse.csr_matrix(adata.X)
            adata.layers["raw_expr"] = adata.X
            # add most of the information & chr_arm as it is required
            adata = genomic_position_from_gtf(Path("gencode.v38.annotation.gtf"), adata)
            dict_accepted_files[k] = adata
    return dict_accepted_files


def run_xclone(path_target: Path, path_out_data: Path, kwargs: dict = {}):

    default_kwargs = {
        "smart_transform": False,
        "filter_ref_ave": 0.5,
        "min_gene_keep_num": 3000,
        "multi_refcelltype": False,
        "get_marker_genes": True,
        "top_n_marker": 15,
        "remove_marker": True,
        "fit_GLM_libratio": False,
        "select_normal_chr_num": 4,
        "WMA_window_size": 40,
        "ab_k_neighbors": 10,
        "ab_pseudo_count": 1e-6,
        "denoise_sd_amplifier": 1.5,
        "low_rank": False,
        "low_rank_n_components": 10
    }

    # adding keys with default values if missing
    for default_key, val in default_kwargs.items():
        if default_key not in kwargs.keys():
            kwargs[default_key] = val


    dict_files = csvs_to_adatas(path_target)
    for tag_dataset, adata in dict_files.items():
        str_kwargs = ";".join([f"{list(x)[0]},{y}" for x, y in kwargs.items()])
        file_name = f"{tag_dataset}__{str_kwargs}__Xclone_RDR"

        path_out = path_out_data / file_name
        path_out.mkdir(exist_ok=True, parents=True)
        print(adata.var)
        @benchmark_method(str(path_out))
        def run_rdr_xclone(adata, file_name, kwargs):
            rdrconfig = xclone.XCloneConfig(dataset_name = file_name, module = "RDR")

            # base settings based on the tutorial
            rdrconfig.set_figure_params(xclone= True, fontsize = 18)
            rdrconfig.outdir = path_out_data / file_name
            rdrconfig.cell_anno_key = "ref_cells"  # obs column with annotation key
            rdrconfig.ref_celltype = "N"
            rdrconfig.marker_group_anno_key = "ref_cells"
            rdrconfig.xclone_plot= True

            # advanced settings
            # --------------------------------------------------------------------------------------------------------------
            rdrconfig.smart_transform = kwargs["smart_transform"]
            rdrconfig.filter_ref_ave = kwargs["filter_ref_ave"]
            rdrconfig.min_gene_keep_num = kwargs["min_gene_keep_num"]
            rdrconfig.multi_refcelltype = kwargs["multi_refcelltype"]
            rdrconfig.get_marker_genes = kwargs["get_marker_genes"]
            rdrconfig.top_n_marker = kwargs["top_n_marker"]
            rdrconfig.remove_marker = kwargs["remove_marker"]
            rdrconfig.fit_GLM_libratio = kwargs["fit_GLM_libratio"]
            rdrconfig.select_normal_chr_num = kwargs["select_normal_chr_num"]
            ## smoothing
            rdrconfig.WMA_window_size = kwargs["WMA_window_size"]
            # adaptive baseline
            rdrconfig.ab_k_neighbors = kwargs["ab_k_neighbors"]
            rdrconfig.ab_pseudo_count = kwargs["ab_pseudo_count"]
            # denoise
            rdrconfig.denoise_sd_amplifier = kwargs["denoise_sd_amplifier"]
            # low rank
            rdrconfig.low_rank = kwargs["low_rank"]
            rdrconfig.low_rank_n_components = kwargs["low_rank_n_components"]
            # --------------------------------------------------------------------------------------------------------------

            # apply config to the RDR model
            # -----------------------------
            xclone.model.run_RDR(adata, config_file = rdrconfig)

        # run model
        run_rdr_xclone(adata, file_name, kwargs)

        # fetch and transform h5ad to csv standard format
        path_out_data = path_out / "data"
        path_data = [p for p in path_out_data.glob("*.h5ad")][0]
        adata = sc.read_h5ad(path_data)
        df_out = pd.DataFrame(adata.layers["WMA_smoothed_log_ratio_ab"], index=list(adata.obs.index), columns=list(adata.var.index)).T.astype(int)
        var_select = adata.var[["chr", "start", "stop"]].rename({"chr":"CHR", "start":"START", "stop":"END"}, axis=1)
        var_select["CHR"] = var_select["CHR"].apply(lambda x: f"chr{x}")
        df_concat = pd.concat([var_select, df_out], axis=1).set_index("CHR")
        df_concat.to_csv(path_out / f"{file_name}__pred.csv")


if __name__ == "__main__":
    """
    under xclone --> _config.py there are the relevant configuration settings with std parameters!
    HMM:
        trans_t = 1e-6  (general)
        min_iter = 1    (RDR only)
        max_iter = 2    (RDR only)
        HMM_Config seems to be broken, let's leave it be
    """

    kwargs_gridsearch = {
        "smart_transform": [True, False],
        "filter_ref_ave": [0.2, 0.5, 1.0],
        "min_gene_keep_num": [1000, 3000, 10000, 20000],
        "multi_refcelltype": [True, False],
        "get_marker_genes": [True, False],
        "top_n_marker": [5, 15, 50],
        "remove_marker": [True, False],
        "fit_GLM_libratio": [True, False],
        "select_normal_chr_num": [1, 4, 10],
        "WMA_window_size": [20, 40, 100],
        "ab_k_neighbors": [3, 10, 30],
        "ab_pseudo_count": [1e-8, 1e-6, 1e-4],
        "denoise_sd_amplifier": [1.0, 1.5, 2, 4],
        "low_rank": [True, False],
        "low_rank_n_components": [5, 10, 20]
    }

    path_in, path_out = val_build_project()
    run_xclone(path_in, path_out)  # standard params
    list_kwargs = grid_by_dict(kwargs_gridsearch)
    for kwarg_opt in list_kwargs:
        print(f"InferCNVpy running with hyperparameters: {kwarg_opt}")
        run_xclone(path_in, path_out, kwargs=kwarg_opt)