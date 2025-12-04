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
import anndata as ad
from utility import benchmark_method
import pandas as pd
import scanpy as sc
# ----------------------------------------------------------------------------------------------------------------------

# paths
general_out_path = Path.cwd() / "app" / "out"
general_out_path.mkdir(parents=True, exist_ok=True)


def run_xclone(adata:ad.AnnData = None, tag:str = None, annot_key:str = "cluster.pred"):
    if adata is None:
        adata = xclone.data.tnbc1_rdr()

    # configure
    if tag is None:
        tag = "test_data"

    path_out = general_out_path/f"{tag}__RDR"
    path_out.mkdir(exist_ok=True, parents=True)

    @benchmark_method(str(general_out_path/f"{tag}__RDR"))
    def run_rdr_xclone(adata, tag, annot_key):
        rdrconfig = xclone.XCloneConfig(dataset_name = f"{str(tag)}__RDR_", module = "RDR")
        rdrconfig.set_figure_params(xclone= True, fontsize = 18)
        rdrconfig.outdir = general_out_path / f"{str(tag)}__RDR"
        rdrconfig.cell_anno_key = annot_key  # obs column with annotation key
        rdrconfig.ref_celltype = "N"
        rdrconfig.marker_group_anno_key = annot_key
        rdrconfig.xclone_plot= True
        xclone.model.run_RDR(adata, config_file = rdrconfig)
    run_rdr_xclone(adata, tag, annot_key)
    path_out_data = path_out / "data"
    path_data = [p for p in path_out_data.glob("*.h5ad")][0]
    adata = sc.read_h5ad(path_data)
    df_out = pd.DataFrame(adata.X.todense(), index=list(adata.obs.index), columns=list(adata.var.index)).T.astype(int)
    var_select = adata.var[["chr", "start", "stop"]].rename({"chr":"CHR", "start":"START", "stop":"END"}, axis=1)
    var_select["CHR"] = var_select["CHR"].apply(lambda x: f"chr{x}")
    df_concat = pd.concat([var_select, df_out], axis=1).set_index("CHR")
    df_concat.to_csv(general_out_path/f"{tag}__RDR"/f"{tag}__Xclone__RDR__csv_pred.csv")


if __name__ == "__main__":
    """
    under xclone --> _config.py there are the relevant configuration settings with std parameters!
    HMM:
        trans_t = 1e-6  (general)
        min_iter = 1    (RDR only)
        max_iter = 2    (RDR only)
        HMM_Config seems to be broken, let's leave it be
    """
    run_xclone()