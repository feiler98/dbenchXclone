# general utility functions for benchmarking x_clone_rdr

# imports
from collections.abc import Callable
from pathlib import Path
import time
import pandas as pd
from configparser import ConfigParser
import logging
import numpy as np
import gtfparse


def get_project_dir(cwd: Path, parent_path_str: str) -> Path:
    """
    Recursive path search until getting parent-path.

    Parameters
    ----------
    cwd: pathlib.Path
        Path that needs to be searched (named after current working directory).
    parent_path_str : str
        String which is contained in the path.

    Returns
    -------
    pathlib.Path
        Query parent-path.
    """

    if cwd.name != parent_path_str:
        parent_path = cwd.parent
        # always return the recursion, otherwise output is None
        return get_project_dir(cwd=parent_path, parent_path_str=parent_path_str)
    else:
        return cwd



class GetConfig:
    """
    Class handling configuration-file [.ini] setup.

    Attributes
    ----------
    cfg_obj : configparser.ConfigParser()
        Loaded config object.
    path_cfg : str
        Absolute-path of ini-file.
    """

    def __init__(self, cfg_obj: ConfigParser = None, path_cfg: str = None):
        self.cfg_obj = cfg_obj
        self.path_cfg = path_cfg

    @classmethod
    def get_config(cls, config_path: str):
        """
        Method for initializing the class GetConfig.

        Parameters
        ----------
            config_path : str
                Path of the config-file [.ini].
        """

        if not Path(config_path).exists():
            raise FileExistsError("Config File could not be found!")
        cfg_obj = ConfigParser()
        cfg_obj.read(str(config_path))
        return cls(cfg_obj=cfg_obj, path_cfg=str(config_path))

    def __check_loaded(self):
        for keys, values in self.__dict__.items():
            if values is None:
                raise ValueError("""
#####################################################
No config-file has been loaded.
Initialize the class by using '.get_config(config_name)'
#####################################################
                                 """)

    def return_section(self, section: str) -> dict:
        """
        Returns the config-section as dictionary.

        Parameters
        ----------
        section : str
            Name of the target-section.

        Returns
        -------
        dict
            Section.
        """

        # check if initialized
        # --------------------
        GetConfig.__check_loaded(self)

        return dict(self.cfg_obj.items(section))

    def return_all(self) -> dict:
        """
        Returns all config-sections as dictionary.

        Returns
        -------
        dict
            All Sections.
        """

        # check if initialized
        # --------------------
        GetConfig.__check_loaded(self)

        return dict(self.cfg_obj.items())

    def get_cfg_sections(self) -> list:
        """
        Returns all config-sections as list.

        Returns
        -------
        list
            List of all section names of the config.
        """

        # check if initialized
        # --------------------
        GetConfig.__check_loaded(self)

        return list(self.cfg_obj.keys())

    def get_repair_config_section(self, section: str, dict_params: dict) -> dict:
        """
        Special function that makes the retrieval of configuration files more robust and 'self-repairing'.
        If you think the user will be a monkey with a hammer, this function is for you!

        Parameters
        ----------
        section : str
            Name of the target-section.
        dict_params : dict
            Dictionary of hardcoded key-value pairs that are necessary for the downstream code to work.

        Returns
        -------
        dict
            Returns section as dictionary + automatically repairs configuration-file with hardcoded key-value pairs
            if needed.
        """

        # check if initialized
        # --------------------
        GetConfig.__check_loaded(self)

        # update config if section is missing
        list_sections_available = GetConfig.get_cfg_sections(self)
        if section not in list_sections_available:
            self.cfg_obj.add_section(section)
            with open(self.path_cfg, 'w') as update_ini:
                self.cfg_obj.write(update_ini)

        section_dict = GetConfig.return_section(self, section)

        # update config if in the section params are missing
        return_section = {}
        for keys, values in dict_params.items():
            if keys not in list(section_dict.keys()):
                self.cfg_obj.set(section, keys, values)
                return_section[keys] = values
                with open(self.path_cfg, 'w') as update_ini:
                    self.cfg_obj.write(update_ini)
            else:
                return_section[keys] = section_dict[keys]
        return return_section


# decorator for benchmarking + saving metadata of func input
# ----------------------------------------------------------
def benchmark_method(path_out: str = str(Path.cwd())) -> Callable:
    """
    Very handy DECORATOR for benchmarking different python functions + saving timing and all parameters in an
    excel-sheet.
    Must be written as '@benchmark_method()' above a function.
    Output excel sheet will have the following naming convention: benchmark__funcName__startTimestamp.xlsx

    Parameters
    ----------
    path_out: str
        Target directory absolute path as string for the excel_workbook with the summary of the benchmarking.
        If no path is supplied or the given path does not exist, the current working directory is set as standard.

    Returns
    -------
    Callable
        Returns the decorated function.

    """
    path_out = Path(str(path_out)).resolve()
    if not path_out.exists():
        print(f"\033[31mGiven path '{str(path_out)}' does not exist!\033[0m")
        path_out = Path.cwd()
    print(f"\033[31mSet output path | {str(path_out)}\033[0m")

    # time the method
    # ---------------------------------------------------------------------------
    def time_wrapper(func):
        func_name = func.__name__
        parm = func.__code__.co_varnames[0:func.__code__.co_argcount]
        def wrapper(*args, **kwargs):
            kwargs.update(dict(zip(parm, args)))
            start_time = time.time()
            target_func = func(**kwargs)
            end_time = time.time()

            # note the time ranges
            # ---------------------------------------------------------------------------
            t_diff = round(abs(start_time - end_time), 2)  # in seconds, rounded to 2 decimal places
            start_strf = time.strftime("%D %T", time.gmtime(start_time))
            end_strf = time.strftime("%D %T", time.gmtime(end_time))

            print(f"""\033[31m
------------------------------------
|     Time results of function     |
------------------------------------
startpoint | {start_strf}
endpoint   | {end_strf}
====================================
t-diff     | {t_diff} s
------------------------------------
                \033[0m""")

            # prepare the DataFrame
            # ---------------------
            dict_df = {"bench_start": start_strf,
                       "bench_end": end_strf,
                       "run_duration [seconds]": t_diff,
                       "func_name": func_name}
            print(kwargs)
            reformed_kwargs = {key:str(value) for key, value in kwargs.items()}
            dict_df.update(reformed_kwargs)
            export_df = pd.DataFrame(dict_df, index=["func_info"]).T

            export_name = f"benchmark__{func_name}__{int(start_time)}.xlsx"

            # create excel workbook
            # ---------------------
            export_df.to_excel(Path(path_out)/export_name)

            return target_func
        return wrapper
    return time_wrapper


# taken from infercnvpy --> shoutout to their clean work
def genomic_position_from_gtf(
    gtf_file: Path,
    adata = None,
    *,
    gtf_gene_id = "gene_name",
    adata_gene_id = None):
    """
    Get genomic gene positions from a GTF file.

    The GTF file needs to match the genome annotation used for your single cell dataset.
    Modified to change always inplace

    .. warning::
        Currently only tested with GENCODE GTFs.

    Parameters
    ----------
    gtf_file
        Path to the GTF file
    adata
        Adds the genomic positions to `adata.var`. If adata is None, returns
        a data frame with the genomic positions instead.
    gtf_gene_id
        Use this GTF column to match it to anndata
    adata_gene_id
        Match this column to the gene ids from the GTF file. Default: use
        adata.var_names.
    """
    try:
        import gtfparse
    except ImportError:
        raise ImportError(
            "genomic_position_from_gtf requires gtfparse as optional dependency. Please install it using `pip install gtfparse`."
        ) from None
    gtf = gtfparse.read_gtf(
        gtf_file, usecols=["seqname", "feature", "start", "end", "gene_id", "gene_name"]
    ).to_pandas()
    gtf = (
        gtf.loc[
            gtf["feature"] == "gene",
            ["seqname", "start", "end", "gene_id", "gene_name"],
        ]
        .drop_duplicates()
        .rename(columns={"seqname": "chromosome"})
    )
    # remove ensembl versions
    gtf["gene_id"] = gtf["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    gene_ids_adata = (adata.var_names if adata_gene_id is None else adata.var[adata_gene_id]).values
    gtf = gtf.loc[gtf[gtf_gene_id].isin(gene_ids_adata), :]

    missing_from_gtf = len(set(gene_ids_adata) - set(gtf[gtf_gene_id].values))
    if missing_from_gtf:
        logging.warning(f"GTF file misses annotation for {missing_from_gtf} genes in adata.")

    duplicated_symbols = np.sum(gtf["gene_name"].duplicated())
    if duplicated_symbols:
        logging.warning(f"Skipped {duplicated_symbols} genes because of duplicate identifiers in GTF file.")
        gtf = gtf.loc[~gtf[gtf_gene_id].duplicated(keep=False), :]

    tmp_var = adata.var.copy()
    orig_index_name = tmp_var.index.name
    TMP_INDEX_NAME = "adata_var_index"
    tmp_var.index.name = TMP_INDEX_NAME
    tmp_var.reset_index(inplace=True)
    var_annotated = tmp_var.merge(
        gtf,
        how="left",
        left_on=TMP_INDEX_NAME if adata_gene_id is None else adata_gene_id,
        right_on=gtf_gene_id,
        validate="one_to_one",
    )
    var_annotated.set_index(TMP_INDEX_NAME, inplace=True)
    var_annotated.index.name = orig_index_name
    var_annotated.rename({'chromosome':"chr", 'end':"stop", 'gene_id':"GeneID", 'gene_name':"GeneName"}, axis=1, inplace=True)
    var_annotated["chr"] = var_annotated["chr"].map(lambda x : str(x).replace("chr", ""))
    dict_chr_arm_bp = {
        "1": 123400000,
        "2": 93900000,
        "3": 90900000,
        "4": 50000000,
        "5": 48800000,
        "6": 59800000,
        "7": 60100000,
        "8": 45200000,
        "9": 43000000,
        "10": 39800000,
        "11": 53400000,
        "12": 35500000,
        "13": 17700000,
        "14": 17200000,
        "15": 19000000,
        "16": 36800000,
        "17": 25100000,
        "18": 18500000,
        "19": 26200000,
        "20": 28100000,
        "21": 12000000,
        "22": 15000000,
        "X": 61000000,
        "Y": 10400000,
    }

    var_annotated["chr_arm"] = var_annotated.index.map(lambda x : f"{var_annotated.loc[x, 'chr']}{'p' if dict_chr_arm_bp[var_annotated.loc[x, 'chr']] > var_annotated.loc[x, 'stop'] else 'q'}" if var_annotated.loc[x, 'chr'] in dict_chr_arm_bp.keys() else np.nan)

    adata.var = var_annotated
    list_var_na_dropped = list(var_annotated.dropna(how="any").index)
    return adata[:, list_var_na_dropped]
