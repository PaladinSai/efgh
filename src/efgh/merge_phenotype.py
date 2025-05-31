import pandas as pd
import logging

def run_merge(config, ds):
    """
    运行数据处理流程，输入zarr，输出处理后的ds。
    参数:
        config: 配置对象 / config object
        vcz_path: zarr文件路径 / zarr file path
    返回:
        处理后的ds / processed ds
    """
    # 从配置对象获取参数 / Get parameters from config object
    pheno_path = config.input.pheno_path

    logging.info("Merge phenotype started.")
    logging.info(f"Loading phenotype file from {pheno_path} ...")
    try:
        df = pd.read_csv(pheno_path, sep=",", index_col=0)
        df.index.name = "samples"
    except Exception:
        logging.error("Failed to load phenotype file. Please check your phenotype file path and format.")
        raise RuntimeError("Failed to load phenotype file.") from None
    try:
        ds_annotations = df.to_xarray()
        ds = ds.set_index({"samples": "sample_id"})
        ds = ds.merge(ds_annotations, join="left")
        ds = ds.reset_index("samples").reset_coords("samples")
        ds = ds.rename_vars({"samples": "sample_id"})
    except Exception:
        logging.error("Failed to merge phenotype data with genotype data.")
        raise RuntimeError("Failed to merge phenotype data with genotype data.") from None
    logging.info("Merge phenotype finished.")
    return ds
