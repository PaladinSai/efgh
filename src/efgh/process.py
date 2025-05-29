import pandas as pd
import logging

def run_process(config, ds):
    """
    运行数据处理流程，输入zarr，输出处理后的ds。
    Run data processing workflow, input zarr, output processed ds.
    参数:
        config: 配置对象 / config object
        vcz_path: zarr文件路径 / zarr file path
    返回:
        处理后的ds / processed ds
    """
    # 从配置对象获取参数 / Get parameters from config object
    pheno_path = config.input.pheno_path

    logging.info("Data preprocessing started.")
    ds["variant_contig_name"] = ds.contig_id[ds.variant_contig]

    # 列出所有变量（包括坐标变量）
    # List all variables (including coordinate variables)
    # logging.debug(f"All columns after preprocessing: {list(ds.variables.keys())}")

    logging.info("Loading phenotype file...")
    try:
        df = pd.read_csv(pheno_path, sep=",", index_col=0)
    except Exception:
        logging.error("Failed to load phenotype file. Please check your phenotype file path and format.")
        raise RuntimeError("Failed to load phenotype file.") from None

    df.index.name = "samples"
    ds_annotations = df.to_xarray()
    try:

        ds = ds.set_index({"samples": "sample_id"})
        ds = ds.merge(ds_annotations, join="left")
        ds = ds.reset_index("samples").reset_coords("samples")
        ds = ds.rename_vars({"samples": "sample_id"})
    except Exception:
        logging.error("Failed to merge phenotype data with genotype data.")
        raise RuntimeError("Failed to merge phenotype data with genotype data.") from None

    # 计算等位基因剂量（后续GWAS用）
    # Calculate allele dosage (for downstream GWAS)
    ds["call_dosage"] = ds.call_genotype.sum(dim="ploidy")
    logging.info("Data preprocessing finished.")
    return ds
