import dask.array as da
import sgkit as sg

def run_process(config, ds, process_path):
    """
    运行数据处理流程，输入zarr，输出处理后的ds。
    Run data processing workflow, input zarr, output processed ds.

    参数:
        config: 配置对象 / config object
        ds: 输入数据集 / input dataset
        process_path: 处理后的zarr文件路径 / processed zarr file path

    返回:
        处理后文件存储路径 / path to processed file storage
    """
    # 从配置对象获取参数 / Get parameters from config object
    traits = config.gwas.traits
    covariates = config.gwas.covariates

    ds["variant_contig_name"] = ds.contig_id[ds.variant_contig]
    # 计算等位基因剂量（后续GWAS用）
    ds["call_dosage"] = ds.call_genotype.sum(dim="ploidy")

    # 去除性状列和协变量列的空值
    mask = da.ones(ds.sizes["samples"], dtype=bool)
    for col in list(traits) + list(covariates):
        mask &= ~ds[col].isnull().persist()
    ds = ds.sel(samples=mask.compute())
    sg.save_dataset(ds, process_path)
    return process_path