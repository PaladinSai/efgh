"""
PCA分析模块，使用sgkit对数据集进行主成分分析（PCA）。
PCA analysis module, perform PCA on dataset using sgkit.
"""

import sgkit as sg
import logging
from .utils import mask_to_numpy_in_chunks

def run_pca(config, ds):
    """
    执行PCA分析流程，输入数据集，输出带主成分的ds。
    Run PCA workflow, input dataset, output ds with principal components.
    参数:
        config: 配置对象 / config object
        ds: 输入数据集 / input dataset
    返回:
        带主成分的ds / ds with principal components
    """
    try:
        pcs = config.pca.pcs
        chunk_size = getattr(getattr(config, "performance", None), "chunk_size", 10000)

        # 计算等位基因计数
        # Calculate alternate allele counts
        ds_pca = sg.stats.pca.count_call_alternate_alleles(ds)
        variant_mask = (((ds_pca.call_alternate_allele_count < 0).any(dim="samples")) |
                        (ds_pca.call_alternate_allele_count.std(dim="samples") <= 0.0))
        mask = mask_to_numpy_in_chunks(variant_mask, chunk_size)
        ds_pca = ds_pca.sel(variants=~mask)
        ds_pca = sg.pca(ds_pca)
        for i in range(pcs):
            ds[f"sample_pca_projection_{i}"] = ds_pca.sample_pca_projection[:, i]
        logging.info(f"PCA analysis completed. {pcs} principal components added.")
        return ds
    except Exception:
        logging.error("Failed to perform PCA analysis. Please check your input data and configuration.")
        raise RuntimeError("Failed to perform PCA analysis.") from None
