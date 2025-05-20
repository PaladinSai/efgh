"""
PCA分析模块，使用sgkit对数据集进行主成分分析（PCA）。
"""

import sgkit as sg

def run_pca(config,ds):

    pcs = config.pca.pcs

    ds_pca = sg.stats.pca.count_call_alternate_alleles(ds)
    variant_mask = (((ds_pca.call_alternate_allele_count < 0).any(dim="samples")) |
                    (ds_pca.call_alternate_allele_count.std(dim="samples") <= 0.0)).compute()
    ds_pca =  ds_pca.sel(variants=~variant_mask)
    ds_pca = sg.pca(ds_pca)
    for i in range(pcs):
        ds[f"sample_pca_projection_{i}"] = ds_pca.sample_pca_projection[:, i]
    return ds