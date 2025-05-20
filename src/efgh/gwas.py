import os
import sgkit as sg
import pandas as pd
import numpy as np

def run_gwas(ds, config):
    """
    执行GWAS分析，使用sgkit进行基因组关联分析。
    Run GWAS analysis using sgkit.
    """
    # 创建输出目录（如果不存在）
    # Create output directory if it does not exist
    os.makedirs(config.output.outdir, exist_ok=True)
    gwas_results_path = os.path.join(config.output.outdir, "gwas_results.csv")
    # 从配置文件中读取性状列
    # Read trait column from config
    trait = config.gwas.traits
    # 去除性状列空值
    # Remove samples with missing trait values
    ds = ds.sel(samples=~ds[trait].isnull())



    # 加入PCA协变量
    covariates = [f"sample_pca_projection_{i}" for i in range(config.pca.pcs)]

    for cov in covariates:
        assert cov in ds, f"{cov} not in dataset"

    # 运行GWAS分析
    # Run GWAS analysis
    print("Running GWAS analysis...")
    ds_lr = sg.gwas_linear_regression(
        ds,
        add_intercept=True,
        dosage='call_dosage',  # 剂量变量名称 / dosage variable name
        covariates=covariates,         # 协变量名称列表 / list of covariate names
        traits=[trait]         # 性状变量名称 / trait variable name
    )
    print("GWAS analysis completed.")

    print("Preparing GWAS results file...")
    chr_ = ds_lr['variant_contig_name'].values
    pos = ds_lr['variant_position'].values
    locus = [f"{c}:{p}" for c, p in zip(chr_, pos)]
    alleles = ds_lr['variant_allele'].values
    alleles_str = [str(list(a)) for a in alleles]
    # n = ds_lr['variant_n_called'].values
    n = ds_lr.sizes['samples']
    # ds_lr['call_dosage'] 形状为 (variants, samples)
    call_dosage = ds_lr['call_dosage'].values
    sum_x = np.nansum(call_dosage, axis=1)
    # 获取性状值，顺序与 call_dosage 匹配
    y = ds[trait].values
    # 计算 y_transpose_x
    y_transpose_x = np.nansum(call_dosage * y, axis=1)
    beta = ds_lr['variant_linreg_beta'][:, 0].values
    t_stat = ds_lr['variant_linreg_t_value'][:, 0].values
    p_value = ds_lr['variant_linreg_p_value'][:, 0].values

    df = pd.DataFrame({
        "locus": locus,
        "alleles": alleles_str,
        "n": n,
        "sum_x": sum_x,
        "y_transpose_x": y_transpose_x,
        "beta": beta,
        "t_stat": t_stat,
        "p_value": p_value
    })
    df.to_csv(gwas_results_path, index=False)
    print(f"GWAS results saved to: {gwas_results_path}")

    return ds_lr

