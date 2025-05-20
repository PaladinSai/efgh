import os
import sgkit as sg
import pandas as pd
from .util import out_test_csv

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

    n_samples = ds.sizes['samples']
    chr_ = ds_lr['variant_contig_name'].values
    rs = ds_lr['variant_id'].values
    ps = ds_lr['variant_position'].values
    n_obs = ds_lr['variant_n_called'].values
    n_mis = n_samples - n_obs
    af = ds_lr['variant_allele_frequency'].values
    beta = ds_lr['variant_linreg_beta'][:, 0].values
    #se = ds_lr['variant_linreg_standard_error'][:, 0].values
    p_wald = ds_lr['variant_linreg_p_value'][:, 0].values

    # 动态处理等位基因
    alleles = ds_lr['variant_allele'].values  # shape: (variants, alleles)
    max_alleles = alleles.shape[1]
    allele_cols = {}
    for i in range(max_alleles):
        allele_cols[f'allele{i}'] = alleles[:, i]

    # 动态处理等位基因频率
    af_cols = {}
    for i in range(max_alleles):
        af_cols[f'af{i}'] = af[:, i]

    # 组装DataFrame
    data = {
        "chr": chr_,
        "rs": rs,
        "ps": ps,
        "n_mis": n_mis,
        "n_obs": n_obs,
        **allele_cols,
        **af_cols,
        "beta": beta,
        #"se": se,
        "p_wald": p_wald
    }
    df = pd.DataFrame(data)
    df.to_csv(gwas_results_path, index=False)

    return ds_lr

