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

    # print(f"call_dosage 有效值范围: [{ds['call_dosage'].min().compute()}, {ds['call_dosage'].max().compute()}]")
    # print(f"使用的性状列: {trait}")
    # print(f"数据集中性状列: {list(ds.data_vars)}")
    # # 检查samples维度大小
    # print(f"Samples维度大小: {ds.sizes['samples']}")
    # print(f"Variants维度大小: {ds.sizes['variants']}")
    # print(f"call_dosage 形状: {ds['call_dosage'].shape}")
    # print(f"call_dosage 有效值范围: [{ds['call_dosage'].min().compute()}, {ds['call_dosage'].max().compute()}]")
    # out_test_csv(ds, config, "data_inspection.csv")

    # 运行GWAS分析
    # Run GWAS analysis
    print("Running GWAS analysis...")
    ds_lr = sg.gwas_linear_regression(
        ds,
        add_intercept=True,
        dosage='call_dosage',  # 剂量变量名称 / dosage variable name
        covariates=[],         # 协变量名称列表 / list of covariate names
        traits=[trait]         # 性状变量名称 / trait variable name
    )
    print("GWAS analysis completed.")

    # out_test_csv(ds, config, "data_result.csv")

    selected_vars = [
        "variant_contig_name",
        "variant_contig",
        "variant_position",
        "variant_linreg_p_value"
    ]
    df = ds_lr[selected_vars].to_dataframe()
    df.to_csv(gwas_results_path)

    return ds_lr

