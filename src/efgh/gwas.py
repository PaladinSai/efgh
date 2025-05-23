import os
import sgkit as sg
import pandas as pd
import numpy as np
import logging
from .plotting import manhattan_plot, qq_plot

def run_gwas(ds, config, chunk_size=10000):
    """
    执行GWAS分析，使用sgkit进行基因组关联分析。
    Run GWAS analysis using sgkit.
    """
    # 创建输出目录（如果不存在）
    # Create output directory if it does not exist
    try:
        os.makedirs(config.output.outdir, exist_ok=True)
    except Exception:
        logging.error("Failed to create output directory. Please check your output path settings.")
        raise RuntimeError("Failed to create output directory.") from None

    # 从配置文件中读取性状列
    # Read trait column from config
    traits = config.gwas.traits

    # 组合自定义协变量和PCA主成分
    # Combine user covariates and PCA principal components
    user_covariates = getattr(config.gwas, "covariates", None)
    if user_covariates and len(user_covariates) > 0:
        if isinstance(user_covariates, str):
            user_covariates = [user_covariates]
    else:
        user_covariates = []
    pca_covariates = [f"sample_pca_projection_{i}" for i in range(config.pca.pcs)]
    covariates = list(dict.fromkeys(user_covariates + pca_covariates))  # 保持顺序去重 / keep order and deduplicate

    # 去除性状列和协变量列的空值
    # Remove samples with missing trait or covariate values
    mask = np.ones(ds.sizes["samples"], dtype=bool)
    for col in list(traits) + list(covariates):
        mask &= ~ds[col].isnull().values
    ds = ds.sel(samples=mask)

    for cov in covariates:
        if cov not in ds:
            logging.error(f"Covariate '{cov}' not found in dataset.")
            raise RuntimeError(f"Covariate '{cov}' not found in dataset.")

    models = config.gwas.models

    for model in models:
        logging.info(f"Running GWAS model: {model} ...")
        if model == "linear_regression":
            try:
                ds_lr = sg.gwas_linear_regression(
                    ds,
                    add_intercept=True,
                    dosage='call_dosage',  # 剂量变量名称 / dosage variable name
                    covariates=covariates,  # 协变量名称列表 / list of covariate names
                    traits=traits  # 性状变量名称 / trait variable name
                )
            except Exception:
                logging.error("Failed to run linear regression GWAS. Please check your input data and configuration.")
                raise RuntimeError("Failed to run linear regression GWAS.") from None

            n_variants = ds_lr.dims["variants"]
            for i, trait in enumerate(traits):
                gwas_lr_results_path = os.path.join(config.output.outdir, f"gwas_results_linear_regression_{trait}.csv")
                header_written = False
                for start in range(0, n_variants, chunk_size):
                    end = min(start + chunk_size, n_variants)
                    chr_ = ds_lr['variant_contig_name'].values[start:end]
                    pos = ds_lr['variant_position'].values[start:end]
                    locus = [f"{c}:{p}" for c, p in zip(chr_, pos)]
                    alleles = ds_lr['variant_allele'].values[start:end]
                    alleles_str = [str(list(a)) for a in alleles]
                    n = ds_lr.sizes['samples']
                    call_dosage = ds_lr['call_dosage'].values[start:end]
                    sum_x = np.nansum(call_dosage, axis=1)
                    y = ds[trait].values
                    y_transpose_x = np.nansum(call_dosage * y, axis=1)
                    beta = ds_lr['variant_linreg_beta'][start:end, 0].values
                    t_stat = ds_lr['variant_linreg_t_value'][start:end, 0].values
                    p_value = ds_lr['variant_linreg_p_value'][start:end, 0].values
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
                    try:
                        df.to_csv(
                            gwas_lr_results_path,
                            mode="a",
                            header=not header_written,
                            index=False
                        )
                    except Exception:
                        logging.error("Failed to save GWAS results to CSV file.")
                        raise RuntimeError("Failed to save GWAS results.") from None
                    header_written = True
                logging.info(f"GWAS results saved to: {gwas_lr_results_path}")
                logging.info("Generating Manhattan plot...")
                try:
                    manhattan_plot(ds_lr, config, trait, i)
                    logging.info("Manhattan plot generated.")
                except Exception:
                    logging.error("Failed to generate Manhattan plot.")
                logging.info("Generating QQ plot...")
                try:
                    qq_plot(ds_lr, config, trait, i)
                    logging.info("QQ plot generated.")
                except Exception:
                    logging.error("Failed to generate QQ plot.")
        elif model == "regenie":
            try:
                ds_rg = sg.regenie(ds, add_intercept=True, dosage='call_dosage', covariates=covariates, traits=traits)
            except Exception:
                logging.error("Failed to run regenie GWAS. Please check your input data and configuration.")
                raise RuntimeError("Failed to run regenie GWAS.") from None
            for i, trait in enumerate(traits):
                gwas_rg_results_path = os.path.join(config.output.outdir, f"gwas_results_regenie_{trait}.csv")
                meta_pred = ds_rg['regenie_meta_prediction'].values
                df = pd.DataFrame({
                    "sample": ds.samples.values,
                    "meta_prediction": meta_pred[:, i] if meta_pred.ndim == 2 else meta_pred
                })
                try:
                    df.to_csv(gwas_rg_results_path, index=False)
                except Exception:
                    logging.error("Failed to save regenie GWAS results to CSV file.")
                    raise RuntimeError("Failed to save regenie GWAS results.") from None
                logging.info(f"GWAS regenie meta prediction saved to: {gwas_rg_results_path}")
        else:
            logging.error(f"Unknown GWAS model: {model}")
            raise RuntimeError(f"Unknown GWAS model: {model}")

