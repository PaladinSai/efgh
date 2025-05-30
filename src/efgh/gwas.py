import os
import sgkit as sg
import pandas as pd
import numpy as np
import dask.array as da
import logging
from .plotting import manhattan_plot_chunked, qq_plot
from .utils import mask_to_numpy_in_chunks

def run_gwas(config, ds):
    """
    执行GWAS分析，使用sgkit进行基因组关联分析。
    Run GWAS analysis using sgkit.
    """
    # 创建输出目录（如果不存在）
    try:
        os.makedirs(config.output.outdir, exist_ok=True)
    except Exception:
        logging.error("Failed to create output directory. Please check your output path settings.")
        raise RuntimeError("Failed to create output directory.") from None

    # 从配置文件中读取配置
    chunk_size = getattr(getattr(config, "performance", None), "chunk_size", 10000)
    traits = config.gwas.traits
    models = config.gwas.models
    user_covariates = getattr(config.gwas, "covariates", None)

    # 组合自定义协变量和PCA主成分
    if user_covariates and len(user_covariates) > 0:
        if isinstance(user_covariates, str):
            user_covariates = [user_covariates]
    else:
        user_covariates = []
    pca_covariates = [f"sample_pca_projection_{i}" for i in range(config.pca.pcs)]
    covariates = list(dict.fromkeys(user_covariates + pca_covariates))  # 保持顺序去重

    logging.info(f"Running GWAS analysis using {traits} traits.")
    # 去除性状列和协变量列的空值
    mask = da.ones(ds.sizes["samples"], dtype=bool)
    for col in list(traits) + list(user_covariates):
        mask &= ~ds[col].isnull().persist()
    mask_np = mask_to_numpy_in_chunks(mask, chunk_size)
    ds = ds.sel(samples=mask_np)

    for model in models:
        logging.info(f"Running GWAS model: {model} ...")
        if model == "linear_regression":
            try:
                ds_lr = sg.gwas_linear_regression(
                    ds,
                    add_intercept=True,
                    dosage='call_dosage',
                    covariates=covariates,
                    traits=traits
                )
            except Exception:
                logging.error("Failed to run linear regression GWAS. Please check your input data and configuration.")
                raise RuntimeError("Failed to run linear regression GWAS.") from None
            logging.info("Linear regression GWAS completed successfully. Start saving results...")

            n_variants = ds_lr.sizes["variants"]

            # 预先以 dask array 提取所有变量，避免 run_spec 警告
            contig_name_da = ds_lr['variant_contig_name'].data
            pos_da = ds_lr['variant_position'].data
            alleles_da = ds_lr['variant_allele'].data
            call_dosage_da = ds_lr['call_dosage'].data
            beta_da = ds_lr['variant_linreg_beta'].data
            t_stat_da = ds_lr['variant_linreg_t_value'].data
            p_value_da = ds_lr['variant_linreg_p_value'].data

            n = ds_lr.sizes['samples']

            for i, trait in enumerate(traits):
                gwas_lr_results_path = os.path.join(config.output.outdir, f"gwas_results_linear_regression_{trait}.csv")
                # 分析前若文件已存在则删除，防止重复写入
                if os.path.exists(gwas_lr_results_path):
                    os.remove(gwas_lr_results_path)
                header_written = False
                # trait变量y一次性load（样本数一般都不大）
                y = ds[trait].values
                for start in range(0, n_variants, chunk_size//2):
                    end = min(start + chunk_size, n_variants)
                    idx = slice(start, end)
                    chr_ = contig_name_da[idx].compute()
                    pos = pos_da[idx].compute()
                    locus = [f"{c}:{p}" for c, p in zip(chr_, pos)]
                    alleles = alleles_da[idx].compute()
                    alleles_str = [str(list(a)) for a in alleles]
                    call_dosage = call_dosage_da[idx].compute()
                    sum_x = np.nansum(call_dosage, axis=1)
                    y_transpose_x = np.nansum(call_dosage * y, axis=1)
                    beta = beta_da[idx, i].compute()
                    t_stat = t_stat_da[idx, i].compute()
                    p_value = p_value_da[idx, i].compute()
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
                    manhattan_plot_chunked(ds_lr, config, trait, i, chunk_size=chunk_size//2)
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
            # 这段有点疑问，暂时不提供支持
            try:
                ds_rg = sg.regenie(ds, add_intercept=True, dosage='call_dosage', covariates=covariates, traits=traits)
            except Exception:
                logging.error("Failed to run regenie GWAS. Please check your input data and configuration.")
                raise RuntimeError("Failed to run regenie GWAS.") from None
            for i, trait in enumerate(traits):
                gwas_rg_results_path = os.path.join(config.output.outdir, f"gwas_results_regenie_{trait}.csv")
                if os.path.exists(gwas_rg_results_path):
                    os.remove(gwas_rg_results_path)
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