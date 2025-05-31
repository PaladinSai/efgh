import logging
import os
import pandas as pd
import sgkit as sg
from .plotting import manhattan_plot_chunked, qq_plot

def run_result(config, result_path):
    """
    运行结果处理流程，输入处理后的zarr，输出结果文件。
    Run result processing workflow, input processed zarr, output result file.

    参数:
        config: 配置对象 / config object
        result_path: 处理后的zarr文件路径 / processed zarr file path

    返回:
        结果文件存储路径 / path to result file storage
    """
    # 从配置对象获取参数 / Get parameters from config object
    traits = config.gwas.traits
    covariates = config.gwas.covariates
    chunk_size = getattr(getattr(config, "performance", None), "chunk_size", 10000)
    ds = sg.load_dataset(result_path)

    n_variants = ds.sizes["variants"]

    # 预先以 dask array 提取所有变量，避免 run_spec 警告
    contig_name_da = ds['variant_contig_name'].data
    pos_da = ds['variant_position'].data
    alleles_da = ds['variant_allele'].data
    call_dosage_da = ds['call_dosage'].data
    beta_da = ds['variant_linreg_beta'].data
    t_stat_da = ds['variant_linreg_t_value'].data
    p_value_da = ds['variant_linreg_p_value'].data

    n = ds.sizes['samples']

    for i, trait in enumerate(traits):
        gwas_lr_results_path = os.path.join(config.output.outdir, f"gwas_results_linear_regression_{trait}.csv")
        # 分析前若文件已存在则删除，防止重复写入
        if os.path.exists(gwas_lr_results_path):
            os.remove(gwas_lr_results_path)
        header_written = False
        # trait变量y一次性load（样本数一般都不大）
        y = ds[trait].values
        for start in range(0, n_variants, chunk_size // 2):
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
            manhattan_plot_chunked(ds, config, trait, i, chunk_size=chunk_size // 2)
            logging.info("Manhattan plot generated.")
        except Exception:
            logging.error("Failed to generate Manhattan plot.")
        logging.info("Generating QQ plot...")
        try:
            qq_plot(ds, config, trait, i)
            logging.info("QQ plot generated.")
        except Exception:
            logging.error("Failed to generate QQ plot.")