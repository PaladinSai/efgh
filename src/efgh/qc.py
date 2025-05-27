"""
质量控制模块，使用sgkit对zarr文件进行QC，输出QC后的zarr文件。
Quality control module, use sgkit to perform QC on zarr file and output QC'ed zarr.
"""
import sgkit as sg
import logging
import dask.array as da

def run_qc(config, ds):
    """
    执行质量控制流程，输入zarr，输出QC后的zarr。
    Run quality control workflow, input zarr, output QC'ed zarr.
    参数:
        config: Config对象，包含所有配置参数
        config: Config object, contains all configuration parameters
    """
    # 从配置对象获取参数 / Get parameters from config object
    qc_cfg = config.qc

    try:
        # 计算样本统计信息 / Calculate sample statistics
        logging.info("Calculating sample statistics...")
        ds = sg.sample_stats(ds)

        # 根据检出率过滤样本（样本缺失率过滤）/ Filter samples by call rate (missingness)
        if qc_cfg.sample_missing is not None:
            logging.info(f"Filtering samples with call rate less than {qc_cfg.sample_missing} ...")
            logging.info(f"Sample count before filtering: {ds.sizes['samples']}")
            simples_missing_filter = (ds.sample_call_rate > (1 - qc_cfg.sample_missing))
            ds = ds.sel(samples=simples_missing_filter.compute())
            logging.info(f"Sample count after call rate filtering: {ds.sizes['samples']}")

        # 计算变异位点统计信息 / Calculate variant statistics
        logging.info("Calculating variant statistics...")
        ds = sg.variant_stats(ds)

        variant_mask = da.ones(ds.sizes['variants'], dtype=bool)
        # 根据检出率过滤变异位点（变异位点缺失率过滤）/ Filter variants by call rate (missingness)
        if qc_cfg.variant_missing is not None:
            logging.info(f"Filtering variants with call rate less than {qc_cfg.variant_missing} ...")
            variant_mask &= (ds.variant_call_rate > (1 - qc_cfg.variant_missing))
        # 根据 MAF 过滤变异位点 / Filter variants by MAF
        if qc_cfg.maf is not None:
            logging.info("Calculating minor allele frequency (MAF)...")
            variant_mask &= (ds.variant_allele_frequency[:, 1] > qc_cfg.maf)
        # HWE过滤 / HWE filtering
        if qc_cfg.hwe is not None:
            logging.info(f"Filtering variants with HWE less than {qc_cfg.hwe} ...")
            ds = sg.hardy_weinberg_test(ds)
            variant_mask &= (ds.variant_hwe_p_value > float(qc_cfg.hwe))
        logging.info(f"Variant count before filtering: {ds.sizes['variants']}")
        ds = ds.sel(variants=variant_mask.compute())
        logging.info(f"Variant count after filtering: {ds.sizes['variants']}")

        logging.info(f"Sample count after QC: {ds.sizes['samples']}")
        logging.info(f"Variant count after QC: {ds.sizes['variants']}")
        return ds

    except Exception:
        logging.error("Failed during quality control. Please check your input data and QC configuration.")
        raise RuntimeError("Failed during quality control.") from None
