"""
质量控制模块，使用sgkit对zarr文件进行QC，输出QC后的zarr文件。
Quality control module, use sgkit to perform QC on zarr file and output QC'ed zarr.
"""
import sgkit as sg

def run_qc(config,ds):
    """
    执行质量控制流程，输入zarr，输出QC后的zarr。
    Run quality control workflow, input zarr, output QC'ed zarr.
    参数:
        config: Config对象，包含所有配置参数
        config: Config object, contains all configuration parameters
    """
    # 从配置对象获取参数
    # Get parameters from config object
    qc_cfg = config.qc



    # 计算样本统计信息
    # Calculate sample statistics
    print("Calculating sample statistics...")
    ds = sg.sample_stats(ds)

    # 根据检出率过滤样本（样本缺失率过滤）
    # Filter samples by call rate (missingness)
    if qc_cfg.sample_missing is not None:
        print(f"Filtering samples with call rate less than {qc_cfg.sample_missing} ...")
        print(f"Sample count before filtering: {ds.sizes['samples']}")
        simples_missing_filter = (ds.sample_call_rate > (1 - qc_cfg.sample_missing))
        ds = ds.sel(samples=simples_missing_filter.compute())
        print(f"Sample count after call rate filtering: {ds.sizes['samples']}")

    # 计算变异位点统计信息
    # Calculate variant statistics
    print("Calculating variant statistics...")
    ds = sg.variant_stats(ds)

    # 根据检出率过滤变异位点（变异位点缺失率过滤）
    # Filter variants by call rate (missingness)
    if qc_cfg.variant_missing is not None:
        print(f"Filtering variants with call rate less than {qc_cfg.variant_missing} ...")
        print(f"Variant count before filtering: {ds.sizes['variants']}")
        variants_missing_filter = (ds.variant_call_rate > (1 - qc_cfg.variant_missing))
        ds = ds.sel(variants=variants_missing_filter.compute())
        print(f"Variant count after call rate filtering: {ds.sizes['variants']}")

    # 根据 MAF 过滤变异位点
    # Filter variants by MAF
    if qc_cfg.maf is not None:
        print("Calculating minor allele frequency (MAF)...")
        print(f"Variant count before MAF filtering: {ds.sizes['variants']}")
        maf_filter = (ds.variant_allele_frequency[:,1] > qc_cfg.maf)
        ds = ds.sel(variants=maf_filter.compute())
        print(f"Variant count after MAF filtering: {ds.sizes['variants']}")

    # HWE过滤
    # HWE filtering
    if qc_cfg.hwe is not None:
        print(f"Variant count before HWE filtering: {ds.sizes['variants']}")
        ds =  sg.hardy_weinberg_test(ds)
        print(f"Filtering variants with HWE less than {qc_cfg.hwe} ...")
        hwe_filter = (ds.variant_hwe_p_value > float(qc_cfg.hwe))
        ds = ds.sel(variants=hwe_filter.compute())
        print(f"Variant count after HWE filtering: {ds.sizes['variants']}")

    # print(f"QC后样本数量: {ds.sizes['samples']}")
    # print(f"QC后变异位点数量: {ds.sizes['variants']}")
    print(f"Sample count after QC: {ds.sizes['samples']}")
    print(f"Variant count after QC: {ds.sizes['variants']}")

    return ds

