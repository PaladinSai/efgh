import time
import logging
import sgkit as sg
from .pca import run_pca
from .gwas import run_gwas

def run_pca_and_gwas(config, process_path, result_path, step_times):
    """PCA分析和GWAS分析的统一入口"""
    ds = sg.load_dataset(process_path)

    # PCA
    t0 = time.time()
    ds = run_pca(config, ds)
    t1 = time.time()
    step_times["PCA"] = t1 - t0
    logging.info(f"Step 'PCA' finished in {step_times['PCA']:.2f} seconds.")

    # GWAS
    t0 = time.time()
    run_gwas(config, ds, result_path)
    t1 = time.time()
    step_times["GWAS"] = t1 - t0
    logging.info(f"Step 'GWAS' finished in {step_times['GWAS']:.2f} seconds.")
    return result_path