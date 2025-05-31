import sgkit as sg
import logging

def run_gwas(config, process_path, result_path):
    """
    执行GWAS分析，使用sgkit进行基因组关联分析。
    Run GWAS analysis using sgkit.
    """
    # 从配置文件中读取配置
    traits = config.gwas.traits
    models = config.gwas.models
    user_covariates = getattr(config.gwas, "covariates", None)
    logging.info(f"Running GWAS analysis using {traits} traits.")

    # 组合自定义协变量和PCA主成分
    if user_covariates and len(user_covariates) > 0:
        if isinstance(user_covariates, str):
            user_covariates = [user_covariates]
    else:
        user_covariates = []
    pca_covariates = [f"sample_pca_projection_{i}" for i in range(config.pca.pcs)]
    covariates = list(dict.fromkeys(user_covariates + pca_covariates))  # 保持顺序去重

    ds = sg.load_dataset(process_path)

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
                sg.save_dataset(ds_lr, config.output.outdir, overwrite=True)
            except Exception:
                logging.error("Failed to run linear regression GWAS. Please check your input data and configuration.")
                raise RuntimeError("Failed to run linear regression GWAS.") from None
            logging.info("Linear regression GWAS completed successfully")
        # elif model == "regenie":
        #     # 这段有点疑问，暂时不提供支持
        #     try:
        #         ds_rg = sg.regenie(ds, add_intercept=True, dosage='call_dosage', covariates=covariates, traits=traits)
        #     except Exception:
        #         logging.error("Failed to run regenie GWAS. Please check your input data and configuration.")
        #         raise RuntimeError("Failed to run regenie GWAS.") from None
        #     for i, trait in enumerate(traits):
        #         gwas_rg_results_path = os.path.join(config.output.outdir, f"gwas_results_regenie_{trait}.csv")
        #         if os.path.exists(gwas_rg_results_path):
        #             os.remove(gwas_rg_results_path)
        #         meta_pred = ds_rg['regenie_meta_prediction'].values
        #         df = pd.DataFrame({
        #             "sample": ds.samples.values,
        #             "meta_prediction": meta_pred[:, i] if meta_pred.ndim == 2 else meta_pred
        #         })
        #         try:
        #             df.to_csv(gwas_rg_results_path, index=False)
        #         except Exception:
        #             logging.error("Failed to save regenie GWAS results to CSV file.")
        #             raise RuntimeError("Failed to save regenie GWAS results.") from None
        #         logging.info(f"GWAS regenie meta prediction saved to: {gwas_rg_results_path}")
        else:
            logging.error(f"Unknown GWAS model: {model}")
            raise RuntimeError(f"Unknown GWAS model: {model}")