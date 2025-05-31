import click
import time
import importlib.resources
import logging
import sgkit as sg
import os
from .config import load_config, get_default_cli_options, get_cpu_cores
from .merge_phenotype import run_merge
from .qc import run_qc
from .pca import run_pca
from .gwas import run_gwas
from .preprocess import run_process
from .result import run_result
from .process import run_pca_and_gwas

DEFAULT_CONFIG_PKG = "efgh.configs"
DEFAULT_CONFIG_FILE = "default.yaml"

def add_dynamic_options(default_yaml_path):
    """
    动态为click命令添加参数，参数名与yaml配置一致，支持多级（如input.vcf_path）
    Dynamically add click command options, parameter names are consistent with YAML config, support nested (e.g., input.vcf_path)
    """
    def decorator(f):
        options = get_default_cli_options(default_yaml_path)
        for key, default in reversed(list(options.items())):
            param_name = key.replace('.', '_')
            param_cli = '--' + key.replace('.', '-')
            param_kwargs = {
                "default": None,
                "show_default": False,
                "help": f"(override config) default: {default}"
            }
            # 根据类型自动推断 / Automatically infer type
            if isinstance(default, bool):
                param_kwargs["is_flag"] = True
            elif isinstance(default, list):
                param_kwargs["multiple"] = True
            f = click.option(param_cli, param_name, **param_kwargs)(f)
        return f
    return decorator

@click.group()
def cli():
    """efgh 命令行工具 / efgh command line tool"""
    # 初始化 logging / Initialize logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[logging.StreamHandler()]
    )

@cli.command()
@add_dynamic_options(lambda: importlib.resources.files(DEFAULT_CONFIG_PKG) / DEFAULT_CONFIG_FILE)
@click.option('--config', 'user_config', default=None, help="User config yaml file path")
def run(user_config, **kwargs):
    """
    运行GWAS流程，自动加载默认配置。
    Run GWAS workflow, automatically load default config.
    """
    total_start = time.time()
    step_times = {}

    # 提取命令行参数，去除None值 / Extract CLI arguments, remove None values
    cli_args = {}
    for k, v in kwargs.items():
        if v is not None:
            # 支持多级key（如input_vcf_path -> input.vcf_path）/ Support nested keys (e.g., input_vcf_path -> input.vcf_path)
            cli_key = k.replace('_', '.')
            cli_args[cli_key] = v

    t0 = time.time()
    config = load_config(
        default_path=None,
        user_path=user_config,
        cli_args=cli_args,
        default_pkg=DEFAULT_CONFIG_PKG,
        default_file=DEFAULT_CONFIG_FILE
    )
    t1 = time.time()
    step_times["Load config"] = t1 - t0
    logging.info(f"Step 'Load config' finished in {step_times['Load config']:.2f} seconds.")

    # 输出配置信息 / Log configuration info
    logging.info("Configuration loaded. Starting GWAS workflow...")
    import pprint
    logging.info("Current configuration:\n" + pprint.pformat(config.to_dict()))



    # 步骤1：VCF转Zarr / Step 1: VCF to Zarr
    # t0 = time.time()
    # vcz_path = vcf_to_zarr(config)
    # t1 = time.time()
    # step_times["VCF to Zarr"] = t1 - t0
    # logging.info(f"Step 'VCF to Zarr' finished in {step_times['VCF to Zarr']:.2f} seconds.")

    vcz_path = config.input.vcz_path
    temp_path = os.path.join(config.output.outdir, "temp")
    try:
        os.makedirs(temp_path, exist_ok=True)
    except Exception:
        logging.error("Failed to create output directory. Please check your output path settings.")
        raise RuntimeError("Failed to create output directory.") from None
    process_path = os.path.join(temp_path, "process.vcz")
    result_path = os.path.join(temp_path, "result.vcz")

    if not os.path.exists(process_path):
        try:
            ds = sg.load_dataset(vcz_path)
        except Exception:
            logging.error("Failed to load Zarr dataset. Please check your Zarr file path and format.")
            raise RuntimeError("Failed to load Zarr dataset.") from None

        # 合并表型文件
        t0 = time.time()
        ds = run_merge(config, ds)
        t1 = time.time()
        step_times["merge"] = t1 - t0
        logging.info(f"Step 'merge' finished in {step_times['merge']:.2f} seconds.")

        # 质量控制 / Quality Control
        t0 = time.time()
        ds = run_qc(config, ds)
        t1 = time.time()
        step_times["QC"] = t1 - t0
        logging.info(f"Step 'QC' finished in {step_times['QC']:.2f} seconds.")

        # 数据处理
        t0 = time.time()
        run_process(config, ds, process_path)
        t1 = time.time()
        step_times["Process"] = t1 - t0
        logging.info(f"Data preprocessing completed in {step_times['Process']:.2f} seconds.")

    else:
        logging.info("Data preprocessing completed, proceeding directly to analysis workflow.")
    if not os.path.exists(result_path):
        # 根据 performance.hpc 启动 Dask 集群
        if getattr(getattr(config, "performance", None), "hpc", False):
            from dask.distributed import Client
            from dask_jobqueue import SLURMCluster
            perf = config.performance
            cpu_cores = get_cpu_cores(config)
            memory = getattr(getattr(config.performance, "memory", None), "memory", "16GB")
            processes = getattr(perf, "processes", 1)
            walltime = getattr(perf, "walltime", "24:00:00")
            job_extra = getattr(perf, "job_extra", ["--exclusive"])
            jobs = getattr(perf, "jobs", 10)
            if isinstance(memory, str):
                mem_str = memory
            else:
                mem_str = "16GB"
            cluster = SLURMCluster(
                cores=cpu_cores,
                memory=mem_str,
                processes=processes,
                walltime=walltime,
                job_extra_directives=job_extra,
            )
            cluster.scale(jobs=jobs)
            with Client(cluster) as client:
                logging.info(f"Dask HPC cluster started with {jobs} jobs, {cpu_cores} cores, and {mem_str} memory per job.")
                run_pca_and_gwas(config, process_path, result_path, step_times)
        else:
            logging.info("HPC cluster not enabled, running in local mode.")
            run_pca_and_gwas(config, process_path, result_path, step_times)
    else:
        logging.info("GWAS results already exist, skipping analysis step.")

    # csv结果生成，画曼哈顿图和QQ图
    t0 = time.time()
    run_result(config, result_path)
    t1 = time.time()
    step_times["Result"] = t1 - t0
    logging.info(f"Results processed and saved in {config.output.outdir}.")

    total_end = time.time()
    for step in ["merge", "QC", "PCA", "Process", "GWAS", "Result"]:
        if step in step_times:
            logging.info(f"Step '{step}' finished in {step_times[step]:.2f} seconds.")
        else:
            logging.info(f"Step '{step}' skipped.")
    logging.info(f"Total time: {total_end - total_start:.2f} seconds")

if __name__ == "__main__":
    cli()
