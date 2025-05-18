import click
from .config import load_config, get_default_cli_options
import time
from .vcf2zarr import vcf_to_zarr
from .process import run_process
from .qc import run_qc
from .pca import run_pca
from .gwas import run_gwas
from .plotting import manhattan_plot, qq_plot


DEFAULT_CONFIG_PATH = "configs/default.yaml"

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
            # 根据类型自动推断
            # Automatically infer type
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
    pass

@cli.command()
@add_dynamic_options(DEFAULT_CONFIG_PATH)
@click.option('--config', 'user_config', default=None, help="User config yaml file path")
def run(user_config, **kwargs):
    """
    运行GWAS流程，自动加载默认配置。
    Run GWAS workflow, automatically load default config.
    """
    total_start = time.time()
    step_times = {}

    # 提取命令行参数，去除None值
    # Extract CLI arguments, remove None values
    cli_args = {}
    for k, v in kwargs.items():
        if v is not None:
            # 支持多级key（如input_vcf_path -> input.vcf_path）
            # Support nested keys (e.g., input_vcf_path -> input.vcf_path)
            cli_key = k.replace('_', '.')
            cli_args[cli_key] = v

    t0 = time.time()
    config = load_config(DEFAULT_CONFIG_PATH, user_config, cli_args)
    t1 = time.time()
    step_times["Load config"] = t1 - t0

    print("Configuration loaded. Starting GWAS workflow...")
    print("Current configuration:")
    import pprint
    pprint.pprint(config.to_dict())

    # 步骤1：VCF转Zarr
    # Step 1: VCF to Zarr
    t0 = time.time()
    vcz_path = vcf_to_zarr(config)
    t1 = time.time()
    step_times["VCF to Zarr"] = t1 - t0

    # 步骤2： 数据预处理
    t0 = time.time()
    ds = run_process(config, vcz_path)
    t1 = time.time()
    step_times["process"] = t1 - t0

    # 步骤3：质量控制
    # Step 3: Quality Control
    t0 = time.time()
    ds = run_qc(config, ds)
    t1 = time.time()
    step_times["QC"] = t1 - t0

    # 步骤4：pca
    t0 = time.time()
    ds = run_pca(config, ds)
    t1 = time.time()
    step_times["PCA"] = t1 - t0

    # 步骤5：gwas
    # Step 5: GWAS
    t0 = time.time()
    ds_lr = run_gwas(ds, config)
    t1 = time.time()
    step_times["GWAS"] = t1 - t0

    # 步骤6：绘图
    # Step 6: Plotting
    t0 = time.time()
    print("Generating Manhattan plot...")
    manhattan_plot(ds_lr, config)
    print("Manhattan plot generated.")
    print("Generating QQ plot...")
    qq_plot(ds_lr, config)
    print("QQ plot generated.")
    t1 = time.time()
    step_times["Plotting"] = t1 - t0

    total_end = time.time()
    print("\nStep timing (seconds):")
    for step, sec in step_times.items():
        print(f"{step}: {sec:.2f}")
    print(f"Total time: {total_end - total_start:.2f} seconds")

if __name__ == "__main__":
    cli()
