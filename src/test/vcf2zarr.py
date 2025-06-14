"""
VCF转Zarr模块，通过命令行调用bio2zarr CLI工具，将VCF文件转换为zarr格式，供sgkit后续分析使用。
VCF to Zarr module. Uses the bio2zarr CLI tool to convert VCF files to Zarr format for downstream sgkit analysis.
"""

import os
import subprocess
import logging
import pysam

def vcf_to_zarr(vcf_path, result_path):
    """
    使用bio2zarr CLI将VCF文件转换为zarr格式，支持多进程。
    Convert VCF file to Zarr format using bio2zarr CLI, support multi-process.
    参数:
        config: Config对象，包含所有配置参数
        config: Config object containing all configuration parameters
    返回zarr文件路径
    Returns the path to the zarr file.
    """
    basename = "test"
    icf_path = os.path.join(result_path, f"{basename}.icf")
    zarr_path = os.path.join(result_path, f"{basename}.vcz")

    # 获取CPU核心数，调用config中的方法
    # Get CPU core count from config utility
    cpu_cores = 8

    # 创建输出目录
    # Create output directory
    try:
        os.makedirs(result_path, exist_ok=True)
    except Exception:
        logging.error("Failed to create output directory for Zarr conversion.")
        #raise RuntimeError("Failed to create output directory for Zarr conversion.") from None

    # 检查VCF索引是否存在
    # Check if VCF index (.tbi) exists
    index_path = vcf_path + ".tbi"
    if not os.path.exists(index_path):
        pysam.tabix_compress(vcf_path, f"{vcf_path}.gz", force=True)
        pysam.tabix_index(f"{vcf_path}.gz", preset="vcf", force=True)

    # 如果zarr文件已存在，直接返回
    # If zarr file already exists, return directly
    if os.path.exists(zarr_path):
        logging.info(f"Zarr file already exists: {zarr_path}")
        return zarr_path

    # 1. vcf2zarr explode --force --worker-processes N vcf icf
    # 1. vcf2zarr explode --force --worker-processes N vcf icf
    try:
        logging.info(f"Running: vcf2zarr explode --force --worker-processes {cpu_cores} {vcf_path} {icf_path}")
        subprocess.run(
            ["vcf2zarr", "explode", "--force", "--worker-processes", str(cpu_cores), vcf_path, icf_path],
            check=True
        )
    except Exception:
        logging.error("Failed to run 'vcf2zarr explode'. Please check your VCF file and bio2zarr installation.")
        #raise RuntimeError("Failed to run 'vcf2zarr explode'.") from None

    # 2. vcf2zarr encode --force --worker-processes N icf vcz
    # 2. vcf2zarr encode --force --worker-processes N icf vcz
    try:
        logging.info(f"Running: vcf2zarr encode --force --worker-processes {cpu_cores} {icf_path} {zarr_path}")
        subprocess.run(
            ["vcf2zarr", "encode", "--force", "--worker-processes", str(cpu_cores), icf_path, zarr_path],
            check=True
        )
    except Exception:
        logging.error("Failed to run 'vcf2zarr encode'. Please check your intermediate files and bio2zarr installation.")
        raise RuntimeError("Failed to run 'vcf2zarr encode'.") from None

    logging.info(f"VCF to Zarr conversion completed. Output path: {zarr_path}")
    return zarr_path


vcf_to_zarr("/mnt/g/code/python/efgh/data/test.vcf.gz", "/mnt/g/code/python/efgh/data")