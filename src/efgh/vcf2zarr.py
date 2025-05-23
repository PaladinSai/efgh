"""
VCF转Zarr模块，通过命令行调用bio2zarr CLI工具，将VCF文件转换为zarr格式，供sgkit后续分析使用。
VCF to Zarr module. Uses the bio2zarr CLI tool to convert VCF files to Zarr format for downstream sgkit analysis.
"""

import os
import subprocess
import pysam

def vcf_to_zarr(config):
    """
    使用bio2zarr CLI将VCF文件转换为zarr格式。
    Convert VCF file to Zarr format using bio2zarr CLI.
    参数:
        config: Config对象，包含所有配置参数
        config: Config object containing all configuration parameters
    返回zarr文件路径
    Returns the path to the zarr file.
    """
    vcf_path = config.input.vcf_path
    outdir = os.path.join(config.output.outdir, "temp")
    basename = "genotype_raw"
    icf_path = os.path.join(outdir, f"{basename}.icf")
    zarr_path = os.path.join(outdir, f"{basename}.vcz")

    # 创建输出目录
    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # 检查VCF索引是否存在
    # Check if VCF index exists
    index_path = vcf_path + ".tbi"
    if not os.path.exists(index_path):
        raise RuntimeError(
            f"VCF index file not found: {index_path}\n"
            f"Please compress your VCF with bgzip and create the index with tabix first, e.g.:\n"
            f"  bgzip {vcf_path.rstrip('.gz')}\n"
            f"  tabix -p vcf {vcf_path}\n"
        )

    # 如果zarr文件已存在，直接返回
    # If zarr file already exists, return directly
    if os.path.exists(zarr_path):
        print(f"Zarr file already exists: {zarr_path}")
        return zarr_path

    # 1. vcf2zarr explode --force vcf icf
    print(f"Running: vcf2zarr explode --force {vcf_path} {icf_path}")
    subprocess.run(
        ["vcf2zarr", "explode", "--force", vcf_path, icf_path],
        check=True
    )

    # 2. vcf2zarr encode --force icf vcz
    print(f"Running: vcf2zarr encode --force {icf_path} {zarr_path}")
    subprocess.run(
        ["vcf2zarr", "encode", "--force", icf_path, zarr_path],
        check=True
    )
    print(f"VCF to Zarr conversion completed. Output path: {zarr_path}")
    return zarr_path

