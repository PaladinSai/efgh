"""
VCF转Zarr模块，通过命令行调用bio2zarr CLI工具，将VCF文件转换为zarr格式，供sgkit后续分析使用。
VCF to Zarr module. Uses the bio2zarr CLI tool to convert VCF files to Zarr format for downstream sgkit analysis.
"""

import os
import subprocess
import pysam

def create_vcf_index(vcf_path):
    """
    使用pysam为VCF文件生成索引（.tbi），如果已存在则跳过。
    Create index (.tbi) for VCF file using pysam. Skip if already exists.
    """
    # 判断是否为bgzip压缩VCF
    # Check if VCF is bgzip compressed
    if vcf_path.endswith(('.gz', '.bgz')):
        index_path = vcf_path + '.tbi'
        if os.path.exists(index_path):
            print(f"VCF index already exists: {index_path}")
            return
        print(f"Creating VCF index: {index_path}")
        pysam.tabix_index(vcf_path, preset="vcf", force=True)
    else:
        raise RuntimeError("VCF file must be bgzip compressed (.vcf.gz) to create index.")

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

    # 先生成VCF索引
    # Generate VCF index first
    create_vcf_index(vcf_path)

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

