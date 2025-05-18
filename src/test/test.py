import os
import zarr
import numpy as np

def zarr_to_vcf(zarr_path, out_path):
    """
    将QC后zarr文件转换为VCF格式
    """
    z = zarr.open(zarr_path, mode='r')
    geno = z['genotype'][:]  # (num_snps, num_samples)
    snp_info = z['snp_info'][:]  # (num_snps, 5) [chrom, pos, vid, ref, alt]
    sample_ids = z['sample_ids'][:]  # (num_samples,)
    num_snps, num_samples = geno.shape

    # 写VCF文件
    with open(out_path, "w") as f:
        # 写VCF头部
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=efgh_qc\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(sample_ids)
        f.write("\t".join(header) + "\n")
        # 写每个SNP
        for i in range(num_snps):
            chrom, pos, vid, ref, alt = snp_info[i]
            # 生成每个样本的GT字段
            gt_row = []
            for j in range(num_samples):
                gt = geno[i, j]
                if gt == -1:
                    gt_str = "./."
                elif gt == 0:
                    gt_str = "0/0"
                elif gt == 1:
                    gt_str = "0/1"
                elif gt == 2:
                    gt_str = "1/1"
                else:
                    gt_str = "./."
                gt_row.append(gt_str)
            row = [chrom, str(pos), vid, ref, alt, ".", "PASS", ".", "GT"] + gt_row
            f.write("\t".join(row) + "\n")

def compare_vcf(file1, file2, diff_out="vcf_diff.txt"):
    """
    对比两个VCF文件，输出不同的SNP或基因型信息到diff_out文件。
    只比较CHROM、POS、ID、REF、ALT和GT字段。
    """
    def parse_vcf(path):
        snp_dict = {}
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                key = tuple(fields[:5])  # (CHROM, POS, ID, REF, ALT)
                gt = tuple(fields[9:])   # 所有样本的GT
                snp_dict[key] = gt
        return snp_dict

    vcf1 = parse_vcf(file1)
    vcf2 = parse_vcf(file2)
    diff_lines = []

    # 找出file1有但file2没有的SNP
    for key in vcf1:
        if key not in vcf2:
            diff_lines.append(f"Only in {file1}: {key}")
        elif vcf1[key] != vcf2[key]:
            diff_lines.append(f"Genotype diff at {key}:\n  {file1}: {vcf1[key]}\n  {file2}: {vcf2[key]}")
    # 找出file2有但file1没有的SNP
    for key in vcf2:
        if key not in vcf1:
            diff_lines.append(f"Only in {file2}: {key}")

    with open(diff_out, "w") as f:
        for line in diff_lines:
            f.write(str(line) + "\n")
    print(f"VCF对比完成，差异输出到 {diff_out}")

if __name__ == "__main__":
    # 示例用法
    # 假设zarr文件路径和输出路径
    zarr_path = "../../results/CD8_qc_only.zarr"
    out_path = "../../results/CD8_qc_only.vcf"
    zarr_to_vcf(zarr_path, out_path)
    print("Zarr转VCF完成。")
    # 比较两个VCF文件
    # compare_vcf("../../results/CD8_qc_only.vcf", "../../results/other.vcf")
    pass
