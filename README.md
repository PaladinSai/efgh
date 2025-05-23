# efGH: Efficient GWAS Pipeline

## 项目简介 | Project Introduction

efGH 是一个简单且高效的全基因组关联分析（GWAS）软件，支持VCF文件到Zarr格式的转换，并集成了sgkit等工具，适用于基因型数据的分析和处理。

## 主要功能 | Main Features

- VCF文件到Zarr格式的自动转换
- 集成bio2zarr和sgkit等主流基因组分析工具
- 简单易用的命令行接口

## 快速开始 | Quick Start

1. 安装依赖 pip install efgh（示例，暂未提供此方法）
2. 使用命令行运行主流程脚本，输入VCF文件，输出gwas分析结果图
```shell
# 使用默认配置（将使用configs/default.yaml文件中的默认配置运行）
efgh run
# 指定配置文件（推荐）
efgh run --config your_config.yaml
```
3. 运行完成后，结果文件将保存在指定的输出目录中。
4. 配置文件参考configs/default.yaml，用户可以根据需要修改配置文件中的参数。

## 适用场景 | Application Scenarios

- 基因型数据的预处理
- 基因型-表型关联分析
- 生物信息学高性能计算流程

## 联系方式 | Contact

如需帮助或有建议，请联系项目维护者。

---

# efGH: Efficient GWAS Pipeline (English Version)

## Project Introduction

efGH is a simple and efficient Genome-Wide Association Study (GWAS) software that supports conversion from VCF files to Zarr format. It integrates tools such as sgkit and is suitable for the analysis and processing of genotype data.

## Main Features

- Automatic conversion from VCF files to Zarr format
- Integration with mainstream genomic analysis tools such as bio2zarr and sgkit
- Easy-to-use command-line interface

## Quick Start

1. Install dependencies: `pip install efgh` (example, this method is not yet available)
2. Run the main workflow script from the command line, input VCF files, and output GWAS analysis result plots
```shell
# Use default configuration (runs with configs/default.yaml)
efgh run
# Specify a configuration file (recommended)
efgh run --config your_config.yaml
```
3. After completion, the result files will be saved in the specified output directory.
4. Refer to configs/default.yaml for configuration file format. Users can modify parameters as needed.

## Application Scenarios

- Preprocessing of genotype data
- Genotype-phenotype association analysis
- High-performance computing workflows in bioinformatics

## Contact

For help or suggestions, please contact the project maintainer.

`