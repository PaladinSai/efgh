import sgkit as sg
import pandas as pd
import xarray as xr

def run_process(config, vcz_path):
    """
    运行数据处理流程，输入zarr，输出处理后的ds。
    Run data processing workflow, input zarr, output processed ds.
    """
    # 从配置对象获取参数
    pheno_path = config.input.pheno_path


    ds = sg.load_dataset(vcz_path)

    print("Data preprocessing")
    ds["variant_contig_name"] = ds.contig_id[ds.variant_contig]

    # 列出所有变量（包括坐标变量）
    print(f"完成数据预处理后的全部列：{list(ds.variables.keys())}")
    print(ds)
    print(sg.display_genotypes(ds, max_variants=10, max_samples=5))

    print("Loading phenotype file...")
    df = pd.read_csv(pheno_path, sep=",", index_col= 0)
    df.index.name = "samples"
    ds_annotations = df.to_xarray()

    print(f"表型文件的全部列：{list(ds_annotations.variables.keys())}")
    print(ds_annotations)

    ds = ds.set_index({"samples": "sample_id"})

    print(f"合并前ds的全部列：{list(ds.variables.keys())}")
    print(ds)

    ds = ds.merge(ds_annotations, join="left")
    ds = ds.reset_index("samples").reset_coords("samples")
    ds = ds.rename_vars({"samples": "sample_id"})
    # 计算等位基因剂量（后续GWAS用）
    ds["call_dosage"] = ds.call_genotype.sum(dim="ploidy")

    # 列出所有变量（包括坐标变量）
    print(f"加载表型文件后的全部列：{list(ds.variables.keys())}")
    print(ds)
    print(sg.display_genotypes(ds, max_variants=10, max_samples=5))
    return ds