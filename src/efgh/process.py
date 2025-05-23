import sgkit as sg
import pandas as pd

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
    #print(f"完成数据预处理后的全部列：{list(ds.variables.keys())}")

    print("Loading phenotype file...")
    df = pd.read_csv(pheno_path, sep=",", index_col= 0)
    df.index.name = "samples"
    ds_annotations = df.to_xarray()
    ds = ds.set_index({"samples": "sample_id"})
    ds = ds.merge(ds_annotations, join="left")
    ds = ds.reset_index("samples").reset_coords("samples")
    ds = ds.rename_vars({"samples": "sample_id"})
    # 计算等位基因剂量（后续GWAS用）
    ds["call_dosage"] = ds.call_genotype.sum(dim="ploidy")
    return ds