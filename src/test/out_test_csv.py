import os
import pandas as pd

def out_csv(ds, config, path):
    # 输出所有变量的前10行数据到CSV文件
    print("正在导出数据集中的所有变量的前10行用于检查...")
    inspection_path = os.path.join(config.output.outdir,path)

    # 创建一个空的DataFrame来存储所有变量
    all_data = {}

    # 循环遍历所有数据变量
    for var_name in ds.data_vars:
        try:
            # 获取变量的维度信息
            dims = ds[var_name].dims
            print(f"变量 {var_name} 的维度: {dims}")

            # 对于不同维度的变量，我们可能需要不同的处理方式
            if len(dims) == 1:  # 一维变量
                # 取前10个值
                values = ds[var_name].values[:10]
                all_data[var_name] = pd.Series(values)
            elif len(dims) == 2:  # 二维变量
                # 对于二维变量，我们可以展平或者只取第一行/列
                # 这里我们选择只取第一行的前10个值
                values = ds[var_name].values[0, :10]
                all_data[var_name + "_first_row"] = pd.Series(values)
            else:  # 更高维度的变量
                # 对于更高维度的变量，我们可以展示其形状
                all_data[var_name + "_shape"] = pd.Series([str(ds[var_name].shape)])
        except Exception as e:
            print(f"处理变量 {var_name} 时出错: {e}")
            all_data[var_name + "_error"] = pd.Series(["Error processing this variable"])

    # 将所有数据合并到一个DataFrame中
    inspection_df = pd.DataFrame(all_data)
    inspection_df.to_csv(inspection_path)
    print(f"数据检查文件已保存到：{inspection_path}")