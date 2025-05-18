import os
import pandas as pd

def out_test_csv(ds, config, path):
    # 输出所有变量的前10行数据到CSV文件
    # Output the first 10 rows of all variables in the dataset to a CSV file
    print("Exporting the first 10 rows of all variables in the dataset for inspection...")
    inspection_path = os.path.join(config.output.outdir, path)

    # 创建一个空的DataFrame来存储所有变量
    # Create an empty DataFrame to store all variables
    all_data = {}

    # 循环遍历所有数据变量
    # Iterate over all data variables
    for var_name in ds.data_vars:
        try:
            # 获取变量的维度信息
            # Get the dimension information of the variable
            dims = ds[var_name].dims
            print(f"Variable {var_name} dimensions: {dims}")

            # 对于不同维度的变量，我们可能需要不同的处理方式
            # Different handling for variables with different dimensions
            if len(dims) == 1:  # 一维变量 / 1D variable
                # 取前10个值 / Take the first 10 values
                values = ds[var_name].values[:10]
                all_data[var_name] = pd.Series(values)
            elif len(dims) == 2:  # 二维变量 / 2D variable
                # 只取第一行的前10个值 / Take the first 10 values of the first row
                values = ds[var_name].values[0, :10]
                all_data[var_name + "_first_row"] = pd.Series(values)
            else:  # 更高维度的变量 / higher-dimensional variable
                # 展示其形状 / Show its shape
                all_data[var_name + "_shape"] = pd.Series([str(ds[var_name].shape)])
        except Exception as e:
            print(f"Error processing variable {var_name}: {e}")
            all_data[var_name + "_error"] = pd.Series(["Error processing this variable"])

    # 将所有数据合并到一个DataFrame中
    # Combine all data into a single DataFrame
    inspection_df = pd.DataFrame(all_data)
    inspection_df.to_csv(inspection_path)
    print(f"Inspection file saved to: {inspection_path}")
