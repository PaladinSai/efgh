import numpy as np

def mask_to_numpy_in_chunks(mask, chunk_size=100000):
    """
    将大掩码（dask array 或 xarray.DataArray）分块转为numpy数组，避免内存溢出。
    返回拼接后的完整numpy布尔数组。
    """
    n = mask.shape[0]
    result = []
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        chunk = mask[start:end].compute() if hasattr(mask, "compute") else np.asarray(mask[start:end])
        # 保证是一维
        chunk = np.asarray(chunk).reshape(-1)
        result.append(chunk)
    return np.concatenate(result)