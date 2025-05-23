import seaborn as sns
import numpy as np
import math
import pandas as pd
import os
import matplotlib.pyplot as plt
import logging

def manhattan_plot(ds_lr, config, trait, trait_idx):
    """
    绘制Manhattan图，用于展示GWAS结果。
    Draw Manhattan plot for GWAS results.
    """
    try:
        # 读取配置 / Read config
        man_cfg = getattr(config.gwas, "manhattan", None)
        chrom_colors = getattr(man_cfg, "chrom_colors", ["#1f77b4", "#ff7f0e"]) if man_cfg else ["#1f77b4", "#ff7f0e"]
        threshold_lines = getattr(man_cfg, "threshold_lines",
                                  [{"value": 5, "style": "solid", "color": "#000000"}]) if man_cfg else [
            {"value": 5, "style": "solid", "color": "#000000"}]
        point_colors = getattr(man_cfg, "point_colors", {}) if man_cfg else {}
        color_below = getattr(point_colors, "below", chrom_colors[0])
        color_above = getattr(point_colors, "above", "#d62728")
        color_between = getattr(point_colors, "between", "#2ca02c")
        # 数据准备 / Data preparation
        df = pd.DataFrame({
            "variant_contig_name": ds_lr["variant_contig_name"].values,
            "variant_contig": ds_lr["variant_contig"].values,
            "variant_position": ds_lr["variant_position"].values,
            "variant_linreg_p_value": ds_lr["variant_linreg_p_value"][:, trait_idx].values
        })
        df["variant_linreg_log_p_value"] = -np.log10(df["variant_linreg_p_value"])
        df = df.astype({"variant_position": np.int64})

        # 染色体累计坐标 / Chromosome cumulative position
        running_pos = 0
        cumulative_pos = []
        for chrom, group_df in df.groupby("variant_contig"):
            cumulative_pos.append(group_df["variant_position"] + running_pos)
            running_pos += group_df["variant_position"].max()
        df["cumulative_pos"] = pd.concat(cumulative_pos)
        df["chrom_color"] = df["variant_contig"] % 2

        # 阈值线处理 / Threshold line handling
        thresholds = sorted([line["value"] for line in threshold_lines])
        # 分类 / Categorization
        if len(thresholds) == 1:
            df["point_cat"] = np.where(df["variant_linreg_log_p_value"] >= thresholds[0], "above", "below")
        elif len(thresholds) == 2:
            t1, t2 = thresholds
            df["point_cat"] = np.select(
                [df["variant_linreg_log_p_value"] >= t2,
                 df["variant_linreg_log_p_value"] >= t1],
                ["above", "between"],
                default="below"
            )
        else:
            df["point_cat"] = "below"

        # 绘图 / Plotting
        fig, ax = plt.subplots(figsize=(16, 6))
        for chrom, color_idx in zip(df["variant_contig"].unique(), [0, 1] * len(df["variant_contig"].unique())):
            for cat, color in [("below", color_below), ("above", color_above), ("between", color_between)]:
                sub = df[(df["variant_contig"] == chrom) & (df["point_cat"] == cat)]
                if not sub.empty:
                    ax.scatter(sub["cumulative_pos"], sub["variant_linreg_log_p_value"],
                                c=color if cat != "below" else chrom_colors[color_idx % len(chrom_colors)],
                                s=10, label=f"{chrom}-{cat}" if cat != "below" else None, alpha=0.8)

        # 阈值线 / Threshold lines
        for i, line in enumerate(threshold_lines):
            y = line.get("value", 5)
            style = line.get("style", "solid")
            color = line.get("color", "#d62728")
            ax.axhline(y, color=color, linestyle="-" if style == "solid" else "--", linewidth=1.5, zorder=0)

        # 坐标轴与标签 / Axis and labels
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("-log10(p-value)")
        ax.set_xticks(df.groupby("variant_contig")["cumulative_pos"].median())
        ax.set_xticklabels(df["variant_contig_name"].unique())
        ax.set_title(f"Manhattan plot for {trait}")

        output_dir = config.output.outdir
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"manhattan_plot_{trait}.png")
        fig.savefig(output_path)
        plt.close()
        logging.info(f"Manhattan plot saved to {output_path}")
    except Exception:
        logging.error("Failed to generate Manhattan plot.")
        # 不抛出异常，继续流程 / Do not raise, just log

def manhattan_plot_chunked(ds_lr, config, trait, trait_idx, chunk_size=10000):
    """
    低内存分块流式绘制Manhattan图，适合xarray+zarr大数据。不刷Dask run_spec警告。
    """
    try:
        man_cfg = getattr(config.gwas, "manhattan", None)
        chrom_colors = getattr(man_cfg, "chrom_colors", ["#1f77b4", "#ff7f0e"]) if man_cfg else ["#1f77b4", "#ff7f0e"]
        threshold_lines = getattr(man_cfg, "threshold_lines",
                                  [{"value": 5, "style": "solid", "color": "#000000"}]) if man_cfg else [
            {"value": 5, "style": "solid", "color": "#000000"}]
        point_colors = getattr(man_cfg, "point_colors", {}) if man_cfg else {}
        color_above = getattr(point_colors, "above", "#d62728")
        color_between = getattr(point_colors, "between", "#2ca02c")

        # 1. 获取基础信息（直接dask array，不要xarray DataArray）
        contig = ds_lr["variant_contig"].data
        pos = ds_lr["variant_position"].data
        contig_name_arr = ds_lr["variant_contig_name"].data
        pval = ds_lr["variant_linreg_p_value"][:, trait_idx].data

        n = ds_lr.sizes.get("variants") or ds_lr.dims.get("variants")

        # 2. 统计unique contig及offsets和tick label
        offsets = {}
        tick_pos = []
        tick_label = []
        running_offset = 0

        # 2.1 找出unique_contigs及每个contig的第一个索引
        contig_first_idx = {}
        unique_contigs = []
        for start in range(0, n, chunk_size):
            idx = slice(start, min(start + chunk_size, n))
            contig_chunk = contig[idx].compute()
            for i, c in enumerate(contig_chunk):
                if c not in contig_first_idx:
                    contig_first_idx[c] = start + i
                    unique_contigs.append(c)
        unique_contigs = sorted(unique_contigs)

        # 2.2 统计offsets和tick
        for c in unique_contigs:
            max_pos = 0
            min_pos = None
            for start in range(0, n, chunk_size):
                idx = slice(start, min(start + chunk_size, n))
                contig_chunk = contig[idx].compute()
                pos_chunk = pos[idx].compute()
                mask = contig_chunk == c
                if np.any(mask):
                    max_pos = max(max_pos, pos_chunk[mask].max())
                    _min = pos_chunk[mask].min()
                    min_pos = _min if min_pos is None else min(min_pos, _min)
            offsets[c] = running_offset
            tick_pos.append((min_pos + max_pos) // 2 + running_offset)
            idx0 = contig_first_idx[c]
            tick_label.append(str(contig_name_arr[idx0].compute().item()))
            running_offset += max_pos

        # 3. 分块绘图
        fig, ax = plt.subplots(figsize=(16, 6))
        for start in range(0, n, chunk_size):
            idx = slice(start, min(start + chunk_size, n))
            contig_chunk = contig[idx].compute()
            pos_chunk = pos[idx].compute()
            pval_chunk = pval[idx].compute()
            log_pval_chunk = -np.log10(pval_chunk)
            cumulative_pos_chunk = np.array([pos_chunk[j] + offsets[contig_chunk[j]] for j in range(len(pos_chunk))])
            for i, c in enumerate(unique_contigs):
                ind = np.where(contig_chunk == c)[0]
                if ind.size == 0:
                    continue
                thresholds = sorted([line["value"] for line in threshold_lines])
                if len(thresholds) == 1:
                    cat = np.where(log_pval_chunk[ind] >= thresholds[0], "above", "below")
                elif len(thresholds) == 2:
                    t1, t2 = thresholds
                    cat = np.full(ind.shape, "below", dtype=object)
                    cat[log_pval_chunk[ind] >= t1] = "between"
                    cat[log_pval_chunk[ind] >= t2] = "above"
                else:
                    cat = np.full(ind.shape, "below", dtype=object)
                for cat_name, color in [("below", chrom_colors[i % len(chrom_colors)]),
                                        ("above", color_above),
                                        ("between", color_between)]:
                    sub_idx = ind[cat == cat_name]
                    if sub_idx.size > 0:
                        ax.scatter(cumulative_pos_chunk[sub_idx], log_pval_chunk[sub_idx],
                                   c=color, s=10, label=None, alpha=0.8)

        # 4. 阈值线
        for line in threshold_lines:
            y = line.get("value", 5)
            style = line.get("style", "solid")
            color = line.get("color", "#d62728")
            ax.axhline(y, color=color, linestyle="-" if style == "solid" else "--", linewidth=1.5, zorder=0)

        # 5. 轴标签
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_label)
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title(f"Manhattan plot for {trait}")

        output_dir = config.output.outdir
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"manhattan_plot_{trait}.png")
        fig.savefig(output_path)
        plt.close()
        logging.info(f"Manhattan plot saved to {output_path}")
    except Exception as e:
        logging.error(f"Failed to generate Manhattan plot (chunked): {e}")

def qq_plot(ds_lr, config, trait, trait_idx):
    """
    绘制QQ图，用于展示GWAS结果的p值分布。
    Draw QQ plot for GWAS p-value distribution.
    """
    try:
        p = ds_lr["variant_linreg_p_value"][:, trait_idx].squeeze().values
        p.sort()
        n = len(p)
        expected_p = -np.log10(np.arange(1, n + 1) / n)
        observed_p = -np.log10(p)
        max_val = math.ceil(max(np.max(expected_p), np.max(observed_p)))

        df = pd.DataFrame({"Expected -log10(p)": expected_p, "Observed -log10(p)": observed_p})

        fig, ax = plt.subplots(figsize=(12, 12))
        g = sns.scatterplot(data=df, x="Expected -log10(p)", y="Observed -log10(p)", ax=ax, linewidth=0)

        x_pred = np.linspace(0, max_val, 50)
        sns.lineplot(x=x_pred, y=x_pred, ax=ax)

        g.set(xlim=(0, max_val), ylim=(0, max_val))

        output_dir = config.output.outdir
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"qq_plot_{trait}.png")
        fig.savefig(output_path)
        plt.close()
        logging.info(f"QQ plot saved to {output_path}")
    except Exception:
        logging.error("Failed to generate QQ plot.")
        # 不抛出异常，继续流程 / Do not raise, just log

def plot_manhattan_and_qq_mp(ds_lr, config, trait, trait_idx):
    """
    多进程绘制Manhattan图和QQ图。
    Draw Manhattan and QQ plots for GWAS results in a subprocess.
    """
    try:
        # 子进程可单独设置日志格式，避免日志混乱
        logging.basicConfig(
            level=logging.INFO,
            format=f'[%(process)d][%(asctime)s] %(levelname)s: %(message)s'
        )
        logging.info(f"Start plotting for trait: {trait} (index={trait_idx})")
        manhattan_plot_chunked(ds_lr, config, trait, trait_idx, chunk_size=10000)
        qq_plot(ds_lr, config, trait, trait_idx)
        logging.info(f"Plots for trait: {trait} finished.")
    except Exception as e:
        logging.error(f"Plotting failed for trait {trait}: {e}")