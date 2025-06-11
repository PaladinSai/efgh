import sgkit as sg
import matplotlib.pyplot as plt
from matplotlib.patches import  Patch
import seaborn as sns
import re
import pandas as pd
import numpy as np
import os
import logging
import matplotlib.patches as mpatches
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statannotations.Annotator import Annotator

# 设置日志记录格式
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

def get_gene_region(
    gff3_path, gene_name, feature_type="CDS", gene_id_key=("gene_name", "Name", "gene", "ID")
):
    """
    功能：根据GFF3文件和目标基因名，解析该基因的染色体、起止位置、CDS绝对/相对坐标。
    入参：
        gff3_path: GFF3注释文件路径
        gene_name: 目标基因名称
        feature_type: 特征类型，默认"CDS"
        gene_id_key: 用于匹配gene_name的属性键集合
    出参：
        chrom: 染色体名
        gene_start: 基因起始位置
        gene_end: 基因终止位置
        cds_coords_abs: 所有CDS绝对坐标列表
        cds_coords_rel: 所有CDS相对坐标列表
    """
    df = pd.read_csv(
        gff3_path,
        sep="\t",
        comment="#",
        header=None,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )
    gene_df = df[df["type"] == "gene"].copy()
    gene_df["gene_name_found"] = None
    for key in gene_id_key:
        m = gene_df["attributes"].str.extract(rf"{key}=([^;]+)")
        mask = (m[0] == gene_name)
        gene_df.loc[mask, "gene_name_found"] = m[0]
    matched = gene_df[gene_df["gene_name_found"] == gene_name]
    if matched.empty:
        raise ValueError(f"Gene '{gene_name}' not found by keys {gene_id_key} in GFF3.")

    gene_row = matched.iloc[0]
    chrom = gene_row["seqid"]
    gene_start, gene_end = int(gene_row["start"]), int(gene_row["end"])

    gene_id_match = None
    for key in gene_id_key:
        m = re.search(rf"{key}=([^;]+)", gene_row["attributes"])
        if m:
            gene_id_match = m.group(1)
            break
    if not gene_id_match:
        raise ValueError(f"Cannot extract gene_id for '{gene_name}'")

    tx_df = df[df["type"] == "transcript"].copy()
    tx_df["gid"] = None
    for key in gene_id_key:
        tx_df["gid"] = tx_df["gid"].combine_first(tx_df["attributes"].str.extract(rf"{key}=([^;]+)")[0])
    tx_df = tx_df[tx_df["gid"] == gene_id_match]
    if tx_df.empty:
        raise ValueError(f"No transcripts found for gene '{gene_name}'.")

    transcript_ids = tx_df["attributes"].str.extract(r"ID=([^;]+)")[0].unique().tolist()

    cds_df = df[df["type"] == feature_type].copy()
    cds_df["parent"] = cds_df["attributes"].str.extract(r"Parent=([^;]+)")
    cds_df = cds_df[cds_df["parent"].isin(transcript_ids)]
    if cds_df.empty:
        raise ValueError(f"No {feature_type} found for gene '{gene_name}'.")

    cds_coords_abs = [(int(start), int(end)) for start, end in zip(cds_df["start"], cds_df["end"])]
    cds_coords_rel = [(int(start) - gene_start, int(end) - gene_start) for start, end in zip(cds_df["start"], cds_df["end"])]

    return chrom, gene_start, gene_end, cds_coords_abs, cds_coords_rel

def load_genotype_data(ds, chrom, start, end, keep_absolute_position=True):
    """
    功能：在sgkit数据集ds中根据染色体和区段起止位置筛选变异位点，返回对应的子集。
    入参：
        ds: sgkit读取的原始Dataset
        chrom: 染色体名
        start: 区段起始坐标
        end: 区段终止坐标
        keep_absolute_position: 是否保留绝对坐标
    出参：
        ds_region: 筛选后的sgkit Dataset子集
    """
    chrom_query = chrom.lower().replace("chr", "")
    contig_ids = [str(cid).lower().replace("chr", "") for cid in ds['contig_id'].values]
    if chrom_query not in contig_ids:
        raise ValueError(f"Chromosome {chrom} not found in dataset contig_id {ds['contig_id'].values}")
    contig_index = contig_ids.index(chrom_query)

    variant_contig = ds["variant_contig"].values
    variant_position = ds["variant_position"].values
    variant_filter = (variant_contig == contig_index) & \
                     (variant_position >= start) & (variant_position <= end)
    if not np.any(variant_filter):
        raise ValueError(f"No variants found in region {chrom}:{start}-{end}")

    ds_region = ds.sel(variants=variant_filter)
    if keep_absolute_position:
        ds_region = ds_region.assign(region_variant_position=ds_region.variant_position)
    ds_region["variant_position"] = ds_region["variant_position"] - start
    return ds_region

def summarize_haplotypes(ds_region):
    """
    功能：统计ds_region中所有单倍型组合及其频率，为下游分析提供统一标签、频率、分组信息。
    入参：
        ds_region: sgkit Dataset，已限定分析区域
    出参：
        summary: dict，包含haplotype_to_label等所有下游分析所需变量
    """
    gt = ds_region["call_genotype"].values  # (variants, samples, ploidy)
    alleles = ds_region["variant_allele"].values  # (variants, alleles)
    sample_ids = ds_region["sample_id"].values.tolist()

    hap_counter = {}
    sample_main_hap = []
    sample_sub_hap = []
    for sample_idx in range(gt.shape[1]):
        h1 = tuple("DEL" if gt[v, sample_idx, 0] < 0 else alleles[v][gt[v, sample_idx, 0]] for v in range(gt.shape[0]))
        h2 = tuple("DEL" if gt[v, sample_idx, 1] < 0 else alleles[v][gt[v, sample_idx, 1]] for v in range(gt.shape[0]))
        hap_counter[h1] = hap_counter.get(h1, 0) + 1
        hap_counter[h2] = hap_counter.get(h2, 0) + 1
        sample_main_hap.append(h1)
        sample_sub_hap.append(h2)
    # 按频率降序、字典序排序
    haplotype_list = sorted(hap_counter, key=lambda h: (-hap_counter[h], h))
    haplotype_to_label = {h: f"H{str(i+1).zfill(3)}" for i, h in enumerate(haplotype_list)}
    label_to_haplotype = {v: k for k, v in haplotype_to_label.items()}
    haplotype_freq = {h: hap_counter[h] for h in haplotype_list}
    sample_main_label = [haplotype_to_label[h] for h in sample_main_hap]
    sample_sub_label = [haplotype_to_label[h] for h in sample_sub_hap]
    alleles_label = [f"{al[0]}/{al[1]}" for al in alleles]

    summary = dict(
        haplotype_to_label=haplotype_to_label,
        label_to_haplotype=label_to_haplotype,
        haplotype_freq=haplotype_freq,
        haplotype_list=haplotype_list,
        sample_main_label=sample_main_label,
        sample_sub_label=sample_sub_label,
        alleles_label=alleles_label,
        gt=gt,
        alleles=alleles,
        sample_ids=sample_ids,
        ds_region=ds_region
    )
    return summary

def calculate_ld_matrix(ds_region):
    """
    功能：计算区域内所有变异位点的LD矩阵，返回R²矩阵。
    入参：
        ds_region: sgkit Dataset，已限定分析区域
    出参：
        ld_matrix: numpy二维数组，LD值矩阵
    """
    ds_windowed = sg.window_by_variant(ds_region, size=ds_region.dims["variants"])
    ld_df = sg.ld_matrix(ds_windowed, dosage='call_dosage', threshold=0.0).compute()
    max_idx = max(ld_df.i.max(), ld_df.j.max()) + 1
    ld_matrix = np.zeros((max_idx, max_idx))
    for row in ld_df.itertuples():
        ld_matrix[row.i, row.j] = row.value
        ld_matrix[row.j, row.i] = row.value
    return ld_matrix

def plot_haplotypes(summary, output_path):
    """
    功能：绘制单倍型矩阵和各单倍型频率表（不再绘制柱状图），位置行用相对变异坐标。
    入参：
        summary: summarize_haplotypes返回的字典
        output_path: 输出图片路径
    出参：无（保存图片）
    """
    haplotype_list = summary["haplotype_list"]
    haplotype_to_label = summary["haplotype_to_label"]
    haplotype_freq = summary["haplotype_freq"]
    alleles_label = summary["alleles_label"]
    variant_positions = summary["ds_region"]["variant_position"].values  # 使用相对坐标
    label_list = [haplotype_to_label[h] for h in haplotype_list]
    freq_list = [haplotype_freq[h] for h in haplotype_list]
    hap_df = pd.DataFrame([list(h) for h in haplotype_list])
    hap_df["freq"] = freq_list
    hap_df.index = label_list

    # 位置行直接用相对坐标
    pos_row = pd.DataFrame(
        [list(variant_positions)],
        columns=hap_df.columns[:-1],
        index=["Position"]
    )
    allele_row = pd.DataFrame(
        [alleles_label],
        columns=hap_df.columns[:-1],
        index=["Alleles"]
    )
    full_df = pd.concat([pos_row, allele_row, hap_df])
    full_df["freq"] = ["", ""] + hap_df["freq"].astype(str).tolist()
    color_map = {
        "A": "#e41a1c", "T": "#984ea3",
        "C": "#4daf4a", "G": "#377eb8",
        "DEL": "#999999"
    }
    neutral_bg = "#f0f0f0"
    # 图像高度随单倍型数调整
    fig = plt.figure(figsize=(max(8, 1.5 * len(label_list)), 0.6 * (len(hap_df) + 6)))
    gs = fig.add_gridspec(1, 1)

    ax = fig.add_subplot(gs[0])
    ax.axis("off")
    table = ax.table(
        cellText=full_df.values,
        rowLabels=full_df.index,
        cellLoc='center',
        loc='center'
    )
    for key, cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(1.5)
        if key[0] < 2:
            cell.set_facecolor(neutral_bg)
        elif key[1] == len(full_df.columns) - 1:
            cell.set_facecolor(neutral_bg)
        else:
            val = full_df.iloc[key[0], key[1]]
            cell.set_facecolor(color_map.get(str(val).split('/')[0], "white"))
    legend_patches = [Patch(facecolor=color_map[nt], edgecolor='white', label=nt)
                      for nt in ["A", "T", "C", "G", "DEL"]]
    fig.legend(
        handles=legend_patches,
        title="Alleles",
        bbox_to_anchor=(1.01, 0.7),
        loc='upper left'
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_variant_positions(gene_name, gene_start, gene_end, cds_rel, variant_positions, alleles, output_path):
    """
    功能：优化后的变异位点分布图。主线为基因区段（细黑线），蓝色块为CDS（更矮），红色竖线为变异位点，图例不与主体重叠。
    入参：
        gene_name: 基因名称
        gene_start: 基因起始坐标（绝对坐标）
        gene_end: 基因终止坐标（绝对坐标）
        cds_rel: CDS相对坐标区间（相对gene_start的偏移）
        variant_positions: 变异位点相对坐标（相对于gene_start）
        alleles: 每个位点的等位基因字符串
        output_path: 输出图片路径
    出参：无（保存图片）
    """


    fig, ax = plt.subplots(figsize=(12, 2.3))
    region_len = gene_end - gene_start
    y_main = 0.5
    line_width = 4   # 主体线条变细
    cds_height = 0.20  # CDS高度降低

    # 主体线条（基因全长）更细
    ax.plot([0, region_len], [y_main, y_main], color='black', linewidth=line_width, solid_capstyle='round', zorder=1)

    # CDS绘制为更矮蓝色块（无边框）
    for cds_start, cds_end in cds_rel:
        ax.add_patch(
            mpatches.Rectangle((cds_start, y_main - cds_height/2), cds_end-cds_start, cds_height,
                               facecolor='#3182bd', edgecolor='none', linewidth=0, zorder=2)
        )
    # 变异位点竖线与标注
    for i, pos in enumerate(variant_positions):
        ax.plot([pos, pos], [y_main + cds_height/2, y_main + cds_height/2 + 0.22], color='red', linewidth=2, zorder=3)
        # 位置号标注在上方
        ax.text(pos, y_main + cds_height/2 + 0.25, f"{int(pos)}", ha='center', va='bottom', fontsize=10, color='red', rotation=60)
        # 等位基因信息标注在下方（可注释掉可隐藏）
        ax.text(pos, y_main - cds_height/2 - 0.13, alleles[i], ha='center', va='top', fontsize=10, color='blue')

    # 图例移到右上角且不叠加内容
    legend_handles = [
        mpatches.Patch(facecolor='black', edgecolor='none', label='Gene region', linewidth=2),
        mpatches.Patch(facecolor='#3182bd', edgecolor='none', label='CDS', linewidth=0),
        mpatches.Patch(facecolor='red', edgecolor='none', label='Variant', linewidth=2)
    ]
    ax.legend(handles=legend_handles, loc='upper right', frameon=False, bbox_to_anchor=(1.01, 1.01))

    # 其他美化
    ax.set_xlim(-region_len * 0.02, region_len * 1.02)
    ax.set_ylim(y_main - 0.55, y_main + 0.7)
    ax.set_yticks([])
    ax.set_xlabel(gene_name)
    ax.set_xticks([0, region_len])
    ax.set_xticklabels([f"{gene_start}", f"{gene_end}"])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_variant_lollipop(gene_name, gene_start, gene_end, cds_rel, variant_positions, alleles, output_path):
    """
    功能：用棒棒图风格展示变异位点分布。主线为基因区段，蓝色块为CDS（更矮），每个位点为红色lollipop（竖线+圆点），上方标注坐标，下方标注等位基因。
    入参：
        gene_name: 基因名称
        gene_start: 基因起始坐标（绝对坐标）
        gene_end: 基因终止坐标（绝对坐标）
        cds_rel: CDS相对坐标区间（相对gene_start的偏移，list of (start, end)）
        variant_positions: 变异位点相对坐标（相对于gene_start，list）
        alleles: 每个位点的等位基因字符串（list）
        output_path: 输出图片路径
    出参：无（保存图片）
    """
    fig, ax = plt.subplots(figsize=(12, 2.5))
    region_len = gene_end - gene_start
    y_base = 0.5
    lollipop_height = 0.32
    cds_height = 0.18

    # 主体线条（基因区段，细黑线）
    ax.hlines(y=y_base, xmin=0, xmax=region_len, color='black', linewidth=2, zorder=1)

    # CDS蓝色块，无边框
    for cds_start, cds_end in cds_rel:
        ax.add_patch(
            plt.Rectangle((cds_start, y_base - cds_height/2), cds_end - cds_start, cds_height,
                          color='#3182bd', ec='none', zorder=2)
        )

    # 棒棒图：每个位点一根红色竖线+顶端红色圆点
    for i, pos in enumerate(variant_positions):
        ax.vlines(x=pos, ymin=y_base, ymax=y_base + lollipop_height, color='red', linewidth=2, zorder=3)
        ax.plot(pos, y_base + lollipop_height, 'o', color='red', markersize=8, zorder=4)
        # 上方标注坐标
        ax.text(pos, y_base + lollipop_height + 0.06, f"{int(pos)}", ha='center', va='bottom', fontsize=10, color='red', rotation=60)
        # 下方标注等位基因
        ax.text(pos, y_base - cds_height/2 - 0.13, alleles[i], ha='center', va='top', fontsize=10, color='blue')

    # 图例移到右上角，不遮挡主体
    legend_handles = [
        mpatches.Patch(facecolor='black', edgecolor='none', label='Gene region', linewidth=2),
        mpatches.Patch(facecolor='#3182bd', edgecolor='none', label='CDS', linewidth=0),
        mpatches.Patch(facecolor='red', edgecolor='none', label='Variant (lollipop)', linewidth=2)
    ]
    ax.legend(handles=legend_handles, loc='upper right', frameon=False, bbox_to_anchor=(1.01, 1.01))

    # 其他美化
    ax.set_xlim(-region_len * 0.01, region_len * 1.01)
    ax.set_ylim(y_base - 0.5, y_base + 0.55)
    ax.set_yticks([])
    ax.set_xlabel(gene_name)
    ax.set_xticks([0, region_len])
    ax.set_xticklabels([f"{gene_start}", f"{gene_end}"])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_ld_heatmap(ld_matrix, positions, output_path):
    """
    方格型LD热图，positions为SNP坐标。
    """
    n = len(positions)
    mask = np.triu(np.ones_like(ld_matrix, dtype=bool), k=1)  # 遮掉上三角
    plt.figure(figsize=(max(6, n*0.5), max(4, n*0.5)))
    sns.heatmap(
        ld_matrix,
        mask=mask,
        cmap="coolwarm",
        vmin=0, vmax=1,
        xticklabels=positions,
        yticklabels=positions,
        square=True,
        cbar_kws={"label": r"$R^2$"}
    )
    plt.xlabel("SNP Position")
    plt.ylabel("SNP Position")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_haplotype_phenotype_boxplot(summary, phenotype, output_path):
    """
    绘制单倍型分组的表型分布箱线图，风格对标BMC Bioinformatics 2023 Fig 5，并美化横坐标显示、避免标签重叠。
    """
    # 组装数据
    box_df = pd.DataFrame({
        "sample": summary["sample_ids"],
        "haplotype": summary["sample_main_label"],
        "phenotype": phenotype
    })
    box_df = box_df.dropna(subset=["phenotype"])

    # 只保留真正有样本的单倍型标签，按样本数降序
    order = box_df["haplotype"].value_counts().index.tolist()
    palette = sns.color_palette("Set2", n_colors=len(order))

    # 剔除样本数小于2的分组
    group_counts = box_df.groupby("haplotype")["phenotype"].count()
    valid_groups = group_counts[group_counts > 1].index.tolist()
    box_df = box_df[box_df["haplotype"].isin(valid_groups)]
    order = [h for h in group_counts.index if h in valid_groups]

    # 统计各组样本数
    group_counts = box_df.groupby("haplotype")["phenotype"].count()

    # 生成新x轴标签（H001\nn=148）
    xticklabels = [f"{h}\nn={group_counts[h]}" for h in order]

    plt.figure(figsize=(max(7, 1.2*len(order)), 6))
    ax = sns.boxplot(
        x="haplotype", y="phenotype", data=box_df,
        palette=palette, order=order,
        width=0.5, linewidth=2, fliersize=0, boxprops=dict(alpha=0.8)
    )

    # # 叠加透明散点（jitter）
    # sns.stripplot(
    #     x="haplotype", y="phenotype", data=box_df,
    #     color="black", size=5, jitter=0.23, alpha=0.45, ax=ax, order=order, zorder=1
    # )
    #
    # # 添加均值点
    # group_means = box_df.groupby("haplotype")["phenotype"].mean()
    # for i, h in enumerate(order):
    #     if h in group_means.index:
    #         ax.scatter(i, group_means[h], color='red', s=55, marker='D', zorder=3, label="Mean" if i == 0 else "")

    # 横坐标美化，分两行显示且可选旋转
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(xticklabels, fontsize=12, rotation=0)  # 可改rotation=30, ha='right'

    # 显著性标记（Tukey HSD检验+星号）
    if len(order) > 1:
        pairs = []
        for i in range(len(order)):
            for j in range(i+1, len(order)):
                pairs.append((order[i], order[j]))
        annotator = Annotator(ax, pairs, data=box_df, x="haplotype", y="phenotype", order=order)
        # 计算p值
        tukey = pairwise_tukeyhsd(box_df["phenotype"], box_df["haplotype"])
        pval_dict = {}
        for row in tukey._results_table.data[1:]:
            a, b, p = row[0], row[1], row[4]
            pval_dict[(a, b)] = p
            pval_dict[(b, a)] = p  # 确保顺序无关
        # 取出与pairs顺序对应的p值
        pvals = [pval_dict.get((a, b), 1.0) for (a, b) in pairs]

        # === 打印每一对的p值和星号表示 ===
        # def pval_to_stars(p):
        #     if p < 0.001:
        #         return '***'
        #     elif p < 0.01:
        #         return '**'
        #     elif p < 0.05:
        #         return '*'
        #     else:
        #         return ''
        # for (a, b), p in zip(pairs, pvals):
        #     print(f"Comparison: {a} vs {b}, p-value={p:.4g}, significance='{pval_to_stars(p)}'")
        # 你的statannotations版本只支持pvals和num_comparisons
        annotator.set_pvalues_and_annotate(pvals)

    # 坐标轴美化
    ax.set_xlabel("Haplotype", fontsize=14)
    ax.set_ylabel("Phenotype", fontsize=14)
    ax.set_title(f"Phenotype distribution across haplotypes", fontsize=16, pad=10)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    sns.despine(ax=ax)

    # 图例优化
    handles, labels = ax.get_legend_handles_labels()
    if "Mean" in labels:
        ax.legend([handles[labels.index("Mean")]], ["Mean"], loc="upper right", frameon=False)
    else:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def run_haplotype(ds, config):
    """
    功能：主流程，完成单倍型区域解析、数据筛选、单倍型总结、LD和绘图。
    入参：
        ds: sgkit读取的原始Dataset
        config: 配置对象，包含基因名、gff路径、输出目录等信息
    出参：
        result_dict: 各图片文件路径的字典
    """
    logging.info(f"Loading region for gene: {config.haplotype.gene_name}")
    chrom, gene_start, gene_end, cds_abs, cds_rel = get_gene_region(config.input.gff_path, config.haplotype.gene_name)
    logging.info(f"chrom: {chrom}, start: {gene_start}, end: {gene_end}, cds: {cds_abs}")
    logging.info(f"Region: {chrom}:{gene_start}-{gene_end}")

    logging.info("Loading genotype data...")
    ds_region = load_genotype_data(ds, chrom, gene_start, gene_end)
    if ds_region.sizes["variants"] == 0:
        raise ValueError("No variants found in target gene region, unable to perform haplotype analysis.")
    start = 0
    end = gene_end - gene_start

    logging.info("Summarizing haplotypes...")
    summary = summarize_haplotypes(ds_region)

    result_path = os.path.join(config.output.outdir, "haplotype_frequency.png")
    logging.info("Plotting haplotype matrix and frequency...")
    plot_haplotypes(summary, result_path)

    position_plot_path = os.path.join(config.output.outdir, "haplotype_positions.png")
    logging.info("Generating variant position plot...")
    plot_variant_lollipop(
        config.haplotype.gene_name, start, end, cds_rel,
        ds_region["variant_position"].values.tolist(),
        summary["alleles_label"],
        position_plot_path
    )

    logging.info("Calculating LD matrix...")
    ld_matrix = calculate_ld_matrix(ds_region)
    ld_plot_path = os.path.join(config.output.outdir, "haplotype_ld_matrix.png")
    logging.info("Plotting LD matrix...")
    plot_ld_heatmap(ld_matrix, ds_region["variant_position"].values.astype(int), ld_plot_path)

    logging.info("Plotting haplotype-phenotype boxplot...")
    boxplot_path = os.path.join(config.output.outdir, "haplotype_boxplot.png")
    plot_haplotype_phenotype_boxplot(summary, ds_region[config.haplotype.phenotype].values, boxplot_path)

    return {
        "haplotype_plot": result_path,
        "position_plot": position_plot_path,
        "ld_matrix_plot": ld_plot_path,
        "phenotype_boxplot": boxplot_path
    }