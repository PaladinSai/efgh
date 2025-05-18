import seaborn as sns
import numpy as np
import sgkit as sg
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import matplotlib.pyplot as plt

def manhattan_plot(ds_lr, config):
    # df = ds_lr[["variant_contig_name", "variant_contig", "variant_position", "variant_linreg_p_value"]].to_dataframe()
    # df["variant_linreg_log_p_value"] = -np.log10(df["variant_linreg_p_value"])
    # df = df.astype({"variant_position": np.int64})
    #
    # running_pos = 0
    # cumulative_pos = []
    # for chrom, group_df in df.groupby("variant_contig"):
    #     cumulative_pos.append(group_df["variant_position"] + running_pos)
    #     running_pos += group_df["variant_position"].max()
    # df["cumulative_pos"] = pd.concat(cumulative_pos)
    # df["color"] = df["variant_contig"].apply(lambda x: "A" if x % 2 == 0 else "B")
    # g = sns.relplot(
    #     data=df,
    #     x="cumulative_pos",
    #     y="variant_linreg_log_p_value",
    #     hue="color group",
    #     palette=["blue", "orange"],
    #     linewidth=0,
    #     s=10,
    #     legend=None,
    #     aspect=3
    # )
    # g.ax.set_xlabel("Chromosome")
    # g.ax.set_xticks(df.groupby("variant_contig")["cumulative_pos"].median())
    # g.ax.set_xticklabels(df["variant_contig_name"].unique())
    #
    # output_dir = config.output.outdir
    # os.makedirs(output_dir, exist_ok=True)
    # output_path = os.path.join(output_dir, "manhattan_plot.png")
    # g.savefig(output_path)
    # plt.close()
    output_path = config.output.outdir
    os.makedirs(output_path, exist_ok=True)
    df = ds_lr[["variant_contig", "variant_position", "variant_linreg_p_value"]].to_dataframe()
    df["-log10p"] = -np.log10(df["variant_linreg_p_value"])
    df = df.astype({"variant_position": np.int64})
    running_pos = 0
    cumulative_pos = []
    for chrom, group_df in df.groupby("variant_contig"):
        cumulative_pos.append(group_df["variant_position"] + running_pos)
        running_pos += group_df["variant_position"].max()
    df["cumulative_pos"] = pd.concat(cumulative_pos)
    df["color"] = df["variant_contig"] % 2
    plt.figure(figsize=(18, 6))
    sns.scatterplot(data=df, x="cumulative_pos", y="-log10p", hue="color", palette=["blue", "orange"], legend=None,
                    s=10)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p-value)")
    plt.title("Manhattan plot (GWAS)")
    plt.tight_layout()
    out_file = os.path.join(output_path, "manhattan_plot.png")
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Manhattan plot saved to {out_file}")


def qq_plot(ds_lr, config):
    # p = ds_lr["variant_linreg_p_value"].squeeze().values
    # p.sort()
    # n = len(p)
    # expected_p = -np.log10(np.arange(1, n + 1) / n)
    # observed_p = -np.log10(p)
    # max_val = math.ceil(max(np.max(expected_p), np.max(observed_p)))
    #
    # df = pd.DataFrame({"Expected -log10(p)": expected_p, "Observed -log10(p)": observed_p})
    #
    # fig, ax = plt.subplots(figsize=(12, 12));
    # g = sns.scatterplot(data=df, x="Expected -log10(p)", y="Observed -log10(p)", ax=ax, linewidth=0)
    #
    # x_pred = np.linspace(0, max_val, 50)
    # sns.lineplot(x=x_pred, y=x_pred, ax=ax)
    #
    # g.set(xlim=(0, max_val), ylim=(0, max_val))
    #
    # output_dir = config.output.outdir
    # os.makedirs(output_dir, exist_ok=True)
    # output_path = os.path.join(output_dir, "qq_plot.png")
    # fig.savefig(output_path)
    # plt.close()
    output_path = config.output.outdir
    os.makedirs(output_path, exist_ok=True)
    p = ds_lr["variant_linreg_p_value"].squeeze().values
    p = np.sort(p)
    n = len(p)
    expected_p = -np.log10(np.arange(1, n + 1) / n)
    observed_p = -np.log10(p)
    plt.figure(figsize=(8, 8))
    plt.scatter(expected_p, observed_p, s=10)
    plt.plot([0, expected_p.max()], [0, expected_p.max()], color='red', lw=1)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("Q-Q plot (GWAS)")
    plt.tight_layout()
    out_file = os.path.join(output_path, "qq_plot.png")
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"QQ plot saved to {out_file}")
