import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# 参数
window_size = 1_000_000
step_size = 500_000

input_file = "samples_snp_matrix.csv"

# 读取数据
df = pd.read_csv(input_file, index_col=0)

# 拆分染色体和位置
df["Chrom"] = df.index.str.split(":").str[0]
df["Pos"] = df.index.str.split(":").str[1].astype(int)

# 数值转换
num_cols = df.columns.drop(["Chrom", "Pos"])
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")

parent_name = "554"
samples = [c for c in df.columns if c not in ["Chrom", "Pos", parent_name]]

# 样品排序
if samples[0].startswith("homo"):
    samples = sorted(samples, key=lambda x: int(x.replace("homo", "")))
else:
    samples = sorted(samples, key=lambda x: int(x))

# 染色体排序
chroms = sorted(df["Chrom"].unique(), key=lambda x: (x[0], int(x[1:])))

print("样品数:", len(samples))
print("染色体数:", len(chroms))

# 滑窗计算
records = []

for chrom in chroms:
    sub = df[df["Chrom"] == chrom].sort_values("Pos")
    
    if sub.empty:
        continue
    
    max_pos = sub["Pos"].max()
    
    parent_chr = sub[parent_name].to_numpy(dtype=float)
    
    for start in range(0, max_pos, step_size):
        end = start + window_size
        
        mask = (sub["Pos"] >= start) & (sub["Pos"] < end)
        if mask.sum() == 0:
            continue
        
        parent_win = parent_chr[mask]
        
        for i, sample in enumerate(samples):
            arr = sub[sample].to_numpy(dtype=float)[mask]
            
            score = 1 - np.abs(arr - parent_win) / 2
            
            if np.all(np.isnan(score)):
                rec = np.nan
            else:
                rec = np.nanmean(score)
            
            records.append({
                "Sample": f"homo{i+1}",
                "Chrom": chrom,
                "Start": start,
                "End": end,
                "Recovery": rec
            })

win_df = pd.DataFrame(records)

print("窗口总数:", len(win_df))

# 输出滑窗数据
win_df.to_csv("sample_sliding.csv")

# 单个样品作图
single_dir = "sample_sliding_single"
os.makedirs(single_dir, exist_ok=True)

for sample in win_df["Sample"].unique():
    plt.figure(figsize=(14,4))
    
    sub = win_df[win_df["Sample"] == sample]
    
    offset = 0
    xticks = []
    xlabels = []
    
    for chrom in chroms:
        tmp = sub[sub["Chrom"] == chrom]
        if tmp.empty:
            continue
        
        x = (tmp["Start"] + tmp["End"]) / 2 + offset
        
        plt.plot(x, tmp["Recovery"], linewidth=1)
        
        xticks.append(x.mean())
        xlabels.append(chrom)
        
        offset = x.max()
    
    plt.xticks(xticks, xlabels, rotation=90)
    plt.ylim(0,1)
    plt.axhline(0.9, linestyle="--")
    
    plt.xlabel("Chromosome")
    plt.ylabel("Recovery")
    plt.title(f"{sample} Genome Recovery Curve")
    
    plt.tight_layout()
    plt.savefig(os.path.join(single_dir, f"{sample}_curve.png"), dpi=300)
    plt.close()

print("单样品曲线完成")

# 所有样品拼图
samples_list = sorted(win_df["Sample"].unique(), key=lambda x: int(x.replace("homo","")))

n = len(samples_list)
cols = 3
rows = (n + cols - 1) // cols

plt.figure(figsize=(18, 4*rows))

for i, sample in enumerate(samples_list):
    plt.subplot(rows, cols, i+1)
    
    sub = win_df[win_df["Sample"] == sample]
    
    offset = 0
    xticks = []
    xlabels = []
    
    for chrom in chroms:
        tmp = sub[sub["Chrom"] == chrom]
        if tmp.empty:
            continue
        
        x = (tmp["Start"] + tmp["End"]) / 2 + offset
        
        plt.plot(x, tmp["Recovery"], linewidth=0.8)
        
        xticks.append(x.mean())
        xlabels.append(chrom)
        
        offset = x.max()
    
    plt.xticks(xticks, xlabels, rotation=90, fontsize=6)
    plt.ylim(0,1)
    plt.axhline(0.9, linestyle="--")
    
    plt.title(sample, fontsize=10)

# 总标题
plt.suptitle("Genome-wide Recovery Curve", fontsize=18, y=0.995)

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig("sample_sliding_joint.png")
plt.show()

print("拼图完成")
