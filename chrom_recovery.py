import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 读取数据
df = pd.read_csv("samples_snp_matrix.csv", index_col=0)

# 拆分染色体
df["Chrom"] = df.index.str.split(":").str[0]

# 数值类型转换
num_cols = df.columns.drop("Chrom")
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")

# 提取样品
parent_name = "554"

samples = [col for col in df.columns if col not in ["Chrom", parent_name]]

# 自动判断排序方式
if samples[0].startswith("homo"):
    samples = sorted(samples, key=lambda x: int(x.replace("homo", "")))
else:
    samples = sorted(samples, key=lambda x: int(x))

# 染色体排序（A01~D13）
chroms = sorted(df["Chrom"].unique(), key=lambda x: (x[0], int(x[1:])))

print("染色体数:", len(chroms))

# 计算恢复度
result = []

for chrom in chroms:
    sub = df[df["Chrom"] == chrom]
    
    parent_chr = sub[parent_name].to_numpy(dtype=float)
    
    for i, sample in enumerate(samples):
        arr = sub[sample].to_numpy(dtype=float)
        
        score = 1 - np.abs(arr - parent_chr) / 2
        
        result.append({
            "Sample": f"homo{i+1}",   # 这里直接编号
            "Chrom": chrom,
            "Recovery": np.nanmean(score)
        })

rec_df = pd.DataFrame(result)

# 转矩阵（样品 × 染色体）
heatmap_df = rec_df.pivot(index="Sample", columns="Chrom", values="Recovery")

# 👉 确保顺序 homo1 → homoN
heatmap_df = heatmap_df.reindex(
    [f"homo{i+1}" for i in range(len(samples))]
)

# 输出数据表
heatmap_df.to_csv("chrom_recovery_matrix.csv")
print("输出矩阵: chrom_recovery_matrix.csv")

# 作图（热图）
from matplotlib.colors import LinearSegmentedColormap

plt.figure(figsize=(12,8))

# 自定义绿色渐变：浅绿（低恢复）→ 深绿（高恢复）
cmap = LinearSegmentedColormap.from_list("green_grad", ["#d0f0c0", "#006400"])

im = plt.imshow(heatmap_df, aspect='auto', cmap=cmap)

# 颜色条
plt.colorbar(im, label="Recovery")

# 坐标轴
plt.xticks(range(len(heatmap_df.columns)), heatmap_df.columns, rotation=90)
plt.yticks(range(len(heatmap_df.index)), heatmap_df.index)

plt.xlabel("Chromosome")
plt.ylabel("Samples")
plt.title("Chromosome-wise Genome Recovery Heatmap")

plt.tight_layout()
plt.savefig("chrom_recovery_heatmap.png", dpi=300)
plt.show()

print("热图输出: chrom_recovery_heatmap.png")
