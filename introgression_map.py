import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os

# =========================
# 1️⃣ 读取滑窗数据
# =========================
df = pd.read_csv("sample_sliding.csv")

# =========================
# 2️⃣ 状态定义
# =========================
threshold = 0.9

df["State"] = np.where(df["Recovery"] >= threshold, "Parent", "Introgression")

# =========================
# 3️⃣ 排序
# =========================
samples = sorted(df["Sample"].unique(), key=lambda x: int(x.replace("homo","")))
chroms = sorted(df["Chrom"].unique(), key=lambda x: (x[0], int(x[1:])))

chrom_offset = {}
offset = 0

for chrom in chroms:
    sub = df[df["Chrom"] == chrom]
    max_pos = sub["End"].max()   # 用滑窗的End
    
    chrom_offset[chrom] = offset
    offset += max_pos

# =========================
# 4️⃣ 画图（区段图）
# =========================
plt.figure(figsize=(14, 6))

y_pos = {s: i for i, s in enumerate(samples)}
xticks = []
xlabels = []

for chrom in chroms:
    start = chrom_offset[chrom]
    end = start + df[df["Chrom"] == chrom]["End"].max()
    
    xticks.append((start + end) / 2)
    xlabels.append(chrom)

plt.xticks(xticks, xlabels, rotation=90)

for sample in samples:
    sub = df[df["Sample"] == sample]
    
    offset = 0
    
    for chrom in chroms:
        tmp = sub[sub["Chrom"] == chrom].sort_values("Start")
        if tmp.empty:
            continue
        
        for _, row in tmp.iterrows():
            offset = chrom_offset[chrom]

            start = row["Start"] + offset
            end = row["End"] + offset

            if row["State"] == "Parent":
                color = "#e6b800"
            else:
                color = "green"
            
            plt.hlines(y=y_pos[sample],
                       xmin=start,
                       xmax=end,
                       colors=color,
                       linewidth=12)
        
        mid = offset + tmp["End"].max() / 2
        xticks.append(mid)
        xlabels.append(chrom)
        offset += tmp["End"].max()
        
for sample in samples:
    y = y_pos[sample]
    
    for chrom in chroms:
        x = chrom_offset[chrom]
        
        plt.vlines(
            x=x,
            ymin=y - 0.4,
            ymax=y + 0.4,
            color="black",
            linewidth=0.6
        )

# =========================
# 5️⃣ 坐标轴
# =========================
plt.yticks(range(len(samples)), samples, fontsize=10)
plt.gca().invert_yaxis()
plt.xticks(xticks, xlabels, fontsize=10)
plt.xlabel("Genome Position")
plt.ylabel("Samples")
plt.title("Genome-wide Introgression Map")

parent_line = mlines.Line2D([], [], color='#e6b800', linewidth=4, label='Recurrent parent')
intro_line = mlines.Line2D([], [], color='green', linewidth=4, label='Introgression')

plt.legend(
    handles=[parent_line, intro_line],
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0
)

# =========================
# 6️⃣ 保存
# =========================
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("introgression_map.png", dpi=600)
plt.show()

print("✅ 完成：introgression_map.png")