import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 读取矩阵
df = pd.read_csv("samples_snp_matrix.csv", index_col=0)

# 提取亲本
parent = df["554"].values
df_samples = df.drop(columns=["554"])

# 计算恢复度
recovery = []

for sample in df_samples.columns:
    arr = df_samples[sample].values
    score = 1 - np.abs(arr - parent) / 2
    recovery.append(score.mean())

recovery_df = pd.DataFrame({
    "Sample": df_samples.columns,  # 保留原样品名
    "Recovery": recovery
})

# 排序（从高到低）
recovery_df = recovery_df.sort_values("Recovery", ascending=False).reset_index(drop=True)

# 输出表格（去掉重编号列，只保留原样品名 + Recovery）
recovery_df.to_csv("average_recovery_rank.csv", index=False)
print("表格输出: average_recovery_rank.csv")

# 设置颜色（最低3个黄色，其余绿色）
colors = ["green"] * len(recovery_df)
for i in range(len(recovery_df)-3, len(recovery_df)):
    colors[i] = "gold"

# 作图
plt.figure(figsize=(12,5))
plt.bar(recovery_df["Sample"], recovery_df["Recovery"], color=colors)
plt.xticks(rotation=45)
plt.ylabel("Genome Recovery Rate")
plt.xlabel("Samples")
plt.title("Genome-wide Recovery Rate")
plt.axhline(0.9, linestyle="--")  # 可选阈值线
plt.tight_layout()
plt.savefig("average_recovery_rank.png", dpi=300)
plt.show()
print("图已输出: average_recovery_rank.png")
