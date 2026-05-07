import gzip
import numpy as np
import pandas as pd

vcf_file = "all_samples.joint.vcf.gz"

records = []
samples = []


# 读取 VCF
with gzip.open(vcf_file, "rt") as f:
    for line in f:
        if line.startswith("##"):
            continue
        
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            print("样品数:", len(samples))
            continue
        
        parts = line.strip().split("\t")
        
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        
        # 只保留 SNP
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        
        gts = []
        for sample_field in parts[9:]:
            gt = sample_field.split(":")[0]
            
            if gt == "./.":
                gts.append(np.nan)
            else:
                a, b = gt.replace("|", "/").split("/")
                gts.append(int(a) + int(b))
        
        records.append((chrom, pos, gts))

print("原始SNP数:", len(records))


# 过滤
filtered = []

for chrom, pos, gts in records:
    arr = np.array(gts)
    
    # 去缺失
    if np.isnan(arr).any():
        continue
    
    # 去无多态
    if np.all(arr == arr[0]):
        continue
    
    filtered.append((chrom, pos, arr))

print("过滤后SNP数:", len(filtered))


# 构建矩阵（位点为行）
site_names = []
matrix = []

for chrom, pos, arr in filtered:
    site_names.append(f"{chrom}:{pos}")
    matrix.append(arr)

matrix = np.array(matrix)

print("矩阵维度:", matrix.shape)


# DataFrame（样品为列）
df = pd.DataFrame(matrix, columns=samples, index=site_names)


# 亲本排第一列
parent_name = "554"

# 去掉亲本
others = [s for s in samples if s != parent_name]

# 按数字大小排序
others = sorted(others, key=lambda x: int(x))

# 重新排序列（亲本在第一列）
new_cols = [parent_name] + others
df = df[new_cols]

# 重命名列
new_names = {}

for i, sample in enumerate(others):
    new_names[sample] = f"homo{i+1}"

df = df.rename(columns=new_names)

# 输出
df.to_csv("samples_snp_matrix.csv")

print("输出完成")
