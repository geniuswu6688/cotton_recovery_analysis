import gzip
import numpy as np
import pandas as pd

vcf_file = "all_samples.joint.vcf.gz"  #VCF数据文件地址

records = []
samples = []

# =========================
# 1️⃣ 读取 VCF
# =========================
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
        
        # ✅ 只保留 SNP
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
        
        records.append((chrom, pos, ref, alt, gts))

print("原始SNP数:", len(records))

# =========================
# 2️⃣ 设置亲本
# =========================
parent_name = "554"
parent_index = samples.index(parent_name)

# =========================
# 3️⃣ 过滤并保存 SNP
# =========================
hq_snps = []

for chrom, pos, ref, alt, gts in records:
    arr = np.array(gts)
    
    # ❌ 去缺失
    if np.isnan(arr).any():
        continue
    
    # ❌ 去无多态
    if np.all(arr == arr[0]):
        continue
    
    parent_gt = arr[parent_index]
    
    hq_snps.append({
        "Chrom": chrom,
        "Pos": pos,
        "REF": ref,
        "ALT": alt,
        "Parent_GT": parent_gt
    })

print("高质量SNP数:", len(hq_snps))

# =========================
# 4️⃣ 输出文件
# =========================
df_snps = pd.DataFrame(hq_snps)

df_snps = df_snps.sort_values(["Chrom", "Pos"])

df_snps.to_csv("high_quality_snps.csv", index=False)

print("✅ 已输出 high_quality_snps.csv")
