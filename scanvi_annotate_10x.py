import tempfile
import numpy as np
import scanpy as sc
import sys
import scvi
import seaborn as sns
import torch
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import os

if len(sys.argv) < 2:
    print("❗请输入10X输入目录\n示例：python run_scanvi.py /path/to/10x/")
    sys.exit(1)

path = sys.argv[1]###
print(f"📂 读取10X路径: {path}")
###导入ref
adata_ref =sc.read("/data_result/dengys/scRNA/AutoScRNA/Brain_scRNA_ref/adata_ref_aligned_gencode.h5ad")
###导入共有基因
with open("/data_result/dengys/scRNA/AutoScRNA/Brain_scRNA_ref/common_genes.txt") as f:
    common_genes = [line.strip() for line in f]
###导入model
scanvi_model = scvi.model.SCANVI.load("/data_result/dengys/scRNA/AutoScRNA/Brain_scRNA_ref/scanvi_model_gencode", adata=adata_ref)
###导入



adata_query = sc.read_10x_mtx(path)
adata_query.var_names_make_unique()
adata_query = adata_query[:, common_genes].copy()

adata_query.obs['library_label'] = 'query_batch'
adata_query.obs['library_label'] = adata_query.obs['library_label'].astype('category')
adata_query.obs['cell_name'] = 'Unknown'
adata_query.obs['cell_name'] = adata_query.obs['cell_name'].astype('category')

# 注释查询数据
scanvi_query = scvi.model.SCANVI.load_query_data(
    adata_query,
    scanvi_model
)
# ---- 训练 ----
print("🚀 开始SCANVI迁移学习 ...")
scanvi_query.train(
    max_epochs=1000,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10,
    early_stopping = True,
    early_stopping_patience = 50,
    early_stopping_min_delta=0.01)

# 预测
print("🔍 预测细胞类型 ...")
adata_query.obs["predicted_celltype"] = scanvi_query.predict()
adata_query.obs["prediction_confidence"] = scanvi_query.predict(soft=True).max(axis=1)
# 统计每个 predicted_celltype 的细胞数目
print("\n📊 细胞注释统计：")
cell_type_counts = adata_query.obs["predicted_celltype"].value_counts()
print(cell_type_counts)

save_path = os.path.join(path, "meta_cellperdict.csv")
adata_query.obs.to_csv(save_path)
print(f"\n✅ 细胞注释完成，meta_cellperdict.csv 已保存到: {save_path}")

print(f"\n✅ 细胞注释完成，meta_cellperdict.csv 已保存到: {save_path}")
print("🎉 All done!")
