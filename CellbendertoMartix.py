"""
CellBender h5 → 10x mtx 格式转换器（修正版）

用法：
  python cellbender_dir2mtx.py /path/to/cellbender_out/CRRxxxx

效果：
  1) 自动在该目录里挑选"最后一次运行"的主文件：
       output_cellbender_retryN.h5（N 最大）或 output_cellbender.h5
     （自动忽略 *_posterior.h5 与 *_filtered.h5）
  2) 读取为 AnnData
  3) **转置矩阵**为 10x 标准格式（基因 × 细胞）
  4) 导出 10x mtx（matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz）
  5) 导出 metadata.csv
  输出目录固定为：<目录>/10x_mtx/
"""

import sys
import re
import gzip
import shutil
from pathlib import Path
from typing import Optional
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import scanpy as sc
from cellbender.remove_background.downstream import anndata_from_h5

RE_MAIN = re.compile(r"^output_cellbender(?:_retry(\d+))?\.h5$")


def pick_latest_main_file(dirpath: Path) -> Optional[Path]:
    """选出目录中最后一次运行的主结果 h5 文件。"""
    candidates = []
    for p in dirpath.glob("*.h5"):
        name = p.name
        if name.endswith("_posterior.h5") or name.endswith("_filtered.h5"):
            continue
        m = RE_MAIN.match(name)
        if m:
            retry = int(m.group(1)) if m.group(1) is not None else 0
            candidates.append((retry, p.resolve()))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def ensure_gene_symbols(adata) -> None:
    """确保导出 features 第二列有 gene_symbols；若无则用 var_names 填充。"""
    if "gene_symbols" not in adata.var.columns:
        adata.var["gene_symbols"] = adata.var_names


def gzip_if_needed(path: Path) -> None:
    """压缩为 .gz（若已是 .gz 则跳过）。"""
    if path.suffix == ".gz":
        return
    gz_path = Path(str(path) + ".gz")
    with open(path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    path.unlink()


def write_10x_mtx_gz(outdir: Path, adata) -> None:
    """
    手写 10x 输出（matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz）
    
    关键修改：
    1. CellBender 输出的 adata 是 (细胞 × 基因) 格式
    2. 10x 标准格式是 (基因 × 细胞)
    3. 需要转置矩阵并交换 barcodes 和 features
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # 获取原始维度信息
    print(f"[INFO] 原始 AnnData 维度: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # 转置矩阵：(细胞 × 基因) → (基因 × 细胞)
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    
    # 转置
    X_T = X.T.tocoo()
    print(f"[INFO] 转置后矩阵维度: {X_T.shape[0]} genes × {X_T.shape[1]} cells")

    # 1) 写 barcodes.tsv（细胞条形码）
    barcodes_path = outdir / "barcodes.tsv"
    with open(barcodes_path, "w") as f:
        for bc in adata.obs_names:
            f.write(f"{bc}\n")
    print(f"[INFO] 写入 {len(adata.obs_names)} 个细胞条形码到 barcodes.tsv")

    # 2) 写 features.tsv（基因信息：id, name, type）
    features_path = outdir / "features.tsv"
    var = adata.var.copy()
    
    # gene_ids: 优先 var['gene_ids']，否则用 var_names
    gene_ids = var["gene_ids"] if "gene_ids" in var.columns else adata.var_names
    # gene_symbols: 优先 var['gene_symbols']，否则用 var_names
    gene_symbols = var["gene_symbols"] if "gene_symbols" in var.columns else adata.var_names
    feature_type = "Gene Expression"
    
    with open(features_path, "w") as f:
        for gid, gsym in zip(gene_ids, gene_symbols):
            f.write(f"{gid}\t{gsym}\t{feature_type}\n")
    print(f"[INFO] 写入 {len(gene_ids)} 个基因信息到 features.tsv")

    # 3) 写 matrix.mtx（转置后的稀疏矩阵）
    mtx_path = outdir / "matrix.mtx"
    
    # 10x 期望整数计数；若为浮点，做安全处理（四舍五入转 int）
    if np.issubdtype(X_T.dtype, np.floating):
        print("[INFO] 检测到浮点数据，转换为整数...")
        data = np.rint(X_T.data).astype(np.int64, copy=False)
        X_T = sp.coo_matrix((data, (X_T.row, X_T.col)), shape=X_T.shape)
    else:
        # 确保是整型
        X_T = sp.coo_matrix(
            (X_T.data.astype(np.int64, copy=False), (X_T.row, X_T.col)), 
            shape=X_T.shape
        )

    # 写入 matrix.mtx
    mmwrite(mtx_path, X_T)
    print(f"[INFO] 写入矩阵到 matrix.mtx，维度: {X_T.shape[0]} × {X_T.shape[1]}")

    # 4) 统一压缩为 .gz
    print("[INFO] 压缩文件为 .gz 格式...")
    for p in [barcodes_path, features_path, mtx_path]:
        gzip_if_needed(p)
    print("[INFO] 压缩完成")


def export_metadata_csv(outdir: Path, adata) -> None:
    """导出 metadata.csv（保留条形码为索引列）。"""
    csv_path = outdir / "metadata.csv"
    adata.obs.to_csv(csv_path, index=True)
    print(f"[INFO] 导出 metadata.csv，共 {len(adata.obs)} 个细胞")


def main():
    if len(sys.argv) != 2:
        print("用法：python cellbender_dir2mtx.py /path/to/dir", file=sys.stderr)
        sys.exit(1)

    indir = Path(sys.argv[1])
    if not indir.is_dir():
        print(f"ERROR: 不是有效目录：{indir}", file=sys.stderr)
        sys.exit(1)

    h5file = pick_latest_main_file(indir)
    if h5file is None:
        print("ERROR: 未找到主结果文件（output_cellbender.h5 或 output_cellbender_retryN.h5）。", 
              file=sys.stderr)
        sys.exit(1)

    outdir = indir / "10x_mtx"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"[INFO] 输入目录: {indir}")
    print(f"[INFO] 选择 h5: {h5file}")
    print(f"[INFO] 输出目录: {outdir}")
    print("=" * 70)

    # 读取 AnnData
    print("[INFO] 读取 h5 文件...")
    adata = anndata_from_h5(str(h5file))
    adata.var_names_make_unique()
    ensure_gene_symbols(adata)
    
    print(f"[INFO] 加载完成: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"[INFO] 非零元素数: {adata.X.nnz if sp.issparse(adata.X) else np.count_nonzero(adata.X)}")

    # 导出 10x mtx（会自动转置）
    print("\n[INFO] 开始写出 10x mtx 格式（基因 × 细胞）...")
    write_10x_mtx_gz(outdir, adata)

    # 导出 metadata
    print("\n[INFO] 导出 metadata.csv...")
    export_metadata_csv(outdir, adata)

    print("\n" + "=" * 70)
    print("[DONE] 完成！生成文件：")
    print(f"  - {outdir / 'matrix.mtx.gz'}")
    print(f"  - {outdir / 'barcodes.tsv.gz'}")
    print(f"  - {outdir / 'features.tsv.gz'}")
    print(f"  - {outdir / 'metadata.csv'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
