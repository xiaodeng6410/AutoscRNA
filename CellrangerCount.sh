#!/bin/bash

# 检查参数是否正确传入
if [ $# -ne 2 ]; then
    echo "用法: $0 <CALL_DIR> <Species>"
    echo "示例: $0 /path/to/call_dir HSA"
    exit 1
fi

CALL_DIR=$1
Species=$2

# 创建输出目录
mkdir -p "${CALL_DIR}/cellranger_out/" || {
    echo "错误：无法创建输出目录 ${CALL_DIR}/cellranger_out/"
    exit 1
}
echo "已创建输出目录: ${CALL_DIR}/cellranger_out/"

# mkdir -p ${CALL_DIR}/cellranger_out/

# 根据物种选择参考基因组路径（严格匹配HSA/MMU）
if [ "$Species" = "HSA" ]; then
    # 人类参考基因组路径（请替换为实际路径）
    RFE_path="/home/dengys/ref/Homo_sapiens/GRCh38-2024-A"
    echo "物种为HSA（人类），使用参考基因组: $RFE_path"
elif [ "$Species" = "MMU" ]; then
    # 小鼠参考基因组路径
    RFE_path="/home/dengys/ref/Mus_musculus/Mus_musculus_113"
    echo "物种为MMU（小鼠），使用参考基因组: $RFE_path"
else
    # 未知物种时退出并提示
    echo "错误：未知物种 '$Species'，仅支持 HSA（人类）和 MMU（小鼠）"
    exit 1
fi

cat ${CALL_DIR}/filename.txt | while read id; do
  echo "开始处理样本: ${id}"
  # 进入输出文件夹
  cd ${CALL_DIR}/cellranger_out/
  # 运行 cellranger count，并把日志重定向到每个样本自己的目录
  cellranger count \
    --id=${id} \
    --create-bam=true \
    --fastqs=${CALL_DIR}/rawdata/${id} \
    --sample=${id} \
    --transcriptome="${RFE_path}" \
    > ${CALL_DIR}/cellranger_out/${id}_cellranger.log 2>&1
    
  echo "样本 ${id} 处理完成"
done

