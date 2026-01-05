#!/bin/bash

# 检查是否传入工作目录参数
if [ $# -ne 1 ]; then
    echo "用法: $0 <工作目录路径>"
    exit 1
fi

# 定义工作目录（从参数获取）
WORK_DIR="$1"

# 定义颜色变量用于输出美化
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # 无颜色

# 定义关键文件和目录路径（基于工作目录）
FILENAME_LIST="${WORK_DIR}/filename.txt"
RAW_DIR="${WORK_DIR}/raw"
# 临时文件用于存储原始列表，避免死循环
TEMP_LIST="${WORK_DIR}/.filename.tmp"

# 检查工作目录是否存在
if [ ! -d "$WORK_DIR" ]; then
    echo -e "${YELLOW}错误: 工作目录 $WORK_DIR 不存在${NC}"
    exit 1
fi

# 检查filename.txt是否存在
if [ ! -f "$FILENAME_LIST" ]; then
    echo -e "${YELLOW}错误: 未在工作目录找到 $FILENAME_LIST${NC}"
    exit 1
fi

# 检查raw目录是否存在
if [ ! -d "$RAW_DIR" ]; then
    echo -e "${YELLOW}错误: 未在工作目录找到 $RAW_DIR 目录${NC}"
    exit 1
fi

# 复制原始列表到临时文件，避免因写入原文件导致的死循环
cp "$FILENAME_LIST" "$TEMP_LIST"

# 从临时文件读取目录列表并处理
while IFS= read -r folder_name; do
    # 跳过空行
    if [ -z "$folder_name" ]; then
        continue
    fi

    # 构建完整的目标目录路径
    target_folder="${RAW_DIR}/${folder_name}"
    
    # 检查目录是否存在
    if [ ! -d "$target_folder" ]; then
        echo -e "${YELLOW}警告: 目录 $target_folder 不存在，跳过处理${NC}"
        continue
    fi
    
    echo -e "\n处理文件夹: ${GREEN}$target_folder${NC}"
    # 进入目标文件夹（若失败则跳过）
    cd "$target_folder" || {
        echo -e "${YELLOW}警告: 无法进入目录 $target_folder，跳过处理${NC}"
        continue
    }
    
    # 使用文件夹名称作为样本名
    sample_name="$folder_name"
    
    # 初始化lane计数器（R1和R2分别计数）
    declare -A lane_counter
    lane_counter["R1"]=1
    lane_counter["R2"]=1
    
    # 获取该文件夹下所有的fastq.gz文件，按名称排序
    files=($(find . -maxdepth 1 -type f -name "*.fastq.gz" | sort))
    
    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${YELLOW}警告: $target_folder 中没有找到fastq.gz文件${NC}"
        cd - >/dev/null || exit  # 返回上一级目录
        continue
    fi
    
    # 处理每个文件（直接重命名）
    for file in "${files[@]}"; do
        # 提取文件名
        basename_file=$(basename "$file")
        
        # 修正正则表达式，匹配实际文件名格式（包含中间的其他字符）
        if [[ "$basename_file" =~ R1.*\.fastq\.gz$ ]]; then
            read_type="R1"
            lane=${lane_counter["R1"]}
            lane_counter["R1"]=$((lane + 1))
        elif [[ "$basename_file" =~ R2.*\.fastq\.gz$ ]]; then
            read_type="R2"
            lane=${lane_counter["R2"]}
            lane_counter["R2"]=$((lane + 1))
        else
            echo -e "${YELLOW}警告: 无法识别文件类型: $basename_file${NC}"
            continue
        fi
        
        # 生成新文件名: 样本名_S1_L00x_Ry_001.fastq.gz
        new_name="${sample_name}_S1_L$(printf "%03d" $lane)_${read_type}_001.fastq.gz"
        
        # 直接在当前目录重命名文件
        echo "  重命名: $basename_file -> $new_name"
        mv "$file" "$new_name"
        
        # 检查重命名是否成功
        if [ $? -ne 0 ]; then
            echo -e "${YELLOW}警告: 重命名 $basename_file 失败${NC}"
        fi
    done
    
    echo -e "${GREEN}✓ $target_folder 处理完成${NC}"
    cd - >/dev/null || exit  # 返回上一级目录
    
    # 清除计数器变量
    unset lane_counter
done < "$TEMP_LIST"

# 清理临时文件
rm -f "$TEMP_LIST"

echo -e "\n${GREEN}所有文件夹处理完毕${NC}"
