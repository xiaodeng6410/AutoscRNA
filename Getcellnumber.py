###Getcellnumber.py
#!/usr/bin/env python3
"""
从 CellRanger 的 web_summary.html 中提取细胞数目

使用方法:
    python Getcellnumber.py /path/to/sample_dir
    
示例:
    python Getcellnumber.py /data_result/dengys/scRNA/2511_DR_NIR/cellranger_out/C1/
    
输出:
    11812
"""

import os
import sys
import json


def extract_json_from_html(html: str, varname: str = "data") -> dict:
    """
    从 10x web_summary.html 里提取 `const data = {...}` 这段 JSON 并解析。
    """
    marker = f"const {varname} = "
    
    # 检查marker是否存在
    if marker not in html:
        raise ValueError(f"无法在 HTML 中找到 '{marker}'")
    
    start_pos = html.index(marker) + len(marker)
    s = html[start_pos:]
    
    in_string = False
    escape = False
    depth = 0
    json_start = None
    json_end = None
    
    for i, ch in enumerate(s):
        if escape:
            escape = False
            continue
        if ch == "\\":
            escape = True
            continue
        if ch == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if ch == "{":
            depth += 1
            if depth == 1:
                json_start = i
        elif ch == "}":
            depth -= 1
            if depth == 0 and json_start is not None:
                json_end = i + 1
                break
    
    if json_start is None or json_end is None:
        raise ValueError("无法在 HTML 中找到完整的 JSON 数据")
    
    json_str = s[json_start:json_end]
    return json.loads(json_str)


def get_estimated_number_of_cells(sample_dir: str) -> int:
    """
    从样本目录中提取细胞数目
    
    参数:
        sample_dir (str): 样本目录路径，例如 "/path/to/cellranger_out/C1"
        
    返回:
        int: 估计的细胞数目
    """
    # 构建 web_summary.html 的完整路径
    html_path = os.path.join(sample_dir, "outs", "web_summary.html")
    
    # 检查文件是否存在
    if not os.path.exists(html_path):
        raise FileNotFoundError(f"文件不存在: {html_path}")
    
    # 读取 web_summary.html
    with open(html_path, encoding="utf-8", errors="ignore") as f:
        html = f.read()
    
    # 提取 JSON 数据
    data = extract_json_from_html(html)
    
    # 10x 的 "Estimated Number of Cells" 对应这个路径：
    try:
        metric_str = data["summary"]["summary_tab"]["filtered_bcs_transcriptome_union"]["metric"]
        # 比如 "11,812" -> 11812
        return int(metric_str.replace(",", ""))
    except KeyError as e:
        raise ValueError(f"无法在 JSON 数据中找到细胞数字段: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python Getcellnumber.py <sample_directory>", file=sys.stderr)
        print("示例: python Getcellnumber.py /data_result/dengys/scRNA/2511_DR_NIR/cellranger_out/C1/", file=sys.stderr)
        sys.exit(1)
    
    sample_dir = sys.argv[1]
    
    try:
        n_cells = get_estimated_number_of_cells(sample_dir)
        print(n_cells)  # 只输出数字，方便脚本调用
        sys.exit(0)
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
