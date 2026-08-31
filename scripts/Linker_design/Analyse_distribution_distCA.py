#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import math
import re
from collections import defaultdict

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("错误: 未检测到 matplotlib 库，请运行: pip install matplotlib")
    sys.exit(1)

# =========================================================
# 路径与参数配置
# =========================================================
PDB_DIR = r"D:\Desktop\UL130_UL131A_linker"      # <- PDB 文件夹路径
OUTPUT_IMG = r"D:\Desktop\2linker_2in1_plot.png"  # <- 输出合成图路径

MIN_G = 5
MAX_G = 25

AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def calc_distance(p1, p2):
    """计算两点间欧几里得距离"""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def parse_single_chain_pdb(pdb_path):
    """
    解析单链 PDB 文件，提取每个残基的单字母代码与 CA 坐标
    """
    residues = []
    current_res_id = None
    current_res_name = None
    current_ca_coord = None

    with open(pdb_path, 'r', errors='ignore') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip()
                res_seq = line[22:26].strip()
                i_code = line[26].strip()
                res_id = (chain, res_seq, i_code)
                atom_name = line[12:16].strip()

                if res_id != current_res_id:
                    if current_res_id is not None and current_ca_coord is not None:
                        residues.append((AA_3TO1.get(current_res_name, 'X'), current_ca_coord))
                    current_res_id = res_id
                    current_res_name = line[17:20].strip()
                    current_ca_coord = None

                if atom_name == 'CA':
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        current_ca_coord = (x, y, z)
                    except ValueError:
                        continue

        if current_res_id is not None and current_ca_coord is not None:
            residues.append((AA_3TO1.get(current_res_name, 'X'), current_ca_coord))

    return residues


def main():
    if not os.path.exists(PDB_DIR):
        print(f"错误: 找不到文件夹 -> {PDB_DIR}")
        sys.exit(1)

    pdb_files = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb")))
    if not pdb_files:
        print(f"提示: 在 {PDB_DIR} 中没有找到任何 .pdb 文件！")
        sys.exit(1)

    print("=" * 80)
    print(f"扫描目录: {PDB_DIR} (共 {len(pdb_files)} 个 PDB 文件)")
    print("=" * 80)
    print(f"{'No.':<4} | {'PDB File Name':<32} | {'G-Len':<6} | {'Max CA-CA Dist (Å)':<20}")
    print("-" * 80)

    # 存储每个长度对应的所有结构的最大 CA-CA 间距
    len_to_max_ca_list = defaultdict(list)
    valid_count = 0

    for idx, pdb_file in enumerate(pdb_files, 1):
        filename = os.path.basename(pdb_file)
        residues = parse_single_chain_pdb(pdb_file)
        if not residues:
            continue

        seq = ''.join([r[0] for r in residues])
        ca_coords = [r[1] for r in residues]

        # 匹配 5 到 25 个连续 G
        matches = []
        for m in re.finditer(r'G+', seq):
            g_l = m.end() - m.start()
            if MIN_G <= g_l <= MAX_G:
                matches.append((m.start(), m.end(), g_l))

        if not matches:
            continue

        matches.sort(key=lambda x: x[2], reverse=True)
        start_idx, end_idx, g_len = matches[0]

        # 包含上下游各 1 个氨基酸
        ext_start = max(0, start_idx - 1)
        ext_end = min(len(residues), end_idx + 1)
        segment_ca = ca_coords[ext_start:ext_end]

        if len(segment_ca) >= 2:
            # 计算该片段内所有相邻 CA-CA 的距离
            ca_pair_dists = [
                calc_distance(segment_ca[i], segment_ca[i+1])
                for i in range(len(segment_ca) - 1)
            ]
            # 取最大距离
            max_ca_d = max(ca_pair_dists)
        else:
            max_ca_d = 0.0

        len_to_max_ca_list[g_len].append(max_ca_d)
        valid_count += 1
        print(f"{idx:<4} | {filename:<32} | {g_len:^6} | {max_ca_d:^20.3f}")

    print("-" * 80)
    print(f">> 统计完成: 共处理 {valid_count} 个符合条件的 PDB 结构。\n")

    if not len_to_max_ca_list:
        print("错误: 未匹配到符合条件的连续 G 结构！")
        sys.exit(0)

    # 完整覆盖从最小长度到最大长度的所有整数刻度（即使中间某长度为 0 也保留柱位）
    min_len = min(len_to_max_ca_list.keys())
    max_len = max(len_to_max_ca_list.keys())
    all_lengths = list(range(min_len, max_len + 1))

    # 指标 1: 各长度结构数量
    counts = [len(len_to_max_ca_list[l]) for l in all_lengths]

    # 指标 2: 各长度最大 CA 间距的平均值
    avg_max_ca = [
        sum(len_to_max_ca_list[l]) / len(len_to_max_ca_list[l]) if len(len_to_max_ca_list[l]) > 0 else None
        for l in all_lengths
    ]

    # =========================================================
    # 绘制 2 合 1 统计图
    # =========================================================
    fig_width = max(8, len(all_lengths) * 0.7)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width, 8.5), sharex=True, dpi=300)
    plt.subplots_adjust(hspace=0.2)

    # ------------ 子图 1: 各长度 PDB 数量直方图 (1 长度对应 1 根 Bar) ------------
    bars = ax1.bar(
        all_lengths,
        counts,
        width=0.6,
        color='#2b5c8f',
        edgecolor='black',
        linewidth=1.0,
        align='center'
    )
    max_c = max(counts) if counts else 1
    for x_val, count in zip(all_lengths, counts):
        if count > 0:
            ax1.text(x_val, count + max_c * 0.02, str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax1.set_ylabel('PDB Count', fontsize=11, fontweight='bold')
    ax1.set_title(f'1. PDB Count Distribution by Linker Length (Total: {valid_count})', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, max_c * 1.18)
    ax1.grid(axis='y', linestyle='--', alpha=0.4)

    # ------------ 子图 2: 最大 CA-CA 间距均值 ------------
    # 过滤掉 count 为 0 的点用于连线
    valid_x = [l for l, v in zip(all_lengths, avg_max_ca) if v is not None]
    valid_y = [v for v in avg_max_ca if v is not None]

    ax2.plot(valid_x, valid_y, color='#d62728', marker='o', linewidth=2.0, markersize=7, label='Mean Max CA-CA Dist')
    ax2.axhline(y=3.80, color='gray', linestyle='--', linewidth=1.2, alpha=0.8, label='Ideal Trans CA-CA (3.80 Å)')

    for x_val, y_val in zip(valid_x, valid_y):
        ax2.text(x_val, y_val + 0.015, f"{y_val:.3f}", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#a01818')

    ax2.set_xlabel('Consecutive Glycine (G) Length', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Max CA-CA Dist (Å)\n(G-run ± 1 AA)', fontsize=11, fontweight='bold')
    ax2.set_title('2. Mean of Max CA-CA Distance in [G-run + Flanking ±1 AA]', fontsize=12, fontweight='bold')
    ax2.set_xticks(all_lengths)
    ax2.set_xticklabels([str(x) for x in all_lengths], fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.legend(loc='upper right', frameon=True)

    y_min = min(valid_y + [3.80]) - 0.05
    y_max = max(valid_y + [3.80]) + 0.08
    ax2.set_ylim(y_min, y_max)

    plt.tight_layout()
    plt.savefig(OUTPUT_IMG)
    print(f"📊 统计图已成功生成并保存至: {OUTPUT_IMG}\n")


if __name__ == '__main__':
    main()