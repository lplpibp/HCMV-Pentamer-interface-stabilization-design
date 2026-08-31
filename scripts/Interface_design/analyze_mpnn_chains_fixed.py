#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import re
import argparse
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

try:
    import logomaker
except ImportError:
    print("错误: 未安装 logomaker 库，请运行: pip install logomaker matplotlib pandas")
    sys.exit(1)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def read_fasta_with_chains(fasta_file):
    """
    读取 FASTA 文件，将每条包含 '/' 的序列拆分为各个链的独立列表。
    返回:
      - headers: 头部注释列表
      - parsed_seqs: 二维列表 [ [chain_seq_0, chain_seq_1, ...], ... ]
    """
    headers = []
    parsed_seqs = []
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    full_str = ''.join(current_seq)
                    # 按照 '/' 拆分成各链序列
                    parsed_seqs.append(full_str.split('/'))
                    current_seq = []
                headers.append(line)
            else:
                current_seq.append(line)
        if current_seq:
            full_str = ''.join(current_seq)
            parsed_seqs.append(full_str.split('/'))

    return headers, parsed_seqs


def parse_chain_order_from_header(header, default_chains):
    """从 FASTA header (如 chain_order='A B C D E' 或 'A, B, C') 提取链顺序"""
    match = re.search(r"chain_order=['\"]([^'\"]+)['\"]", header)
    if match:
        raw_order = match.group(1)
        # 兼容空格或逗号分隔
        chains = [c.strip() for c in re.split(r'[, ]+', raw_order) if c.strip()]
        return chains
    return default_chains


def get_pdb_dict_from_jsonl(jsonl_file):
    """
    从 parsed_pdbs.jsonl 中提取各链的原始 PDB 残基编号和野生型序列
    返回: { 'A': {'wt_seq': '...', 'pdb_res': [40, 41, ...]}, ... }
    """
    with open(jsonl_file, 'r') as f:
        data = json.loads(f.readline().strip())

    chain_dict = {}
    chain_keys = [k for k in data.keys() if k.startswith('seq_chain_')]
    for k in chain_keys:
        chain_id = k.split('_')[-1]
        seq = data.get(f'seq_chain_{chain_id}', '')
        pdb_res = data.get(f'pdb_chain_{chain_id}', [])
        chain_dict[chain_id] = {
            'wt_seq': seq,
            'pdb_res': pdb_res
        }
    return chain_dict


def main():
    parser = argparse.ArgumentParser(description="准确解析带'/'的多链 MPNN FASTA 并绘制各链设计位点 WebLogo")
    parser.add_argument("--fasta", "-f", required=True, help="MPNN 输出的 FASTA 文件 (.fa)")
    parser.add_argument("--jsonl", "-j", required=True, help="parsed_pdbs.jsonl 路径")
    parser.add_argument("--out", "-o", default="design_logo.png", help="输出图片文件名 (默认: design_logo.png)")
    args = parser.parse_args()

    # 1. 解析 FASTA 和 JSONL
    headers, all_parsed_seqs = read_fasta_with_chains(args.fasta)
    pdb_dict = get_pdb_dict_from_jsonl(args.jsonl)

    if len(all_parsed_seqs) < 2:
        print("错误: FASTA 文件中至少需要包含 WT 和 1 条生成序列！")
        sys.exit(1)

    # 提取链顺序 (例如 A, B, C, D, E)
    default_chain_order = sorted(pdb_dict.keys())
    chain_order = parse_chain_order_from_header(headers[0], default_chain_order)

    wt_chains = all_parsed_seqs[0]
    mpnn_designs = all_parsed_seqs[1:]
    num_designs = len(mpnn_designs)

    print("=" * 80)
    print(f"解析成功: 共包含 {len(chain_order)} 条链 (顺序: {', '.join(chain_order)})，{num_designs} 条设计序列")
    print("=" * 80)

    # 2. 逐链比对并收集设计位点
    design_info_list = []  # 存储所有突变位点的映射信息
    freq_matrix = []       # 存储频率数据 (用于绘图)
    x_labels = []          # 存储横坐标名称

    print(f"{'Chain':<7} | {'PDB ResNum':<12} | {'WT AA':<8} | {'Top Designed Amino Acids (Frequency)'}")
    print("-" * 80)

    for chain_idx, chain_id in enumerate(chain_order):
        if chain_idx >= len(wt_chains):
            continue

        wt_seq_chain = wt_chains[chain_idx]
        pdb_res_chain = pdb_dict.get(chain_id, {}).get('pdb_res', [])

        chain_has_mutation = False

        for pos_idx, wt_aa in enumerate(wt_seq_chain):
            # 提取所有设计模型在该链该位点的残基
            col_aas = [des[chain_idx][pos_idx] for des in mpnn_designs if len(des) > chain_idx and len(des[chain_idx]) > pos_idx]

            # 判断是否有突变（即不完全等同于 WT）
            if any(aa != wt_aa for aa in col_aas):
                chain_has_mutation = True
                
                # 获取真实 PDB 编号
                pdb_res_num = pdb_res_chain[pos_idx] if pos_idx < len(pdb_res_chain) else (pos_idx + 1)
                label = f"{chain_id}:{wt_aa}{pdb_res_num}"
                x_labels.append(label)

                # 统计仅来自设计序列的氨基酸频率
                counts = Counter(col_aas)
                probs = {aa: counts.get(aa, 0) / num_designs for aa in AMINO_ACIDS}
                freq_matrix.append(probs)

                top_aas = counts.most_common(4)
                top_str = ", ".join([f"{aa}: {cnt/num_designs:.1%}" for aa, cnt in top_aas])
                print(f"{chain_id:<7} | {pdb_res_num:<12} | {wt_aa:<8} | {top_str}")

        if chain_has_mutation:
            print("-" * 80)

    if not freq_matrix:
        print("未检测到任何突变设计位点。")
        return

    # 3. 绘制 WebLogo 图
    prob_df = pd.DataFrame(freq_matrix, index=range(len(x_labels)))
    fig_width = max(10, len(x_labels) * 0.65)
    fig, ax = plt.subplots(figsize=(fig_width, 4.8), dpi=300)

    logo = logomaker.Logo(
        prob_df,
        ax=ax,
        color_scheme='chemistry',
        vpad=0.02,
        width=0.8
    )

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True, linewidth=1.2)

    ax.set_ylim([0, 1.0])
    ax.set_ylabel("Frequency", fontsize=12, fontweight='bold')
    ax.set_xlabel("Chain : WT Residue & PDB Position", fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title("ProteinMPNN Designed Positions - Amino Acid Frequency", fontsize=14, fontweight='bold', pad=15)

    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9.5, fontweight='bold')
    ax.grid(axis='y', linestyle='--', alpha=0.35)

    plt.tight_layout()
    plt.savefig(args.out)
    print(f"\nWebLogo 图片已保存至: {args.out}\n")


if __name__ == "__main__":
    main()
