#!/bin/bash
set -e

# 1. 基础配置
MPNN_PATH=~/softwares/ProteinMPNN-main
PDB_PATH="./WT.pdb"
OUTPUT_DIR="./mpnn_output"

mkdir -p ${OUTPUT_DIR}/inputs
mkdir -p ${OUTPUT_DIR}/outputs

# 2. 复制 PDB 并解析为张量
cp ${PDB_PATH} ${OUTPUT_DIR}/inputs/

python ${MPNN_PATH}/helper_scripts/parse_multiple_chains.py \
    --input_path=${OUTPUT_DIR}/inputs \
    --output_path=${OUTPUT_DIR}/parsed_pdbs.jsonl

python ${MPNN_PATH}/helper_scripts/make_fixed_positions_dict.py \
    --input_path=${OUTPUT_DIR}/parsed_pdbs.jsonl \
    --output_path=${OUTPUT_DIR}/fixed_positions.jsonl \
    --chain_list="B C D" \
    --position_list="105 108 110 116 118 123 125, 111 112 113 115 116 117, 21 22 26 53 54 55" \
    --specify_non_fixed
# 3. 按 PDB 真实编号生成固定位点字典（仅设计指定位点，A/E 自动全固定）
#python make_fixed_dict_by_pdb_res.py \
   # --jsonl_path=${OUTPUT_DIR}/parsed_pdbs.jsonl \
   # --output_path=${OUTPUT_DIR}/fixed_positions.jsonl \
   # --design_res="B:141,144,146,152,154,159,161 C:139,140,141,143,144,145 D:65,66,70,97,98,99"
   # --design_res="B:105,108,110,116,118,123,125 C:111,112,113,115,116,117 D:21,22,26,53,54,55"

# 4. 运行 ProteinMPNN 设计序列（生成 50 条序列作为示例）
python ${MPNN_PATH}/protein_mpnn_run.py \
    --jsonl_path=${OUTPUT_DIR}/parsed_pdbs.jsonl \
    --fixed_positions_jsonl=${OUTPUT_DIR}/fixed_positions.jsonl \
    --out_folder=${OUTPUT_DIR}/outputs \
    --num_seq_per_target=100 \
    --sampling_temp="0.2" \
    --seed=37 \
    --batch_size=1

# 5. 分析设计位点并在终端打印 + 生成 WebLogo 图片
python analyze_mpnn_chains_fixed.py \
    --fasta ${OUTPUT_DIR}/outputs/seqs/WT.fa \
    --jsonl ${OUTPUT_DIR}/parsed_pdbs.jsonl \
    --out design_logo.png

echo "全部运行完毕！WebLogo 已保存在当前目录下的 design_logo.png"
