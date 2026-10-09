#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :convert_sequences.py
# @Time :2024/8/1 22:11
# @Last time: 2024/8/1 22:11
# python3 convert_sequences.py
# 将单个多序列比对格式文件转化为其他格式。
import argparse
from Bio import SeqIO


def convert_sequences(input_file, output_file, input_format, output_format):
    records = list(SeqIO.parse(input_file, input_format))
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, output_format)
    print(f'Converted {input_file} to {output_file}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='多序列比对格式转换')
    parser.add_argument('-i', '--input', required=True, help='输入的多序列比对文件')  # 改 - 输入的多序列比对文件路径
    parser.add_argument('-o', '--output', required=True, help='输出的结果文件')  # 改 - 输出的结果文件路径
    parser.add_argument('-if', '--input-format', required=True, help='输入文件格式（如 clustal, fasta, phylip 等）')  # 改 - 输入文件的多序列比对格式
    parser.add_argument('-of', '--output-format', required=True, help='输出文件格式（如 fasta, clustal, phylip 等）')  # 改 - 要改成的多序列比对文件格式
    args = parser.parse_args()

    convert_sequences(args.input, args.output, args.input_format, args.output_format)


# 本工具支持的 19 种格式（两个下拉菜单均包含以下选项）：
#   clustal            - ClustalW 多序列比对
#   embl               - EMBL 格式
#   fasta              - FASTA 格式
#   fasta-2line        - FASTA（序列不折行）
#   fastq              - FASTQ 质量分数
#   fastq-illumina     - Illumina FASTQ
#   fastq-sanger       - Sanger FASTQ
#   fastq-solexa       - Solexa FASTQ
#   gb                 - GenBank 格式
#   genbank            - GenBank 格式（别名）
#   nexus              - NEXUS 格式
#   phylip             - PHYLIP 交错格式
#   phylip-relaxed     - PHYLIP 宽松格式
#   phylip-sequential  - PHYLIP 顺序格式
#   pir                - PIR / NBRF 格式
#   qual               - 质量分数格式
#   seqxml             - SeqXML 格式
#   stockholm          - Stockholm 格式
#   tab                - 制表符分隔格式
