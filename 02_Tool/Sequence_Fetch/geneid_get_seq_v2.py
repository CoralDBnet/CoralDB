#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @FileName      : geneid_get_seq.py
# @Time          : 2024-11-02 09:35:06
# @description   : 根据基因id，提取序列
"""
import argparse
from Bio import SeqIO


def read_seq_file(seq_file):
    seq_dict = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        gene_id = record.id
        sequence = str(record.seq)
        seq_dict[gene_id] = sequence
    return seq_dict


def write_sequences_from_gene_ids(gene_ids_file, seq_dict, output_file):
    with open(gene_ids_file, 'r') as gene_file, open(output_file, 'w') as result_file:
        for line in gene_file:
            gene_id = line.strip()
            if gene_id in seq_dict:
                result_file.write(f">{gene_id}\n{seq_dict[gene_id]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='根据基因ID提取序列')
    parser.add_argument('-g', '--genes', required=True, help='基因ID的输入文件')  # 改 - 基因id的输入文件
    parser.add_argument('-p', '--pep', default=None, help='输入的PEP序列文件')  # 改 - 输入的PEP序列文件
    parser.add_argument('-c', '--cds', default=None, help='输入的CDS序列文件')  # 改 - 输入的CDS序列文件
    parser.add_argument('-op', '--output_pep', default=None, help='PEP输出结果文件')  # 改 - PEP输出结果文件名字
    parser.add_argument('-oc', '--output_cds', default=None, help='CDS输出结果文件')  # 改 - CDS输出结果文件名字
    args = parser.parse_args()

    if args.pep is None and args.cds is None:
        parser.error("至少需要提供 -p/--pep 或 -c/--cds 其中一个输入文件")

    if args.pep:
        if args.output_pep is None:
            parser.error("提供 -p/--pep 时必须同时提供 -op/--output_pep")
        seq_dict = read_seq_file(args.pep)
        write_sequences_from_gene_ids(args.genes, seq_dict, args.output_pep)

    if args.cds:
        if args.output_cds is None:
            parser.error("提供 -c/--cds 时必须同时提供 -oc/--output_cds")
        seq_dict = read_seq_file(args.cds)
        write_sequences_from_gene_ids(args.genes, seq_dict, args.output_cds)

