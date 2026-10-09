#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : translate.py
# @Date          : 2026-07-21
# @Description   : DNA序列翻译为蛋白质
"""
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def translate_dna(input_fasta, output_fasta):
    records = []
    for r in SeqIO.parse(input_fasta, 'fasta'):
        protein_seq = r.seq.translate(to_stop=False)
        records.append(SeqRecord(protein_seq, id=r.id, description=''))
    SeqIO.write(records, output_fasta, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA序列翻译为蛋白质')
    parser.add_argument('-i', '--input', required=True, help='输入DNA序列FASTA文件')  # 改 - 输入DNA序列FASTA文件
    parser.add_argument('-o', '--output', required=True, help='输出蛋白质序列FASTA文件')  # 改 - 输出蛋白质序列FASTA文件
    args = parser.parse_args()

    translate_dna(args.input, args.output)
