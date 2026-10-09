#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : reverse_complement.py
# @Date          : 2026-07-21
# @Description   : 单条DNA序列的反向、互补、反向互补
"""
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def reverse_complement(input_fasta, output_fasta):
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    if len(records) != 1:
        raise ValueError('仅支持一条序列，当前输入包含 %d 条' % len(records))

    r = records[0]
    records = [
        SeqRecord(r.seq[::-1], id=r.id + '_reverse', description=''),
        SeqRecord(r.seq.complement(), id=r.id + '_complement', description=''),
        SeqRecord(r.seq.reverse_complement(), id=r.id + '_reverse_complement', description=''),
    ]
    SeqIO.write(records, output_fasta, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA序列反向、互补、反向互补')
    parser.add_argument('-i', '--input', required=True, help='输入单条DNA序列FASTA文件')  # 改 - 输入单条DNA序列FASTA文件
    parser.add_argument('-o', '--output', required=True, help='输出结果FASTA文件')  # 改 - 输出包含反向、互补、反向互补结果的FASTA文件
    args = parser.parse_args()

    reverse_complement(args.input, args.output)
