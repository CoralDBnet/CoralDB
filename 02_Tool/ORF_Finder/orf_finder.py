#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : orf_finder.py
# @Date          : 2026-07-22
# @Description   : DNA序列ORF查找，六框翻译，输出TSV结果
"""
import argparse
from Bio import SeqIO


CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}
START_CODONS = {'ATG', 'GTG', 'TTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
COMPLEMENT = str.maketrans('ATCGatcg', 'TAGCtagc')

MIN_ORF_LEN = 30  # 最少 10 个氨基酸


def _translate(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def find_orfs(input_fasta):
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    if len(records) != 1:
        raise ValueError('仅支持一条DNA序列，当前输入包含 %d 条' % len(records))

    seq = str(records[0].seq).upper()
    rev_seq = seq[::-1].translate(COMPLEMENT)

    results = []

    def _scan(frame_seq, strand, frame):
        for i in range(len(frame_seq) - 2):
            codon = frame_seq[i:i + 3]
            if codon not in START_CODONS:
                continue
            protein = []
            for j in range(i, len(frame_seq) - 2, 3):
                c = frame_seq[j:j + 3]
                if c in STOP_CODONS:
                    if len(protein) * 3 >= MIN_ORF_LEN:
                        abs_start = i + (frame - 1) if strand == '+' else len(seq) - (i + len(protein) * 3 + 3) + frame
                        abs_end = i + len(protein) * 3 + 2 + (frame - 1) if strand == '+' else len(seq) - i - frame + 1
                        results.append({
                            'strand': strand, 'frame': frame,
                            'start': abs_start + 1, 'end': abs_end + 1,
                            'length_aa': len(protein),
                            'protein': ''.join(protein),
                        })
                    break
                protein.append(CODON_TABLE.get(c, 'X'))

    for f in range(3):
        _scan(seq[f:], '+', f + 1)
        _scan(rev_seq[f:], '-', f + 1)

    results.sort(key=lambda x: -x['length_aa'])
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA序列ORF查找，六框翻译')
    parser.add_argument('-i', '--input', required=True, help='输入DNA序列FASTA文件')  # 改 - 输入DNA序列FASTA文件（仅支持单条）
    parser.add_argument('-o', '--output', required=True, help='输出TSV结果文件')  # 改 - 输出TSV结果文件路径
    args = parser.parse_args()

    orfs = find_orfs(args.input)
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write('\t'.join(['Rank', 'Strand/Frame', 'Start', 'End', 'Length_aa', 'Protein_Sequence']) + '\n')
        for i, o in enumerate(orfs, 1):
            f.write('\t'.join([
                str(i),
                f"{o['strand']}{o['frame']}",
                str(o['start']),
                str(o['end']),
                str(o['length_aa']),
                o['protein'],
            ]) + '\n')
    print(f"共找到 {len(orfs)} 个 ORF (>={MIN_ORF_LEN} bp)，结果已写入 {args.output}")
