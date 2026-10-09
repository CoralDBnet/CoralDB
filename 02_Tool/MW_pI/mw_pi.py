#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : mw_pi.py
# @Date          : 2026-07-22
# @Description   : 计算蛋白质分子量、等电点等理化性质
"""
import argparse
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def _calc_one(record):
    seq = str(record.seq)
    analysis = ProteinAnalysis(seq)
    aa_str = ', '.join(f'{aa}:{round(pct * 100, 1)}%' for aa, pct in sorted(analysis.amino_acids_percent.items()))

    return {
        'id': record.id,
        'length': len(seq),
        'mw_da': round(analysis.molecular_weight(), 2),
        'mw_kda': round(analysis.molecular_weight() / 1000, 2),
        'pi': round(analysis.isoelectric_point(), 2),
        'extinction_280': round(analysis.molar_extinction_coefficient()[0], 2),
        'extinction_280_reduced': round(analysis.molar_extinction_coefficient()[1], 2),
        'instability_index': round(analysis.instability_index(), 2),
        'gravy': round(analysis.gravy(), 4),
        'aromaticity': round(analysis.aromaticity(), 4),
        'aa_percent': aa_str,
    }


def calc_mw_pi(input_fasta, output_tsv):
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    if not records:
        raise ValueError('No sequences found in the input file')
    results = [_calc_one(r) for r in records]

    header = ['ID', 'Length(aa)', 'MW(Da)', 'MW(kDa)', 'pI',
              'ExtCoeff_Ox(M-1cm-1)', 'ExtCoeff_Red(M-1cm-1)',
              'InstabilityIndex', 'GRAVY', 'Aromaticity', 'AA_Percent']
    with open(output_tsv, 'w', newline='') as f:
        # suppress newline conversion on Windows
        f.write('\t'.join(header) + '\n')
        for r in results:
            row = [
                str(r['id']),
                str(r['length']),
                str(r['mw_da']),
                str(r['mw_kda']),
                str(r['pi']),
                str(r['extinction_280']),
                str(r['extinction_280_reduced']),
                str(r['instability_index']),
                str(r['gravy']),
                str(r['aromaticity']),
                r['aa_percent'],
            ]
            f.write('\t'.join(row) + '\n')

    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='计算蛋白质分子量、等电点等理化性质')
    parser.add_argument('-i', '--input', required=True, help='输入蛋白序列FASTA文件')  # 改 - 输入蛋白序列FASTA文件
    parser.add_argument('-o', '--output', required=True, help='输出 TSV 结果文件')    # 改 - 输出 TSV 结果文件
    args = parser.parse_args()

    results = calc_mw_pi(args.input, args.output)
    print(f'Total sequences: {len(results)}')
    for r in results:
        print(f"\nID: {r['id']}")
        print(f"  Length: {r['length']} aa")
        print(f"  MW: {r['mw_da']} Da ({r['mw_kda']} kDa)")
        print(f"  pI: {r['pi']}")
        print(f"  ExtCoeff (oxidized): {r['extinction_280']} M-1cm-1")
        print(f"  ExtCoeff (reduced):  {r['extinction_280_reduced']} M-1cm-1")
        print(f"  Instability index: {r['instability_index']} (<40 stable)")
        print(f"  GRAVY: {r['gravy']}")
        print(f"  Aromaticity: {r['aromaticity']}")
