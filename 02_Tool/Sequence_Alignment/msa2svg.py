#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : msa2svg.py
# @Date          : 2026-07-21 14:41:19
# @Last Modified : 2026-07-21 14:41:19
# @Description   : 将多序列比对结果转换为svg格式
"""
import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# === 氨基酸颜色 (Clustal X 风格) ===
AA_COLORS = {
    'A': '#80A0F0',
    'R': '#F01505', 'K': '#F01505',
    'N': '#00FF00', 'Q': '#00FF00',
    'D': '#C048C0', 'E': '#C048C0',
    'C': '#F08080',
    'G': '#F09048',
    'H': '#15A4A4',
    'I': '#80A0F0', 'L': '#80A0F0', 'M': '#80A0F0', 'V': '#80A0F0',
    'F': '#80A0F0', 'W': '#80A0F0',
    'Y': '#15A4A4',
    'P': '#CCCC00',
    'S': '#00FF00', 'T': '#00FF00',
    '-': '#FFFFFF', '.': '#FFFFFF',
    'B': '#E0E0E0', 'Z': '#E0E0E0', 'X': '#E0E0E0',
}

# === 核苷酸颜色 ===
NT_COLORS = {
    'A': '#64F73F', 'T': '#FFB340', 'U': '#FFB340',
    'C': '#FF4040', 'G': '#EB413C',
    '-': '#FFFFFF', '.': '#FFFFFF',
    'N': '#E0E0E0',
}

CHAR_W = 14
CHAR_H = 18
LEFT_MARGIN = 180


def _detect_seq_type(alignment):
    for record in alignment:
        for letter in record.seq:
            if letter.upper() not in 'ATGCRYKMSWBDHVNU-.X':
                return 'protein'
    return 'dna'


def _generate_svg(aligned, output_svg):
    seq_type = _detect_seq_type(aligned)
    colors = AA_COLORS if seq_type == 'protein' else NT_COLORS

    n_seq = len(aligned)
    aln_len = len(aligned[0].seq)

    ids = []
    for record in aligned:
        rid = record.id
        if len(rid) > 22:
            rid = rid[:20] + '..'
        ids.append(rid)

    svg_w = LEFT_MARGIN + aln_len * CHAR_W + 20
    top_margin = 25
    svg_h = top_margin + n_seq * CHAR_H + 30 + CHAR_H + 20

    parts = []
    # -- SVG 头部 --
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {svg_w} {svg_h}" width="100%" height="100%">'
    )
    parts.append('<style>')
    parts.append('  .cell { cursor: pointer; }')
    parts.append('  .cell:hover { stroke: #000; stroke-width: 1.5; }')
    parts.append('  .seq-label { font-family: monospace; font-size: 12px; fill: #333; }')
    parts.append('  .pos-label { font-family: monospace; font-size: 10px; fill: #999; }')
    parts.append('  .residue { font-family: monospace; font-size: 10px; fill: #000; pointer-events: none; }')
    parts.append('</style>')
    parts.append(f'<rect width="{svg_w}" height="{svg_h}" fill="#fff"/>')

    # -- 顶部刻度 (每 10 个) --
    for i in range(0, aln_len, 10):
        x = LEFT_MARGIN + i * CHAR_W + CHAR_W / 2
        parts.append(
            f'<text x="{x}" y="{top_margin - 5}" '
            f'class="pos-label" text-anchor="middle">{i + 1}</text>'
        )

    # -- 序列行 --
    for row, record in enumerate(aligned):
        y = top_margin + row * CHAR_H
        parts.append(
            f'<text x="5" y="{y + CHAR_H * 0.7}" class="seq-label">{ids[row]}</text>'
        )
        for col, letter in enumerate(str(record.seq)):
            x = LEFT_MARGIN + col * CHAR_W
            fill = colors.get(letter.upper(), '#E0E0E0')
            parts.append(
                f'<rect class="cell" x="{x}" y="{y}" '
                f'width="{CHAR_W}" height="{CHAR_H}" fill="{fill}">'
                f'<title>{ids[row]}  pos:{col + 1}  {letter}</title>'
                f'</rect>'
            )
            parts.append(
                f'<text class="residue" x="{x + CHAR_W / 2}" '
                f'y="{y + CHAR_H * 0.7}" text-anchor="middle">{letter}</text>'
            )

    # -- 保守性行 --
    con_y = top_margin + n_seq * CHAR_H + 30
    parts.append(
        f'<text x="5" y="{con_y + CHAR_H * 0.7}" '
        f'class="seq-label">Conservation</text>'
    )
    for col in range(aln_len):
        col_letters = [str(record.seq)[col].upper() for record in aligned]
        unique = set(col_letters) - {'-', '.'}
        if not unique:
            ratio = 0.0
        else:
            ratio = 1.0 - (len(unique) - 1) / max(len(col_letters), 1)
        v = int(255 * (1 - ratio))
        fill = f'rgb({v},{v},{v})'
        x = LEFT_MARGIN + col * CHAR_W
        parts.append(
            f'<rect class="cell" x="{x}" y="{con_y}" '
            f'width="{CHAR_W}" height="{CHAR_H}" fill="{fill}">'
            f'<title>pos:{col + 1}  保守性:{ratio:.0%}</title>'
            f'</rect>'
        )

    parts.append('</svg>')

    with open(output_svg, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))


def _progressive_msa(records):
    if len(records) <= 1:
        return [SeqRecord(Seq(str(r.seq)), id=r.id, description='') for r in records]

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    original = [str(r.seq) for r in records]
    ids = [r.id for r in records]

    old = [list(original[0])]

    for i in range(1, len(records)):
        center_clean = ''.join(c for c in old[0] if c != '-')
        alns = aligner.align(center_clean, original[i])
        best = alns[0]
        cen_aln = best[0]
        new_aln = best[1]

        clean_to_full = [pos for pos, ch in enumerate(old[0]) if ch != '-']

        new_old = []
        for seq in old:
            new_seq = []
            ci = 0
            for ch in cen_aln:
                if ch == '-':
                    new_seq.append('-')
                else:
                    new_seq.append(seq[clean_to_full[ci]])
                    ci += 1
            new_old.append(new_seq)
        new_old.append(list(new_aln))
        old = new_old

    return [SeqRecord(Seq(''.join(s)), id=rid, description='') for s, rid in zip(old, ids)]


def multiple_sequence_alignment(input_fasta, output_svg):
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    if len(records) < 2:
        raise ValueError('至少需要 2 条序列才能进行多序列比对')

    aligned = _progressive_msa(records)
    _generate_svg(aligned, output_svg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='多序列比对并输出SVG可视化')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件')  # 改 - 输入FASTA文件路径
    parser.add_argument('-o', '--output', required=True, help='输出SVG文件')  # 改 - 输出SVG文件路径
    args = parser.parse_args()

    multiple_sequence_alignment(args.input, args.output)

