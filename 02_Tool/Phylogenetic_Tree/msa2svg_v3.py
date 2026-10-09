#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @File Name     : msa2svg.py
# @Date          : 2026-07-21 14:41:19
# @Last Modified : 2026-07-21 14:41:19
# @Description   : 多序列比对后构建系统发育树，输出可交互SVG
"""
import argparse
import random

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class _Node:
    __slots__ = ('name', 'left', 'right', 'height', 'bootstrap', 'y')

    def __init__(self, name=None, left=None, right=None, height=0.0):
        self.name = name
        self.left = left
        self.right = right
        self.height = height
        self.bootstrap = None


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


def _calc_distance(seq1, seq2):
    matches = 0
    total = 0
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':
            continue
        total += 1
        if a == b:
            matches += 1
    return 1.0 - matches / total if total > 0 else 1.0


def _leaf_names(node):
    if node.name:
        return frozenset([node.name])
    return _leaf_names(node.left) | _leaf_names(node.right)


def _upgma_from_dist(dist, names):
    n = len(names)
    if n == 1:
        return _Node(name=names[0])

    clusters = [(1, {i}, _Node(name=names[i])) for i in range(n)]

    def _cdist(ci, cj):
        li, lj = clusters[ci][1], clusters[cj][1]
        total = 0.0
        for a in li:
            for b in lj:
                total += dist[(min(a, b), max(a, b))]
        return total / (len(li) * len(lj))

    while len(clusters) > 1:
        min_d = float('inf')
        mi = mj = -1
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = _cdist(i, j)
                if d < min_d:
                    min_d = d
                    mi, mj = i, j

        si, leaves_i, node_i = clusters[mi]
        sj, leaves_j, node_j = clusters[mj]
        new_node = _Node(left=node_i, right=node_j, height=min_d / 2)

        del clusters[max(mi, mj)]
        del clusters[min(mi, mj)]
        clusters.append((si + sj, leaves_i | leaves_j, new_node))

    return clusters[0][2]


def _upgma(aligned):
    seqs = [str(r.seq) for r in aligned]
    names = [r.id for r in aligned]
    n = len(seqs)
    if n == 1:
        return _Node(name=names[0])

    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = _calc_distance(seqs[i], seqs[j])

    return _upgma_from_dist(dist, names)


def _bootstrap(aligned, tree, n_replicates=100):
    seqs = [str(r.seq) for r in aligned]
    names = [r.id for r in aligned]
    n_seq = len(seqs)
    aln_len = len(seqs[0])

    if n_seq <= 2 or aln_len == 0:
        return  # 无需 bootstrap，直接返回

    def _get_bips(node):
        if node.name:
            return []
        return [(_leaf_names(node), node)] + _get_bips(node.left) + _get_bips(node.right)

    orig_bips = _get_bips(tree)
    counts = {ls: 0 for ls, _ in orig_bips}

    for _ in range(n_replicates):
        sampled = [random.randint(0, aln_len - 1) for _ in range(aln_len)]
        boot_seqs = [''.join(s[i] for i in sampled) for s in seqs]

        dist = {}
        for i in range(n_seq):
            for j in range(i + 1, n_seq):
                dist[(i, j)] = _calc_distance(boot_seqs[i], boot_seqs[j])

        boot_tree = _upgma_from_dist(dist, names)
        boot_sets = set()

        def _collect(node):
            if node.name:
                return
            boot_sets.add(_leaf_names(node))
            _collect(node.left)
            _collect(node.right)

        _collect(boot_tree)

        for ls in counts:
            if ls in boot_sets:
                counts[ls] += 1

    for ls, node in orig_bips:
        node.bootstrap = round(counts[ls] / n_replicates * 100)


def _layout(root):
    leaves = []

    def collect(node):
        if node.name:
            leaves.append(node)
        else:
            collect(node.left)
            collect(node.right)
    collect(root)

    ROW_H = 22
    for i, leaf in enumerate(leaves):
        leaf.y = (i + 1) * ROW_H

    def assign_y(node):
        if node.name:
            return
        assign_y(node.left)
        assign_y(node.right)
        node.y = (node.left.y + node.right.y) / 2
    assign_y(root)

    return leaves


def _draw_tree(root, leaves, output_svg):
    max_h = root.height
    scale = 400 / max_h if max_h > 0 else 400
    ml = 30
    mt = 20
    mr = 200

    last_leaf = leaves[-1]
    svg_w = ml + int(max_h * scale) + mr
    svg_h = int(last_leaf.y) + mt * 2

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_w} {svg_h}" width="100%" height="100%">',
        '<defs>',
        '<style>',
        '  .branch { stroke: #333; stroke-width: 1.5; fill: none; }',
        '  .label { font-family: monospace; font-size: 10px; fill: #333; }',
        '  .bs-text { font-family: monospace; font-size: 9px; fill: #E65100; text-anchor: end; }',
        '</style>',
        '</defs>',
        f'<rect id="bg" width="{svg_w}" height="{svg_h}" fill="#fff"/>',
        '<g id="tree-g" transform="translate(0,0)">',
    ]

    def draw(node):
        x = ml + int((max_h - node.height) * scale)
        y = node.y + mt

        if node.name:
            rid = node.name
            if len(rid) > 22:
                rid = rid[:20] + '..'
            parts.append(f'<text class="label" x="{x + 8}" y="{y + 4}">{rid}</text>')
        else:
            draw(node.left)
            draw(node.right)
            parts.append(f'<line class="branch" x1="{x}" y1="{node.left.y + mt}" x2="{x}" y2="{node.right.y + mt}"/>')
            lx = ml + int((max_h - node.left.height) * scale)
            parts.append(f'<line class="branch" x1="{lx}" y1="{node.left.y + mt}" x2="{x}" y2="{node.left.y + mt}"/>')
            rx = ml + int((max_h - node.right.height) * scale)
            parts.append(f'<line class="branch" x1="{rx}" y1="{node.right.y + mt}" x2="{x}" y2="{node.right.y + mt}"/>')
            if node.bootstrap is not None:
                parts.append(f'<text class="bs-text" x="{x - 10}" y="{y - 4}">{node.bootstrap}</text>')

    draw(root)

    parts.append('</g>')
    parts.append('<script type="text/javascript"><![CDATA[')
    parts.append('var bg=document.getElementById("bg");')
    parts.append('var g=document.getElementById("tree-g");')
    parts.append('var px=0,py=0,dr=0,sx,sy;')
    parts.append('bg.addEventListener("mousedown",function(e){if(e.target!==bg)return;dr=1;sx=e.clientX-px;sy=e.clientY-py;bg.style.cursor="grabbing";});')
    parts.append('document.addEventListener("mousemove",function(e){if(!dr)return;px=e.clientX-sx;py=e.clientY-sy;var z=g.getAttribute("data-zoom")||1;g.setAttribute("transform","translate("+px+","+py+") scale("+z+")");});')
    parts.append('document.addEventListener("mouseup",function(){dr=0;bg.style.cursor="default";});')
    parts.append('bg.addEventListener("wheel",function(e){e.preventDefault();var z=parseFloat(g.getAttribute("data-zoom")||1);z+=e.deltaY>0?-0.1:0.1;z=Math.max(0.3,Math.min(3,z));g.setAttribute("data-zoom",z);g.setAttribute("transform","translate("+px+","+py+") scale("+z+")");});')
    parts.append(']]></script>')
    parts.append('</svg>')

    with open(output_svg, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))


def multiple_sequence_alignment(input_fasta, output_svg):
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    if len(records) < 2:
        raise ValueError('至少需要 2 条序列才能进行多序列比对')

    aligned = _progressive_msa(records)
    tree = _upgma(aligned)
    _bootstrap(aligned, tree)
    leaves = _layout(tree)
    _draw_tree(tree, leaves, output_svg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='多序列比对后构建系统发育树，输出可交互SVG')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件')  # 改 - 输入FASTA文件路径
    parser.add_argument('-o', '--output', required=True, help='输出系统发育树SVG文件')  # 改 - 输出系统发育树SVG文件路径
    args = parser.parse_args()

    multiple_sequence_alignment(args.input, args.output)
