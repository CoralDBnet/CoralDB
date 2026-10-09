## Description

Align multiple sequences, construct a phylogenetic tree using UPGMA, and export an interactive SVG tree. Bootstrap analysis is supported, with support values displayed at branch nodes.

## Usage

### Command line

```bash
python msa2svg_v3.py -i input.fasta -o output.svg
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input FASTA file containing at least two sequences |
| `-o`, `--output` | Output phylogenetic tree in SVG format |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste sequences in FASTA format.
- Provide between 2 and 50 sequences.

**Output:**

The page displays an SVG phylogenetic tree with:

- **Panning:** Drag a blank area to move the view.
- **Zooming:** Use the mouse wheel to zoom in or out.
- **Bootstrap values:** Orange values at branch nodes.
- Sequence names next to the leaf nodes.

**Download** the SVG image.

## Notes

- At least two sequences are required.
- UPGMA assumes a molecular clock (equal evolutionary rates) and is suitable for closely related species.
- The default bootstrap analysis uses 100 replicates; higher values indicate stronger branch support.
