## Description

Perform progressive multiple sequence alignment (MSA) and visualize the alignment as a colored matrix in SVG format. Both amino acid and nucleotide sequences are supported.

## Usage

### Command line

```bash
python msa2svg.py -i input.fasta -o output.svg
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input FASTA file containing at least two sequences |
| `-o`, `--output` | Output SVG file |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste sequences in FASTA format.
- Provide between 2 and 50 sequences.

**Output:**

The page displays an SVG alignment with:

- A colored matrix in which each residue is represented by a colored square.
- A conservation track at the bottom; darker gray indicates greater conservation.
- A position ruler at the top, labeled every 10 positions.
- Position details available on mouse hover.

**Download** the SVG image.

## Notes

- At least two sequences are required.
- The tool automatically identifies amino acid or nucleotide sequences and applies the corresponding color scheme.
- The global alignment algorithm is suitable for aligning homologous sequences.
