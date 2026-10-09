## Description

Scan all six reading frames of an input DNA sequence to identify open reading frames (ORFs), sorted by length in descending order.

## Usage

### Command line

```bash
python orf_finder.py -i input.fasta -o result.tsv
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input DNA sequence in FASTA format; only one sequence is supported |
| `-o`, `--output` | Path to the output TSV file |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste a DNA sequence in FASTA format.
- Only one DNA sequence is supported.

**Output:**

All identified ORFs (at least 30 bp) are displayed in a table:

| Field | Description |
|-------|-------------|
| Rank | ORFs ranked by length in descending order |
| Strand/frame | +1 to +3 for the forward strand, or -1 to -3 for the reverse strand |
| Start/end | Absolute coordinates on the sequence |
| Length | Number of amino acids |
| Protein sequence | Translated amino acid sequence |

- **Copy** the protein sequence of an individual ORF.
- **Download** the complete results as a TSV file.

## Notes

- Only one DNA sequence is supported.
- Start codons: ATG / GTG / TTG; stop codons: TAA / TAG / TGA.
- The minimum ORF length is 30 bp (10 amino acids).
- Results are empty if no start codon or sufficiently long ORF is found.
