## Description

Translate CDS sequences (protein-coding DNA sequences) into PEP protein sequences (amino acid sequences).

## Usage

### Command line

```bash
python translate.py -i input.fasta -o output.fasta
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input DNA sequences in FASTA format |
| `-o`, `--output` | Output protein sequences in FASTA format |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste CDS sequences in FASTA format.

**Output:**

Translated protein sequences are displayed on the page, each with a **Copy** button. **Download** the complete results as a FASTA file. Sequence IDs are preserved.

## Notes

- Input sequences must be coding DNA sequences (CDS), not genomic sequences.
- Translation does not stop at stop codons (`to_stop=False`), which is suitable for complete CDS sequences.
