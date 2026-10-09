## Description

Apply three transformations to a single input DNA sequence:

- **Reverse:** Reverse the sequence order without changing the bases.
- **Complement:** Replace each base with its complement (A/T and C/G) without reversing the order.
- **Reverse complement:** Complement the sequence and then reverse its order.

## Usage

### Command line

```bash
python reverse_complement.py -i input.fasta -o output.fasta
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input FASTA file containing a single DNA sequence |
| `-o`, `--output` | Output FASTA file containing the reverse, complement, and reverse complement sequences |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste a DNA sequence in FASTA format.
- Only one DNA sequence is supported.

**Output:**

The page displays all three results (reverse, complement, and reverse complement), each with a **Copy** button. **Download** the complete results as a FASTA file.
