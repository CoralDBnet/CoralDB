## Description

Calculate physicochemical properties of protein sequences, including molecular weight, theoretical isoelectric point (pI), extinction coefficients for reduced and oxidized states, instability index, and grand average of hydropathicity (GRAVY). The tool uses the Biopython ProtParam module and supports batch processing of multiple sequences.

## Usage

### Command line

```bash
python mw_pi.py -i input.fasta -o result.tsv
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input protein sequences in FASTA format; multiple sequences are supported |
| `-o`, `--output` | Output TSV file |

### Web tool

**Input options:**

- **File upload:** Upload a FASTA file up to 1 MB.
- **Text input:** Paste protein sequences in FASTA format.

**Output:**

Calculated properties for all sequences are displayed in a table:

| Field | Description |
|-------|-------------|
| ID | Sequence identifier |
| Length(aa) | Sequence length in amino acids |
| MW(Da) / MW(kDa) | Molecular weight |
| pI | Theoretical isoelectric point |
| ExtCoeff_Ox/Red | Extinction coefficients for oxidized and reduced states |
| InstabilityIndex | Instability index; values below 40 indicate a stable protein |
| GRAVY | Grand average of hydropathicity; positive values indicate hydrophobicity and negative values indicate hydrophilicity |
| Aromaticity | Aromaticity |
| AA_Percent | Percentage amino acid composition |

**Download** the complete results as a TSV file or **copy** the table contents.

## Notes

- Sequences containing nonstandard amino acids (such as B, Z, or X) cannot be evaluated accurately.
