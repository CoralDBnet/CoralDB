## Description

Batch-extract sequences from PEP and/or CDS files using a user-provided list of gene IDs.

## Usage

### Command line

```bash
# Extract PEP sequences only
python geneid_get_seq_v2.py -g gene_id.txt -p Adu.pep -op pep_result.fasta

# Extract CDS sequences only
python geneid_get_seq_v2.py -g gene_id.txt -c Adu.cds -oc cds_result.fasta

# Extract both PEP and CDS sequences
python geneid_get_seq_v2.py -g gene_id.txt -p Adu.pep -op pep_result.fasta -c Adu.cds -oc cds_result.fasta
```

| Parameter | Description |
|-----------|-------------|
| `-g`, `--genes` | Required gene ID file, with one ID per line |
| `-p`, `--pep` | Optional input PEP sequence file |
| `-c`, `--cds` | Optional input CDS sequence file |
| `-op`, `--output_pep` | PEP output file; used with `-p` |
| `-oc`, `--output_cds` | CDS output file; used with `-c` |

### Web tool

**Input options:**

- **Gene IDs:** Paste one gene ID per line into the text box.
- **Species:** Select the target species from the dropdown list.
- **Sequence type:** Select CDS or PEP sequences.

**Output:**

Extracted sequences are displayed on the page, each with a **Copy** button. **Download** the complete results as a FASTA file.

## Notes

- Gene IDs must exactly match the FASTA headers in the database.
- CDS and PEP sequences are extracted independently.
- Unmatched gene IDs are reported as not found.
- At least one of `-p` or `-c` must be provided.
