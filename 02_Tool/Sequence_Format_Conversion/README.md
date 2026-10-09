## Description

Convert multiple sequence alignment files from one format to another. The tool uses the Biopython SeqIO module and supports 19 common sequence and alignment formats.

## Usage

### Command line

```bash
python input_seq.py -i input.aln -o output.fasta -if clustal -of fasta
```

| Parameter | Description |
|-----------|-------------|
| `-i`, `--input` | Input multiple sequence alignment file |
| `-o`, `--output` | Output file |
| `-if`, `--input-format` | Input format, such as clustal, fasta, or phylip |
| `-of`, `--output-format` | Output format, such as fasta, clustal, or phylip |

### Web tool

**Input options:**

- **File upload:** Upload an alignment file up to 1 MB.
- **Text input:** Paste the file contents.

The page provides two dropdown menus:

1. **Input format:** Select the format of the input file.
2. **Output format:** Select the desired output format.

**Output:**

The converted result is displayed on the page and can be **copied** or **downloaded**.

## Supported formats (19)

| Format | Description |
|--------|-------------|
| clustal | ClustalW multiple sequence alignment |
| embl | EMBL format |
| fasta | FASTA format |
| fasta-2line | FASTA with unwrapped sequences |
| fastq | FASTQ with quality scores |
| fastq-illumina | Illumina FASTQ |
| fastq-sanger | Sanger FASTQ |
| fastq-solexa | Solexa FASTQ |
| gb | GenBank format |
| genbank | GenBank format alias |
| nexus | NEXUS format |
| phylip | Interleaved PHYLIP |
| phylip-relaxed | Relaxed PHYLIP |
| phylip-sequential | Sequential PHYLIP |
| pir | PIR / NBRF format |
| qual | Quality score format |
| seqxml | SeqXML format |
| stockholm | Stockholm format |
| tab | Tab-delimited format |

## Notes

- One file is converted per run; batch conversion is not supported.
- Some conversions lose information. For example, converting GenBank to FASTA discards annotations. Choose the input and output formats accordingly.
