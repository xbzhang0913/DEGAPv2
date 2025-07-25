# ![DEGAP](DEGAP.png)
## Dynamic Elongation of a Genome Assembly Path

**DEGAP degaps gaps!**

**DEGAP v2.0** is an advanced gap-filling software that resolves gap regions by utilizing the dual advantages of accuracy and length of high-fidelity (HiFi) reads. This version significantly enhances the original DEGAP with new features, improved performance, and better usability.

## Features

### Core Functionality
- **Gap Filling**: Fill gaps between known sequences using HiFi reads
- **Contig Linking**: Link contigs to create scaffolds
- **Telomere Filling**: NEW - Fill telomeric regions at chromosome ends
- **Resume Capability**: NEW - Resume interrupted runs from specific rounds

### Key Improvements in v2.0
- **Enhanced Parameter System**: Migrated from getopt to argparse for better usability
- **K-mer Filtering**: Advanced k-mer based read filtering for improved accuracy
- **Parallel Processing**: Optimized parallel job control with customizable thread counts
- **Index Reuse**: Intelligent reads index caching to avoid rebuilding
- **Batch Processing**: AutoGapfiller.py for automated batch gap filling
- **LSF Integration**: Built-in LSF cluster job submission support

## Installation

### System Requirements

**DEGAP v2.0** has been developed with the Python programming language under a Linux environment. The following packages/tools are required:

- **Biopython** (version 1.80)
- **Pysam** (version 0.20.0)
- **Minimap2** (2.17-r941)
- **Hifiasm** (0.16.1-r375)
- **SAMtools** (version 1.6)
- **Seqkit** （version 2.8.0）
- **MUMmer** (4.0.0beta2) 


### Setup
```bash
git clone https://github.com/xbzhang0913/DEGAPv2.git
cd DEGAPv2
# Ensure all dependencies are installed
```

## Tutorial

**DEGAP v2.0** provides 3 modes (GapFiller, CtgLinker, TelFiller) to fill the gaps.

### GapFiller Mode
GapFiller is designed for one specific gap, of which the exactly left and right sequences between this gap are already known. GapFiller is designed to fill the gap by elongating the sequence from only one direction.

```bash
python DEGAP.py --mode gapfiller \
    --seqleft ./path/gapLeftSequence.fasta \
    --seqright ./path/gapRightSequence.fasta \
    --reads ./path/HiFiReads.fasta \
    -o ./path/Output \
    --flag left 
```

### CtgLinker Mode
CtgLinker is designed for the assembly with gaps. If sequences in an assembly are unordered, CtgLinker can try to not only fill the gaps, but also link the potentially neighbored sequences by elongating the edges.

```bash
python DEGAP.py --mode ctglinker \
    --ctgseq ./path/contigs.fasta \
    --reads ./path/HiFiReads.fasta \
    --out ./path/Output \
    --filterDepth 0.3
```

### TelFiller Mode (NEW in v2.0)
TelFiller is designed for filling telomeric regions at chromosome ends using HiFi reads.

```bash
python DEGAP.py --mode telfiller \
    --reads ./path/HiFiReads.fasta \
    --seqleft ./path/startSequence.fasta \
    --seqright ./path/endSequences.fasta \
    --flag left \
    -o ./path/Output 
```

### Filter Options
**DEGAP v2.0** also provides parameter `--filterDepth num` to filter the HiFi reads and `--filterDepthOnly` to only filter the HiFi reads but not execute the elongation process.

```bash
# Filter reads only
python DEGAP.py --mode ctglinker \
    --ctgseq ./path/contigs.fasta \
    --reads ./path/HiFiReads.fasta \
    --out ./path/Output \
    --filterDepthOnly 0.3

# Filter reads in gapfiller mode
python DEGAP.py --mode gapfiller \
    --seqleft ./path/gapLeftSequence.fasta \
    --seqright ./path/gapRightSequence.fasta \
    --reads ./path/HiFiReads.fasta \
    -o ./path/Output \
    --filterDepthOnly 0.3
```

### Advanced Parameters

#### K-mer Filtering Options
- `--kmer_size` / `-ks`: K-mer size for filtering reads (default: 41)
- `--kmer_num` / `-kn`: Number of k-mers to use for filtering (default: 10)
- `--kmer_length` / `-kl`: Proportion of mean read length for k-mer extraction (default: 0.1)

#### Performance Options
- `-j`: Number of parallel jobs for processing reads (default: 100)
- `-t` / `--thread`: Number of threads (default: 20)
- `--filterDepth`: Filter HiFi reads by mapped depth threshold

#### Resume Options (NEW)
- `--resume <round_number>`: Resume from specific round (e.g., --resume 118)
- `--resume_auto`: Automatically resume from last interrupted round

#### Output Control
- `--remove`: File cleanup level (1: final only, 2: basic results, 3: all files)
- `--edge`: Maximum edge length for misassembly detection (default: 500)
- `--MaximumExtensionLength`: Stop extension at specified length (default: 1000000)

### Batch Processing with AutoGapfiller

For processing multiple gaps automatically:

```bash
python AutoGapfiller.py \
    --reads HiFi_reads.fasta \
    --genome genome.fasta \
    --mode gapfiller \
    --flag all \
    -o output_directory \
    --batch 10
```

#### AutoGapfiller Features
- Automatic gap detection and filtering
- Batch job script generation
- LSF cluster integration
- Parallel and serial execution modes
- Comprehensive result monitoring

## Output Files

### Standard Output
- `*.final.fa`: Final assembled sequences
- `*.log`: Detailed processing logs
- `*.stat`: Assembly statistics and N50 metrics

### Intermediate Files (depending on --remove setting)
- `reads.idx`: Reads index database
- `reads_part/`: Split reads files for parallel processing
- `project/`: Round-by-round extension results
- `jobScripts/`: Generated batch processing scripts (AutoGapfiller)

## Algorithm Overview

DEGAP v2.0 uses an iterative extension approach:

1. **Initialization**: Build reads index and split files for parallel processing
2. **Seed Selection**: Use provided sequences as extension seeds
3. **Read Mapping**: Map HiFi reads to current sequences using minimap2
4. **Extension**: Extend sequences using overlapping reads
5. **Validation**: Validate extensions and detect potential misassemblies
6. **Iteration**: Repeat until no further extension possible

### New in v2.0: Enhanced Overlap Detection
- Improved simple overlap detection
- Tandem repeat handling
- Better identity and length thresholds
- Direct sequence connection for overlapping flanking sequences

## Performance Considerations

- **Memory**: Scales with reads file size; index caching reduces memory usage
- **CPU**: Highly parallelizable; use `-j` and `-t` parameters to optimize
- **Storage**: Intermediate files can be large; use `--remove 1` for minimal storage
- **Cluster**: Use AutoGapfiller with LSF for large-scale processing

## Troubleshooting

### Common Issues
1. **Index Building Fails**: Check reads file format and available disk space
2. **No Extension Found**: Verify sequence quality and HiFi read coverage
3. **Memory Issues**: Reduce parallel jobs (`-j`) or use reads filtering
4. **Resume Fails**: Ensure output directory structure is intact

### Resume Functionality
If a run is interrupted:
```bash
# Resume from specific round
python DEGAP.py --resume 118 [other parameters]

# Auto-resume from last round
python DEGAP.py --resume_auto [other parameters]
```

## How to cite?

If you use DEGAP v2.0 in your research, please cite:

Huang, Y. et al. DEGAP: Dynamic Elongation of a Genome Assembly Path. Briefings in Bioinformatics, Volume 25, Issue 3, May 2024, bbae194. [DOI: 10.1093/bib/bbae194](https://doi.org/10.1093/bib/bbae194)

For questions, comments and suggestions, please contact <jzhang@mail.hzau.edu.cn>.

## Changelog

### v2.0 vs v1.0
- **NEW**: Telomere filling mode
- **NEW**: Resume capability for interrupted runs
- **NEW**: K-mer based read filtering
- **NEW**: AutoGapfiller for batch processing
- **NEW**: LSF cluster integration
- **IMPROVED**: Parameter system with argparse
- **IMPROVED**: Reads index caching and reuse
- **IMPROVED**: Parallel processing optimization
- **IMPROVED**: Overlap detection algorithms
- **IMPROVED**: Error handling and logging
