# CNVAnno - ClinGen CNV Pathogenicity Annotation Tool

A comprehensive CNV (Copy Number Variant) pathogenicity annotation tool based on **ClinGen official cumulative scoring system**. This tool evaluates CNV pathogenicity following the ClinGen CNV Working Group guidelines.

## Features

- **ClinGen Official Scoring**: Implements the 5-section cumulative scoring system
- **Multi-database Integration**: Uses ClinGen HI/TS genes, curated regions, gnomAD CNV, DECIPHER population data, GenCC gene-disease associations, and Gencode exon coordinates
- **Exon-level Analysis**: Determines 5' vs 3' partial overlap for precise Section 2C/2D scoring
- **Evidence Tracking**: Provides detailed evidence fields showing which rules were applied
- **Manual Decision Markers**: Flags subjective evidence types requiring user input (de novo, segregation, inheritance)
- **Flexible Input**: Supports BED, VCF, CNVkit .cnr, and generic TSV formats
- **Dual Output**: Both TSV and JSON formats available

## Installation

### Option 1: Docker (Recommended)

```bash
# Build image
docker build -t cnvanno:latest .

# Run with data mounted
docker run -v /path/to/data:/app/data cnvanno:latest input.bed -o output.tsv -g GRCh38
```

### Option 2: Direct Installation

#### Requirements

```bash
pip install pandas openpyxl intervaltree
```

Or use the provided requirements file:

```bash
pip install -r requirements.txt
```

### Database Files

Download the following databases to the `data/` directory:

| Database | File Name | Description |
|----------|-----------|-------------|
| ClinGen HI Genes | `ClinGen_haploinsufficiency_gene_{build}.bed` | Haploinsufficiency gene scores (0-40) |
| ClinGen TS Genes | `ClinGen_triplosensitivity_gene_{build}.bed` | Triplosensitivity gene scores (0-40) |
| ClinGen Curated Regions | `ClinGen_region_curation_list_{build}..tsv` | Pathogenic/benign curated CNV regions |
| GenCC | `gencc-submissions.xlsx` | Gene-disease associations with inheritance modes |
| Gencode Exons | `Gencode.{build}.cnvkit.target.bed` | Exon coordinates for partial overlap analysis |
| gnomAD CNV | `gnomad.v4.1.cnv.all.{build}.vcf.gz` | Population CNV frequencies |
| DECIPHER | `population_cnv_{build}.txt.gz` | Population CNV frequencies from DECIPHER |
| UCSC Cytoband | `cytoBand_{build}.txt` | Chromosome band coordinates for ISCN formatting |

Replace `{build}` with `GRCh37` or `GRCh38` based on your genome build.

**Download Sources:**
- ClinGen: https://clinicalgenome.org/
- GenCC: https://thegencc.org/
- Gencode: https://www.gencodegenes.org/
- gnomAD CNV: https://gnomad.broadinstitute.org/
- DECIPHER: https://decipher.sanger.ac.uk/
- UCSC Cytoband: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz (GRCh38) or https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz (GRCh37)

## Usage

### Basic Usage

```bash
python cnvanno.py input.bed -o output.tsv
```

### Full Options

```bash
python cnvanno.py input.bed \
    -d data/ \
    -g GRCh38 \
    -o output.tsv \
    -f tsv \
    -l INFO
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-d, --data-dir` | Database directory | `data` |
| `-g, --genome-build` | Genome build (GRCh37/GRCh38) | `GRCh38` |
| `-o, --output` | Output file (required) | - |
| `-f, --format` | Output format (tsv/json) | `tsv` |
| `-l, --log-level` | Logging level | `INFO` |

### Input Format

The simplest input format is a 4-column BED file:

```
chrom   start   end     status
chr1    1020119 1056116 Deletion
chr1    6785453 7769706 Deletion
chr7    5497217 5760091 Deletion
chr21   25734570 26573286 Amplification
chr22   18578203 18649330 Normal
```

**Supported formats:**
- **BED** (4 columns): chrom, start, end, status
- **VCF**: SVTYPE=DEL/DUP parsed from INFO field
- **CNVkit .cnr**: log2 ratio converted to status (log2<-0.6=Deletion, log2>0.4=Amplification)
- **Generic TSV**: Auto-detects column names

**Status values:**
- `Deletion` / `DEL` / `Loss`: For copy number loss
- `Amplification` / `DUP` / `Gain` / `Duplication`: For copy number gain
- `Normal`: No CNV detected

**Chromosome format:**
Both `chr1` and `1` formats are accepted.

## ClinGen Scoring System

### Classification Thresholds

| Total Score | Classification |
|-------------|----------------|
| ≤ -0.99 | Benign |
| -0.98 ~ 0.99 | Uncertain Significance (VUS) |
| ≥ 1.00 | Pathogenic |

### Section Scoring

#### Section 1: Genomic Content
- **1A**: Contains protein-coding genes → **Score: 0**
- **1B**: No protein-coding genes → **Score: -0.60**

#### Section 2: Gene/Region Overlap

**For Loss (Deletion):**
| Evidence | Description | Score |
|----------|-------------|-------|
| 2A | Full overlap with established HI gene/region (HI=3) | +1.00 |
| 2B | Partial overlap with HI gene/region | +0.45 |
| 2C | Partial 5' overlap affecting coding sequence | +0.90 |
| 2D | Partial 3' overlap (score varies by context) | +0.45~0.90 |
| 2E | Both breakpoints within same gene (PVS1-like) | +0.90 |
| 2F | Overlaps benign region (HI/TR=40) | -1.00 |
| 2H | Likely HI gene (HI=2) | +0.15 |

**For Gain (Amplification):**
| Evidence | Description | Score |
|----------|-------------|-------|
| 2A | Full overlap with established TS gene/region (TR=3) | +1.00 |
| 2B | Partial overlap with TS gene | +0.45 |
| 2I | Both breakpoints in HI gene (gene disruption) | +0.90 |
| 2K | One breakpoint in HI gene + phenotype match | +0.45 |
| 2F | Overlaps benign region (TR=40) | -1.00 |
| 2H | Likely TS gene (TR=2) | +0.15 |

#### Section 3: Gene Count

**For Loss:**
| Gene Count | Score |
|------------|-------|
| 0-1 | 0.00 |
| 2-4 | 0.15 |
| 5-9 | 0.30 |
| 10-24 | 0.45 |
| 25-49 | 0.60 |
| ≥50 | 0.90 |

**For Gain:**
| Gene Count | Score |
|------------|-------|
| 0-34 | 0.00 |
| 35-49 | 0.45 |
| ≥50 | 0.90 |

#### Section 4: Frequency/Database Evidence

| Evidence | Description | Score |
|----------|-------------|-------|
| 4O | Frequency ≥ 1% in population databases | -1.00 |
| 4O-like | Frequency ≥ 0.1% | -0.45 |

**Manual Decision Required (marked in output):**
- 4A-4D: De novo variants (requires phenotype confirmation)
- 4F-4H: Segregation analysis
- 4I-4K: Non-segregation evidence (negative scores)
- 4L-4N: Case-control studies

#### Section 5: Inheritance

Requires manual input for:
- 5A: De novo inheritance
- 5B-5C: Inherited from unaffected parent
- 5D: Segregation in family
- 5F-5H: No inheritance information available

## Output Fields

### TSV Output Columns

| Column | Description |
|--------|-------------|
| `#Chromosome`, `Start`, `End`, `Size`, `Status` | CNV basic information |
| `ISCN` | ISCN standard nomenclature (e.g., `seq[GRCh38] 1p36.33p36.11(827144_27147702)del`) |
| `Gene_Count` | Number of overlapping genes |
| `HI_Max` | Maximum HI score among overlapping genes |
| `TR_Max` | Maximum TR score among overlapping genes |
| `Max_Frequency` | Maximum population frequency |
| `Section1-5`, `Total_Score` | ClinGen section scores and total |
| `Evidence_1A-5` | Evidence usage flags (1=used, 0=not used) |
| `Dosage_Genes` | Overlapping dosage-sensitive genes |
| `Pathogenic_Regions` | Overlapping pathogenic curated regions |
| `Benign_Regions_Overlap` | Overlapping benign regions |
| `GenCC_AD_Genes` | GenCC genes with AD inheritance |
| `Classification` | Final classification result |
| `Reason` | Classification reasoning |
| `Evidence_Details` | Detailed evidence explanations |

### Evidence Field Definitions

| Field | Meaning |
|-------|---------|
| `Evidence_1A` | Contains protein-coding genes |
| `Evidence_1B` | No protein-coding genes |
| `Evidence_2A` | Full overlap with established pathogenic gene/region |
| `Evidence_2B` | Partial overlap with pathogenic gene/region |
| `Evidence_2C` | Partial 5' end overlap |
| `Evidence_2D` | Partial 3' end overlap |
| `Evidence_2E` | Both breakpoints within same gene |
| `Evidence_2F` | Overlaps benign region |
| `Evidence_2H` | Likely pathogenic gene (score=2) |
| `Evidence_2K` | One breakpoint in HI gene (Gain only) |
| `Evidence_3` | Gene count scoring applied |
| `Evidence_4O` | High population frequency |
| `Evidence_4A` | De novo evidence (requires manual input) |
| `Evidence_4L` | Case-control evidence (requires manual input) |
| `Evidence_5` | Inheritance evidence (requires manual input) |

### JSON Output

JSON output includes additional structured information:
- `exon_overlap_details`: List of affected exons with strand orientation
- `manual_decision_needed`: Flag for subjective evidence requiring user decision
- `dosage_genes`: List of all overlapping dosage-sensitive genes

## HI/TS Score Interpretation

| Score | Meaning |
|-------|---------|
| 3 | Sufficient evidence for dosage sensitivity (established) |
| 2 | Likely dosage sensitive |
| 1 | Little evidence for dosage sensitivity |
| 0 | No evidence available |
| 30 | No evidence available (gene not curated) |
| 40 | Dosage sensitivity unlikely (benign evidence) |

## Example

```bash
# Annotate CNVs with ClinGen scoring
python cnvanno.py test/test_input.bed -o results.tsv -g GRCh38

# JSON output for programmatic access
python cnvanno.py test/test_input.bed -o results.json -f json -g GRCh38
```

## Limitations

1. **Subjective Evidence**: The following evidence types require manual user decision and are marked in output:
   - De novo variants with phenotype information
   - Segregation analysis
   - Inheritance patterns from family history
   - Case-control study significance

2. **Section 5 Default**: Inheritance scoring defaults to 0 when no data available

3. **Partial Overlap**: Requires Gencode exon data for precise 5'/3' determination

## References

- ClinGen CNV Working Group: https://clinicalgenome.org/cnv/
- ClinGen CNV Calculator: https://cnvcalc.clinicalgenome.org/
- Riggs et al. (2020). Technical standards for the interpretation and reporting of constitutional copy-number variants. *Genetics in Medicine*, 22(2), 245-257.

## License

MIT License
