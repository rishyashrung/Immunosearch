# Immunosearch

A pipeline for filtering and classifying MHC peptides from mass spectrometry data, especially useful for immunopeptidomics research.

## Overview

Immunosearch processes peptide data from PEAKS searches and Gibbs clustering outputs to identify and categorize peptides based on their potential origin. The pipeline can identify:

- Known MHC-presented peptides
- Single Amino Acid Variants (SAAVs)
- Peptides matching the six-frame translated human genome
- Proteasome-catalyzed peptide spliced (PCPS) variants

## Installation

### Requirements

- Python 3.7+
- BLAST+ suite
- SeqKit
- Multiple Python packages (see environment.yml)

### Setup with Conda

```bash
# Clone the repository
git clone https://github.com/rishyashrung/Immunosearch.git
cd Immunosearch

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate immunosearch
```

## Database Setup

Immunosearch requires several databases to be available. These should be placed in a single directory:

1. **APD_Hs_all** - Known human MHC-presented peptides database
   - Download from: https://peptideatlas.org/builds/human/hla/202311/APD_Hs_all.fasta
2. **human_canonical** - Human canonical protein database
3. **sap_db** - Single amino acid polymorphism database
   - Download from: http://119.3.41.228/dbSAP/download.html
4. **human_6FT_m.fasta** - Six-frame translated human genome
   - Generate from human reference genome with:
   ```bash
   seqkit translate db/GRCh38.p14_genomic.fna --append-frame -x -f 6 -M -m 8 --transl-table 1 -s > db/human_6FT_m.fasta
   ```

Each database (except the FASTA file) should be formatted using `makeblastdb`:

```bash
makeblastdb -in APD_Hs_all.fasta -dbtype prot -out APD_Hs_all
makeblastdb -in human_canonical.fasta -dbtype prot -out human_canonical
makeblastdb -in sap_db.fa -dbtype prot -out sap_db
```

## Usage

```bash
python test_immuno.py -w /path/to/workdir \
                     -i /path/to/input.xlsx \
                     -m 1 \
                     -g 1,2,3 \
                     -o /path/to/output \
                     -f output_filename \
                     -d /path/to/database_folder
```

### Parameters

- `-w, --work_dir`: Working directory where runtime files are created
- `-i, --input_file_path`: Input Excel file with PEAKS search result (sheet 1) and Gibbs clustering output (sheet 2)
- `-m, --MHC_class`: MHC class (1, 2, or E)
  - Class 1: selects peptides of 8-11 amino acids
  - Class 2: selects peptides of 12-17 amino acids
  - Class E: selects peptides of 8-15 amino acids
  - leave blank to skip filtering by length
- `-g, --gibbs_cluster`: Comma-separated list of Gibbs cluster values to select
- `-o, --output_file_path`: Output directory
- `-f, --file_name`: Excel file name (without extension)
- `-d, --db_path`: Path to database folder with all database files
- `-c, --cleanup`: Flag to cleanup all files in the working directory from previous runs (optional)

## Input Format

The input Excel file should contain:
- Sheet 1: PEAKS search results with at least a "Peptide" and "Found By" column
- Sheet 2: Gibbs clustering results named "gibbs_clustering" with at least "Sequence" and "Gn" (Gibbs cluster number) columns

## Output

The program produces an Excel file with multiple sheets:
- **Single_AA_variants**: Peptides with single amino acid variations from known proteins
- **Matches_to_six_frame**: Peptides matching regions in the six-frame translated human genome
- **Six_frame_non_matched**: Peptides not matching in the six-frame translated genome
- **cis_PCPS**: Potential proteasome-catalyzed spliced peptide variants (cis splicing)

Additional working files are created in the working directory.

## Pipeline Process

1. Filters peptides from Gibbs clustering based on MHC class and cluster numbers
2. Searches filtered peptides against known MHC peptide database
3. Searches non-matched peptides against human canonical protein database
4. Checks for single amino acid variants
5. Searches remaining peptides in six-frame translated human genome
6. Identifies potential proteasome-catalyzed peptide splicing events

## Logs

The pipeline generates detailed logs in the output directory (pipeline.log).

## Troubleshooting

- Ensure all database files are properly formatted with makeblastdb
- Check that the input Excel file has the correct sheet structure
- Verify file paths in command arguments (use absolute paths if possible)
- Check pipeline.log for detailed information about each step

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added]

## License

[License information to be added]
