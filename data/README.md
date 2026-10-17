# Data

This folder is for input data. Large files are not tracked by Git.

## Required files

Download manually from NCBI Datasets:

1. Go to https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1646
2. Download the following assemblies (Protein (FASTA); RefSeq/GCF):

|   Strain   |     RefSeq      |            Save as             |
|------------|-----------------|--------------------------------|
| ATCC 33209 | GCF_000018885.1 | `data/Proteoma_referencia.faa` |
|     H2     | GCF_003343285.1 |     `data/Proteoma_H2.faa`     |
|    DJ2R    | GCF_003343265.1 |    `data/Proteoma_DJ2R.faa`    |

Place files in this folder before running the scripts.

## Other databases

- VFDB core dataset: http://www.mgc.ac.cn/VFs/
- DEG 15: https://tubic.org/deg/public/index.php/download
- IPD MHC (fish): https://www.ebi.ac.uk/ipd/mhc/group/FISH/
- *Salmo salar* proteome: GCF_905237065.1# Data

This folder should contain input data for the pipeline. Large files are not tracked by Git (see `.gitignore`). Download them manually:

- `GCF_000018885.1_protein.faa` — NCBI RefSeq, ATCC 33209 reference proteome
- `GCF_003343265.1_protein.faa` — NCBI RefSeq, DJ2R field strain
- `GCF_003343285.1_protein.faa` — NCBI RefSeq, H2 field strain
- `vfdb.fasta` — Virulence Factor Database
- `deg15.fasta` — Database of Essential Genes
- `sasa_mhc_i.fasta` — *Salmo salar* MHC-I sequences (IPD)
- `sasa_mhc_ii.fasta` — *Salmo salar* MHC-II sequences (IPD)

Place files in this folder before running the scripts.