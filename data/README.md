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
- DEG 10: https://tubic.org/deg/public/index.php/download
- IPD MHC (fish): https://www.ebi.ac.uk/ipd/mhc/group/FISH/
- *Salmo salar* proteome: GCF_905237065.1