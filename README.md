# RS-MEVC reverse vaccinology pipeline

Scripts and workflows used for the design of a multi-epitope vaccine candidate against *Renibacterium salmoninarum*, the etiological agent of bacterial kidney disease (BKD) in salmonids.

## Associated publication

Morales C, Arriagada V, Bonilla E, Prina B, **Godoy E**, Molina-Chavez B, Osorio A, Gutierrez E, Roa V, Gomez LA, Oñate AA (2026). *Reverse vaccinology and immunoinformatics-driven design of a multi-epitope vaccine candidate against Renibacterium salmoninarum*. Frontiers in Bioinformatics 6:1899416. doi: [10.3389/fbinf.2026.1899416](https://doi.org/10.3389/fbinf.2026.1899416)

## My contribution

I was responsible for:
- Retrieval and curation of the three *R. salmoninarum* proteomes (ATCC 33209, DJ2R, H2).
- Subcellular localization filtering (PSORTb v3.0).
- Virulence factor screening (VFDB).
- Essentiality analysis (DEG 15).
- Host homology exclusion (BLASTp against *Salmo salar*).
- Transmembrane domain prediction (TMHMM v2.0).
- Antigenicity assessment (Vaxijen v2.0).
- B-cell epitope prediction (BepiPred-3.0).
- MHC-I and MHC-II epitope prediction (NetMHCpan-4.1, NetMHCIIpan-4.0) using *S. salar* MHC sequences.
- Multi-epitope construct design (RS-MEVC) with flagellin D0/D1 adjuvant and linkers.
- Physicochemical characterization (ProtParam, AllerTOP, ToxinPred, SOLpro).

Molecular docking, molecular dynamics simulations, and MM-PBSA calculations were performed by a co-author and are not included in this repository.

## Pipeline overview

See `workflow.png` and `docs/` for details.

1. Proteome retrieval from NCBI.
2. Subcellular localization filtering.
3. Virulence factor screening.
4. Essentiality analysis.
5. Host homology exclusion.
6. Transmembrane domain filtering.
7. Antigenicity assessment.
8. Consensus protein selection across three strains.
9. B-cell, MHC-I, and MHC-II epitope prediction.
10. Multi-epitope construct design.

## Repository structure
.
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── workflow.png
├── docs/
├── scripts/
├── notebooks/
├── data/
└── results/


## Requirements

- Python 3.10+
- R 4.3+ (optional, for downstream analyses)
- Biopython, pandas, numpy, matplotlib
- See `environment.yml` for the full list.

## Usage

```bash
conda env create -f environment.yml
conda activate rs-mevc
bash scripts/01_download_proteomes.sh
python scripts/02_filter_localization.py
python scripts/03_predict_epitopes.py
python scripts/04_assemble_construct.py

## Data

All input data are publicly available:

- NCBI RefSeq: GCF_000018885.1, GCF_003343265.1, GCF_003343285.1
- VFDB: http://www.mgc.ac.cn/VFs/
- DEG 15: https://tubic.org/deg/public/index.php
- IPD MHC: https://www.ebi.ac.uk/ipd/mhc/

Scripts expect data in `data/` with specific filenames. See `data/README.md`.

## Limitations

- MHC-I and MHC-II predictions were generated using NetMHCpan-4.1 and NetMHCIIpan-4.0, whose models are trained predominantly on mammalian data. Predictions for salmonid MHC molecules should be interpreted as an in silico prioritization step, not as experimental evidence of peptide-MHC binding.

- The Score_EL ≥ 0.6 threshold is not a salmonid-specific validated cutoff.

- B-cell epitope predictions are candidate antigenic regions, not experimentally validated epitopes.

- No docking, MD, or MM-PBSA analyses are included in this repository.

## License

MIT License. See `LICENSE`.

## Contact

Esteban Godoy Sáez

- LinkedIn: https://www.linkedin.com/in/estebangodoysaez
- ORCID: https://orcid.org/0009-0009-6854-3064
- Email: egodoy2016@udec.cl