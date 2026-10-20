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
```text
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
```

## Requirements

- Python 3.10+
- R 4.3+ (optional, for downstream analyses)
- Biopython, pandas, numpy, matplotlib
- See `environment.yml` for the full list.

## Usage

## Usage

The full pipeline is documented step-by-step in `docs/01_protein_selection.md`. 
All scripts for the candidate selection phase are located in 
`scripts/01_candidate_selection/`.

A typical execution for one strain (e.g., the reference) would be:

```bash
# 1. Subcellular localization post-processing
python scripts/01_candidate_selection/01_psortb_count.py -i PSORTb_proteoma_referencia.txt -o ...
python scripts/01_candidate_selection/02_psortb_split.py -i PSORTb_proteoma_referencia.txt -s referencia
python scripts/01_candidate_selection/03_psortb_extract_sequences.py -f data/Proteoma_referencia.faa -s referencia
python scripts/01_candidate_selection/04_psortb_merge_localizations.py -s referencia -o ...

# 2. Virulence factor screening (VFDB)
bash scripts/01_candidate_selection/05_vfdb_blast.sh
python scripts/01_candidate_selection/06_vfdb_analyze.py -i VFDB_proteoma_referencia_resultados_BLASTp.txt ...

# 3. Essentiality analysis (DEG)
python scripts/01_candidate_selection/07_deg_prepare_sequences.py -f data/Proteoma_referencia.faa ...
bash scripts/01_candidate_selection/08_deg_blast.sh
python scripts/01_candidate_selection/09_deg_analyze.py -i DEG_proteoma_referencia_resultado_BLASTp.txt ...

# 4. Host homology filtering
python scripts/01_candidate_selection/10_host_homology_analyze.py -b Homologia_proteoma_referencia_resultado_BLASTp_nr_2.txt ...

# 5. Transmembrane domain prediction (TMHMM)
python scripts/01_candidate_selection/11_tmhmm_prepare_sequences.py -f data/Proteoma_referencia.faa ...
python scripts/01_candidate_selection/12_tmhmm_analyze.py -i TMHMM_proteoma_referencia_resultado_2.txt ...

# 6. Antigenicity prediction (VaxiJen)
python scripts/01_candidate_selection/13_vaxijen_prepare_sequences.py -f data/Proteoma_referencia.faa ...
python scripts/01_candidate_selection/14_vaxijen_analyze.py -i Antigenicidad_proteoma_referencia_resultado_VaxiJen_2.txt ...

# 7. Consensus selection
python scripts/01_candidate_selection/15_consensus_compare.py --ref ... --dj2r ... --h2 ... -o ...

# 8. Top 15 antigenic proteins
python scripts/01_candidate_selection/16_get_top15.py -c Proteoma_proteinas_candidatas_R_salmoninarum_2.txt ...
```

**Note**: Steps that require web servers (PSORTb, TMHMM, VaxiJen, BLASTp on NCBI) must be performed manually as described in `docs/01_protein_selection.md` before running the post-processing scripts.

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