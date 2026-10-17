# Candidate protein selection

## 1. Proteome retrieval

Three *Renibacterium salmoninarum* proteomes were retrieved from NCBI Datasets:

|   Strain   |   Assembly  |     RefSeq      |         Role         |
|------------|-------------|-----------------|----------------------|
| ATCC 33209 |  ASM1888v1  | GCF_000018885.1 |      Reference       |
|     H2     | ASM334328v1 | GCF_003343285.1 | Chilean field strain |
|    DJ2R    | ASM334326v1 | GCF_003343265.1 | Chilean field strain |

Source: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1646

Download: for each genome, select **Protein (FASTA)**. RefSeq (GCF) versions were preferred over GenBank (GCA).

Files were renamed to:
- `Proteoma_referencia.faa`
- `Proteoma_H2.faa`
- `Proteoma_DJ2R.faa`

These files are **not** included in the repository. See `data/README.md` for download instructions.

## 2. Subcellular localization (PSORTb v3.0)

Subcellular localization was predicted for each proteome using the PSORTb v3.0 web server (https://www.psort.org/psortb/).

**Parameters:**
- Organism: Bacteria
- Gram stain: Gram-positive
- Output format: Long format (tab-delimited)

**Output files (not tracked in repository):**
- `PSORTb_proteoma_referencia.txt`
- `PSORTb_proteoma_H2.txt`
- `PSORTb_proteoma_DJ2R.txt`

**Post-processing scripts:**

1. `01_psortb_count.py` - counts proteins per localization category.
   ```bash
   python scripts/01_candidate_selection/01_psortb_count.py \
       -i PSORTb_proteoma_referencia.txt \
       -o PSORTb_proteoma_referencia_conteo_localizacion.txt
   ```
2. `02_psortb_split.py` - splits SeqIDs by localization into separate files.
    ```bash
    python scripts/01_candidate_selection/02_psortb_split.py \
        -i PSORTb_proteoma_referencia.txt \
        -s referencia
    ```
3. `03_psortb_extract_sequences.py` - retrieves FASTA sequences for each category from the original proteome.
    ```bash
    python scripts/01_candidate_selection/03_psortb_extract_sequences.py \
        -f Proteoma_referencia.faa \
        -s referencia
    ```
4. `04_psortb_merge_localizations.py`- merges all localizations into a single file.
    ```bash
    python scripts/01_candidate_selection/04_psortb_merge_localizations.py \
        -s referencia \
        -o secuencias_proteoma_referencia_TODAS_localizaciones.txt
    ```

**Filtering rule:**
Only proteins with known localization (extracellular, cell wall, cytoplasmic membrane, cytoplasm) were retained. Proteins with unknown localization were excluded.

