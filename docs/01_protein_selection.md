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

## 3. Virulence factor screening (VFDB)

Proteins were screened for similarity to experimentally validated virulence factors using the VFDB core dataset.

### 3.1. Requirements

- **BLAST+**: The standalone NCBI BLAST+ suite must be installed. It can be installed via conda:

  ```bash
  conda install -c bioconda blast
  ```
Or downloaded directly from the site:
https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html

- **VFDB Core Dataset**:
    1. Go to the VFDB download page: https://www.mgc.ac.cn/VFs/download.htm
    2. Download the "VFDB core dataser (setA)" in FASTA format.
    3. Decompress the file (e.g., `gunzip VFDB_setA_pro.fas.gz`) and place it in the working directory as `VFDB_setA_pro.fas`.

### 3.2. Execution

The pipeline uses two scripts located in `scripts/01_candidate_selection/`:

1. `05_vfdb_blast.sh`- Creates the BLAST database (if it doesn't exist) and runs `blastp`for all the strains against the VFDB core dataset.

    ```bash
    bash scripts/01_candidate_selection/05_vfdb_blast.sh
    ```
**Note**: This script assumes the query files (`secuencias_proteoma_<strain>_TODAS_localizaciones.txt`) are in the same directory where the script is executed.

2. `06_vfdb_analyze.py`- Filters the BLASTp results to retain only significant hits (E-value ≤ 1e-15 and identity ≥ 35%) and extracts the unique protein IDs.

    ```bash
    python scripts/01_candidate_selection/06_vfdb_analyze.py \
        -i VFDB_proteoma_referencia_resultados_BLASTp.txt \
        -p VFDB_proteoma_referencia_proteinas_virulentas_e-15.txt \
        -e 1e-15 \
        -p 35
    ```

## 4. Essentiality analysis (DEG)

Proteins that passed the VFDB screening were further filtered to identify those encoded by essential genes, using the Database of Essential Genes (DEG).

### 4.1. Requirements

- **DEG Dataset**:
  1. Go to the DEG download page: https://tubic.org/deg/public/index.php/download
  2. Download the **Bacteria -> amino acid** dataset (e.g., `DEG10.aa`).
  3. Place it in the working directory.

### 4.2. Execution

The workflow uses three scripts located in `scripts/01_candidate_selection/`:

1. `07_deg_prepare_sequences.py` - Extracts the FASTA sequences of the previously identified virulent proteins from the original proteome, creating the query files for BLAST.
   
   ```bash
   python scripts/01_candidate_selection/07_deg_prepare_sequences.py \
       -f Proteoma_referencia.faa \
       -i VFDB_proteoma_referencia_proteinas_virulentas_e-15.txt \
       -o DEG_proteoma_referencia_secuencias_virulentas_2.txt
    ```
2. `08_deg_blast.sh` - Automates the creation of the DEG BLAST database and runs `blastp`for the three strains.

    ```bash
    bash scripts/01_candidate_selection/08_deg_blast.sh
    ```
3. `09_deg_analyze.py`- Filters the BLASTp results (identity ≥ 35%, E-value < 1e-15) and extracts the unique essential protein IDs.

    ```bash
    python scripts/01_candidate_selection/09_deg_analyze.py \
        -i DEG_proteoma_referencia_resultado_BLASTp.txt \
        -o DEG_proteoma_referencia_proteinas_esenciales_2.txt
    ```

## 5. Host homology filtering

To minimize the risk of cross-reactivity, the essential proteins identified in the previous step were screened against the host proteome (*Salmo salar*).

### 5.1. Requirements

- **NCBI BLASTp (Web Server)**: This step is performed online using the NCBI BLAST service.
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch

### 5.2. Execution

The workflow requires two manual steps on the NCBI server followed by a local script execution:

1. **Online BLASTp submission**:
   1. Open the NCBI BLASTp web page.
   2. Paste the FASTA sequences of the essential proteins (or upload the FASTA file) into the query box.
   3. Set the database to **`nr` (Non-redundant protein sequences)**.
   4. In the "Organism" field, type and select **`Salmo salar (taxid:8030)`**.
   5. Click **BLAST** and wait for the results.
   6. Once finished, click **"Download All"** and save the results in **Text** format. 
   7. Rename the downloaded file to match the pipeline convention: `Homologia_proteoma_referencia_resultado_BLASTp_nr_2.txt`.
     Repeat this process for the `H2` and `DJ2R` strains.

2. **Local filtering script**: 
   Run the `10_host_homology_analyze.py` script to parse the downloaded NCBI text file and identify proteins with no significant homology to the host.

    ```
    python scripts/01_candidate_selection/10_host_homology_analyze.py \
       -b Homologia_proteoma_referencia_resultado_BLASTp_nr_2.txt \
       -i DEG_proteoma_referencia_proteinas_esenciales_2.txt \
       -o Homologia_proteoma_referencia_SinHomologia_nr_2.txt
    ```

### 5.3. Filtering rule

The script evaluates the BLAST results based on the following custom logic to ensure only non-homologous proteins are retained:

 - Keep if the protein has 0 hits (NCBI found no matches in the host proteome).
 - Keep if the protein has hits, but ALL of them have an identity < 25%, AND at least one hit has an E-value < 1e-15 (confirming a significant alignment).
 - Discard if any hit has an identity ≥ 25% (indicating potential cross-reactivity with the host).

## 6. Transmembrane domain prediction (TMHMM)

Proteins with zero or one transmembrane (TM) domain are ideal for vaccine development. Proteins with multiple TM domains were excluded.

### 6.1. Sequence extraction

### 6.1. Sequence extraction

First, extract the FASTA sequences of the proteins that passed the host homology filtering.

```
python scripts/01_candidate_selection/11_tmhmm_prepare_sequences.py \
    -f Proteoma_referencia.faa \
    -i Homologia_proteoma_referencia_SinHomologia_nr_2.txt \
    -o TMHMM_proteoma_referencia_secuencias_SinHomologia_2.txt
```

### 6.2. TMHMM web server execution

The analysis is performed on the DTU Health Tech TMHMM server.

1. Go to the TMHMM-2.0 server: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0
2. Upload the generated FASTA file.
3. In the Output format section, select "One line per protein".
4. Click Submit.
5. Once the results are ready, copy the text output and save it as `TMHMM_proteoma_referencia_resultado_2.txt`.

### 6.3. Filtering results

The script `12_tmhmm_analyze.py` parses the TMHMM output and retains proteins with 0 transmembrane domains (PredHel=0), which represent the optimal candidates.

    ```
    python scripts/01_candidate_selection/12_tmhmm_analyze.py \
        -i TMHMM_proteoma_referencia_resultado_2.txt \
        -o0 TMHMM_proteoma_referencia_0TM_secretadas_2.txt \
        -o1 TMHMM_proteoma_referencia_0-1TM_accesibles_2.txt
    ```
