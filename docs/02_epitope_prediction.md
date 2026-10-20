# Epitope prediction

This phase identifies the epitopes used to build the RS-MEVC multiepitope vaccine. It starts from the 15 consensus antigenic proteins generated in Phase 1.

## Important note on salmonid MHC prediction

NetMHCpan-4.1b and NetMHCIIpan-4.3i were trained on mammalian MHC data. When applied to *Salmo salar* MHC sequences, two limitations were observed:

1. **13 of the 48 Sasa-UBA alleles (MHC-I) could not be processed** by NetMHCpan. The pseudo-sequence generated from the salmon allele contained gaps (e.g. `YYV-YRQIGGANFANVAR...`) because the alignment to the human templates fails at those positions. Only **35 alleles** produced usable predictions.
2. For the 35 alleles that did run, the **"Distance to training data"** was >0.4 (ideal <0.1) and **`%Rank_EL` was `-99.9`** for every peptide. The Rank-based classifier is therefore not reliable for salmon MHC.
3. **Consequence**: The selection criterion was changed from `%Rank_EL` to **`Score_EL`**, which is computed independently of the training-data distance. All downstream MHC analysis uses `Score_EL`.

These predictions should be treated as **in silico prioritization**, not experimental evidence of peptide-MHC binding.

---

## 1. Requirements

| Tool | Version | Use |
|------|---------|-----|
| IPD-MHC | - | Allele sequences for *Salmo salar* |
| NetMHCpan | 4.1b | MHC-I prediction |
| NetMHCIIpan | 4.3i | MHC-II prediction |

### 1.1 Allele download (IPD-MHC)

From https://www.ebi.ac.uk/ipd/mhc/group/FISH/:

- **Sasa-UBA** (MHC-I) → `Download → Protein FASTA` (48 alleles) → rename to `Sasa_UBA.txt`.
- **Sasa-DAA** (MHC-II α-chain, 22 alleles) → `Sasa-DAA.txt`.
- **Sasa-DAB** (MHC-II β-chain, 42 alleles) → `Sasa-DAB.txt`.

Place these files in the working directory. They are **not tracked** in the repository.

### 1.2 NetMHCpan / NetMHCIIpan installation (WSL on Windows)

Both tools were installed locally inside WSL (Ubuntu). The pipeline scripts receive the binary and data paths via CLI arguments, so no shell-wide setup is strictly required. Templates for a per-session setup are provided in `scripts/setup/`:

```bash
# Optional: source a setup file before running the pipelines
cp scripts/setup/netmhcpan_setup.sh.example   ~/netmhcpan_setup.sh
cp scripts/setup/netmhciipan_setup.sh.example ~/netmhciipan_setup.sh
# edit the *_HOME variables, then:
source ~/netmhcpan_setup.sh
source ~/netmhciipan_setup.sh
```

## 2. MHC class I epitope prediction (NetMHCpan 4.1b)

### 2.1 Allele preparation

`01_pipeline_netmhcpan.py` splits `Sasa_UBA.txt` (48 sequences) into individual clean FASTA files (`<temp-dir>/allele_01.fsa`, ..., `allele_48.fsa`), removing gaps (`-`), stop codons (`*`) and spaces.

### 2.2 NetMHCpan execution

For each allele, the script runs:

```
netMHCpan -rdir <NMHOME> -hlaseq allele_XX.fsa -BA -l 9 <proteins.fasta>
```

**Key flags**:
- `-BA`: binding affinity mode.
- `-l 9`: 9-mer peptides.
- `-rdir <NMHOME>`: pointer to the NetMHCpan data directory.

```bash
python scripts/02_epitope_prediction/01_pipeline_netmhcpan.py \
    --netmhcpan-bin "$NETMHCpan/bin/netMHCpan" \
    --nmhome       "$NETMHCpan" \
    --proteins     Prediccion_epitopos_TOP15_antigenicas_secuencias_2.txt \
    --alleles      Sasa_UBA.txt \
    --output-dir   resultados_salmon_TOP15_antigenicas_2 \
    --temp-dir     temp_alelos_TOP15_antigenicas_2
```

Each successful allele produces one file `<output-dir>/allele_XX.txt` with the full NetMHCpan detailed table. The 13 alleles that fail (pseudo-sequence gaps) produce `<output-dir>/allele_XX.txt.err`.

### 2.3 Per-protein top-40 `Score_EL` extraction

NetMHCpan outputs predictions for every protein x every peptide, so the number of rows per protein is very different. For a fair comparison, we extract the 40 peptides with the highest `Score_EL` for each protein from every allele output.

The script preserves the NetMHCpan format (header, dashed lines, protein footer) and appends `_TOP40_BY_PROTEIN.txt` to the output name. It accepts either a directory (all `allele_*.txt` files inside) or a single file.

```bash
# Process an entire directory
python scripts/02_epitope_prediction/02_sort_top40_scoreEL.py \
    -i resultados_salmon_TOP15_antigenicas_2 \
    -o resultados_salmon_TOP15_antigenicas_2

# Or a single allele file
python scripts/02_epitope_prediction/02_sort_top40_scoreEL.py \
    -i resultados_salmon_TOP15_antigenicas_2/allele_01.txt
```

> **Note.** Do not pass `allele_*.txt` to `-i`: the shell expands the glob to multiple arguments and `argparse` rejects them. Use a directory or a single file.

### 2.4 Promiscuity and allele coverage

The consolidator reads every `allele_*_TOP40_BY_PROTEIN.txt` file and, for each unique peptide (per protein), computes:

- Number of alleles in which the peptide appears with `Score_EL ≥ 0.60` (`--score-min`, default 0.60).
- Percentage of alleles covered (relative to `--num-alleles`, default 48 - see note below).
- Mean / max / min / std of `Score_EL` across alleles.
- The list of alleles in which the peptide passed the threshold.

```bash
python scripts/02_epitope_prediction/03_consolidate_promiscuity.py \
    -i resultados_salmon_TOP15_antigenicas_2 \
    -o resultados_salmon_TOP15_antigenicas_2
```

## 3. MHC class II epitope prediction (NetMHCIIpan 4.3i)

### 3.1 Allele download and reformatting

Original IPD-MHC headers have the form `FISH08212|Sasa-DAA*13:01|225`. NetMHCIIpan needs just the genotype (`Sasa-DAA*13:01`). `04_reformat_alleles.py` extracts the genotype and writes:

- `Sasa-DAA_reformatted.txt` (22 sequences)
- `Sasa-DAB_reformatted.txt` (42 sequences)

It can also optionally clean a proteins FASTA (reducing each header to the `WP_` accession) via `--proteins`.

```bash
python scripts/02_epitope_prediction/04_reformat_alleles.py \
    --daa      Sasa-DAA.txt \
    --dab      Sasa-DAB.txt \
    --proteins Prediccion_epitopos_TOP15_antigenicas_secuencias_2.txt
```

Output files (in the working directory):

- `Sasa-DAA_reformatted.txt`
- `Sasa-DAB_reformatted.txt`
- `Prediccion_epitopos_TOP15_antigenicas_secuencias_limpio_2.txt`

### 3.2 Prediction pipeline

`05_pipeline_netmhciipan.py` runs every (protein x DAA x DAB) combination:

```
netMHCIIpan-4.3 -f Proteina_X.txt -hlaseqA DAA_Y.txt -hlaseq DAB_Z.txt
```

- `15 x 22 x 42 = 13 860` predictions.
- Parallelised with `multiprocessing.Pool` (`cpu_count() - 1` workers by default).
- One output file per combination: `<output-dir>/prediction_P{prot}_{daa}_{dab}.txt`.
- Per-prediction timeout: 300 s.

The script takes the NetMHCIIpan binary as an absolute path (`--netmhciipan-bin`), so no shell-wide setup is required. Only `TMPDIR` is recommended:

```bash
# Option A: inline (no setup script needed)
export NETMHCIIPAN_HOME="/path/to/netMHCIIpan-4.3"
export TMPDIR="$NETMHCIIPAN_HOME/tmp"
mkdir -p "$TMPDIR"

python scripts/02_epitope_prediction/05_pipeline_netmhciipan.py \
    --netmhciipan-bin "$NETMHCIIPAN_HOME/bin/netMHCIIpan-4.3" \
    --proteins Prediccion_epitopos_TOP15_antigenicas_secuencias_limpio_2.txt \
    --daa      Sasa-DAA_reformatted.txt \
    --dab      Sasa-DAB_reformatted.txt \
    --output-dir Sasa_prediccion_MHCII_TOP15_antigenicas_2

# Option B: using the template setup file (see scripts/setup/README.md)
# source ~/netmhciipan_setup.sh
```

At the end the script also writes two **preliminary** summaries into `<output-dir>`:

- `1_epitopes_ranked_GLOBAL.txt` - peptides sorted by promiscuity across DAAxDAB pairs, filtered to ≥ 3 pairs.
- `1_epitopes_promiscuity_analysis_GLOBAL.txt` - top-10 peptides with per-peptide coverage of DAAxDAB pairs.

### 3.3 Final promiscuity analysis

`06_epitope_promiscuity.py` is the analysis that produces the ranking used in the paper.

For each combination `prediction_P{prot}_{daa}_{dab}.txt`:

1. Read the peptides and their `Score_EL` (9th column of the NetMHCIIpan detailed table; value in `[0, 1]`, higher is better).
2. Sort by `Score_EL` descending and keep the TOP 15 (`--top-n`, default 15).
3. Register each peptide in a global map `{peptide: {combinations: set, scores: list}}`.
4. Keep peptides that appear in ≥ `MIN_COMBINATIONS` combinations (`--min-combinations`, default 1).
5. Sort by `(num_combinations DESC, avg_score_el DESC)`.

```bash
python scripts/02_epitope_prediction/06_epitope_promiscuity.py \
    --input-dir Sasa_prediccion_MHCII_TOP15_antigenicas_2 \
    --top-n 15 \
    --min-combinations 1
```

Output file: `<input-dir>/epitopes_promiscuity_el_analysis_v3.txt`.

## 4. Output files summary

### 4.1 MHC-I

| File | Description | Generated by |
|------|-------------|--------------|
| `<output-dir>/allele_XX.txt` | Raw NetMHCpan output per allele (35 files) | `01_pipeline_netmhcpan.py` |
| `<output-dir>/allele_XX.txt.err` | Error log for failed alleles (13 files) | `01_pipeline_netmhcpan.py` |
| `<output-dir>/allele_XX_TOP40_BY_PROTEIN.txt` | Top-40 peptides per protein per allele | `02_sort_top40_scoreEL.py` |
| `epitopos_promiscuos.txt` | Main MHC-I peptide ranking by promiscuity | `03_consolidate_promiscuity.py` |
| `cobertura_alélica.txt` | Allele coverage with the top-N peptides | `03_consolidate_promiscuity.py` |
| `epitopos_detallado.txt` | Full peptide x allele `Score_EL` matrix | `03_consolidate_promiscuity.py` |

### 4.2 MHC-II

| File | Description | Generated by |
|------|-------------|--------------|
| `Sasa-DAA_reformatted.txt` | 22 DAA sequences with clean headers | `04_reformat_alleles.py` |
| `Sasa-DAB_reformatted.txt` | 42 DAB sequences with clean headers | `04_reformat_alleles.py` |
| `<output-dir>/prediction_P{prot}_{daa}_{dab}.txt` | Raw NetMHCIIpan output per combination (~13 860 files) | `05_pipeline_netmhciipan.py` |
| `<output-dir>/1_epitopes_ranked_GLOBAL.txt` | Preliminary global ranking (≥ 3 DAAxDAB pairs) | `05_pipeline_netmhciipan.py` |
| `<output-dir>/1_epitopes_promiscuity_analysis_GLOBAL.txt` | Preliminary coverage summary | `05_pipeline_netmhciipan.py` |
| `<output-dir>/epitopes_promiscuity_el_analysis_v3.txt` | **Final** MHC-II peptide ranking (TOP-15 per combination) | `06_epitope_promiscuity.py` |

