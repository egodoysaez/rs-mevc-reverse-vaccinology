# Epitope prediction (B-cell, MHC-I, MHC-II)

This phase identifies the epitopes used to build the RS-MEVC multiepitope vaccine. It starts from the 15 consensus antigenic proteins generated in Phase 1.

## Important note on salmonid MHC prediction

NetMHCpan-4.1b and NetMHCIIpan-4.3i were trained on mammalian MHC data. When applied to *Salmo salar* MHC sequences, two limitations were observed:

1. **13 of the 48 Sasa-UBA alleles (MHC-I) could not be processed** by NetMHCpan. The pseudo-sequence generated from the salmon allele contained gaps (e.g. `YYV-YRQIGGANFANVAR...`) because the alignment to the human templates fails at those positions. Only **35 alleles** produced usable predictions.
2. For the 35 alleles that did run, the **"Distance to training data"** was >0.4 (ideal <0.1) and **`%Rank_EL` was `-99.9`** for every peptide. The Rank-based classifier is therefore not reliable for salmon MHC.
3. **Consequence**: The selection criterion was changed from `%Rank_EL` to **`Score_EL`**, which is computed independently of the training-data distance. All downstream MHC-I analysis uses `Score_EL`.

These predictions should be treated as **in silico prioritization**, not experimental evidence of peptide–MHC binding.

---

## 1. Requirements

| Tool | Version | Use |
|------|---------|-----|
| IPD-MHC | — | Allele sequences for *Salmo salar* |
| NetMHCpan | 4.1b | MHC-I prediction |

### 1.1 Allele download (IPD-MHC)

From https://www.ebi.ac.uk/ipd/mhc/group/FISH/:

- **Sasa-UBA** (MHC-I) → `Download → Protein FASTA` (48 alleles) → rename to `Sasa_UBA.txt`.

Place `Sasa_UBA.txt` in the working directory.

### 1.2 NetMHCpan installation (WSL on Windows)

NetMHCpan 4.1b was installed locally inside WSL (Ubuntu) on Windows:

```
tar -xzvf netMHCpan-4.1b.Linux.tar.gz -C ~/netMHCpan-4.1
export NETMHCpan="$HOME/netMHCpan-4.1"
export TMPDIR="$NETMHCpan/tmp"
```

The scripts call the binary via `subprocess.run` and pass `-rdir "$NETMHCpan"` so that NetMHCpan finds its data directory even when it is not installed system-wide.

## 2. MHC class I epitope prediction (NetMHCpan 4.1b)

### 2.1 Allele preparation

`01_pipeline_netmhcpan.py` splits `Sasa_UBA.txt` (48 sequences) into individual clean FASTA files (`<temp-dir>/allele_01.fsa`, …, `allele_48.fsa`), removing gaps (`-`), stop codons (`*`) and spaces.

### 2.2 NetMHCpan execution

For each allele, the script runs:

```
netMHCpan -rdir <NMHOME> -hlaseq allele_XX.fsa -BA -l 9 <proteins.fasta>
```

**Key flags**:
`-BA`: binding affinity mode.
`-l 9`: 9-mer peptides.
`-rdir <NMHOME>`: pointer to the NetMHCpan data directory.

```
python scripts/02_epitope_prediction/mhc_i/01_pipeline_netmhcpan.py \
    --netmhcpan-bin "$NETMHCpan/bin/netMHCpan" \
    --nmhome       "$NETMHCpan" \
    --proteins     Prediccion_epitopos_TOP15_antigenicas_secuencias_2.txt \
    --alleles      Sasa_UBA.txt \
    --output-dir   resultados_salmon_TOP15_antigenicas_2 \
    --temp-dir     temp_alelos_TOP15_antigenicas_2
```

Each successful allele produces one file `<output-dir>/allele_XX.txt` containing the full NetMHCpan detailed table.

### 2.3 Per-protein top-40 `Score_EL`extraction

NetMHCpan outputs predictions for every protein × every peptide, so the number of rows per protein is very different. For a fair comparison, we extract the 40 peptides with the highest Score_EL for each protein from every allele output.

The script preserves the NetMHCpan format (header, dashed lines, protein footer) and appends `_TOP40_BY_PROTEIN.txt` to the output name.

```
python scripts/02_epitope_prediction/mhc_i/02_sort_top40_scoreEL.py \
    -i resultados_salmon_TOP15_antigenicas_2 \
    -o resultados_salmon_TOP15_antigenicas_2
```

Or on a single file:

```
python scripts/02_epitope_prediction/mhc_i/02_sort_top40_scoreEL.py \
    -i resultados_salmon_TOP15_antigenicas_2/allele_*.txt
```

### 2.4 Promiscuity and allele coverage

The consolidator reads every `allele_*_TOP40_BY_PROTEIN.txt` file and, for each unique peptide (per protein), computes:

- Number of alleles in which the peptide appears with `Score_EL ≥ 0.60` (`--score-min`, default 0.60).
- Percentage of alleles covered (relative to `--num-alleles`, default 48).
- Mean / max / min / std of `Score_EL` across alleles.
- The list of alleles in which the peptide passed the threshold.

```
python scripts/02_epitope_prediction/mhc_i/03_consolidate_promiscuity.py \
    -i resultados_salmon_TOP15_antigenicas_2 \
    -o resultados_salmon_TOP15_antigenicas_2
```

## 3. MHC class II epitope prediction (NetMHCIIpan 4.3i)

### 3.1 Allele download and reformatting

From https://www.ebi.ac.uk/ipd/mhc/group/FISH/:

- **Sasa-DAA** (α-chain, 22 alleles) → `Sasa-DAA.txt`.
- **Sasa-DAB** (β-chain, 42 alleles) → `Sasa-DAB.txt`.

Original headers have the form `FISH08212|Sasa-DAA*13:01|225`. NetMHCIIpan needs just the genotype (`Sasa-DAA*13:01`). `04_reformat_alleles.py` extracts the genotype and writes:

- `Sasa-DAA_reformatted.txt` (22 sequences)
- `Sasa-DAB_reformatted.txt` (42 sequences)

It can also optionally clean a proteins FASTA by passing `--proteins`.

```
python scripts/02_epitope_prediction/mhc_ii/04_reformat_alleles.py \
    --daa Sasa-DAA.txt \
    --dab Sasa-DAB.txt \
    --proteins Prediccion_epitopos_TOP15_antigenicas_secuencias_2.txt
```

### 3.2 Prediction pipeline

`05_pipeline_netmhciipan.py` runs every protein x DAA x DAB combination:

```
netMHCIIpan-4.3 -f Proteina_X.txt -hlaseqA DAA_Y.txt -hlaseq DAB_Z.txt
```

```
source ~/netmhciipan_setup.sh   # sets NETMHCIIpan, TMPDIR and PATH
python scripts/02_epitope_prediction/mhc_ii/05_pipeline_netmhciipan.py \
    --netmhciipan-bin "$NETMHCIIpan/bin/netMHCIIpan-4.3" \
    --proteins Prediccion_epitopos_TOP15_antigenicas_secuencias_limpio_2.txt \
    --daa      Sasa-DAA_reformatted.txt \
    --dab      Sasa-DAB_reformatted.txt \
    --output-dir Sasa_prediccion_MHCII_TOP15_antigenicas_2
```

At the end the script also writes two preliminary summaries into `<output-dir>`:

- `1_epitopes_ranked_GLOBAL.txt` - peptides sorted by promiscuity across DAA×DAB pairs, filtered to ≥ 3 pairs.
- `1_epitopes_promiscuity_analysis_GLOBAL.txt` - top-10 peptides, per-peptide coverage of DAA×DAB pairs.

### 3.3 Promiscuity analysis

`06_epitope_promiscuity.py` is the analysis that produces the ranking used in the paper.

For each combination `prediction_P{prot}_{daa}_{dab}.txt`:

1. Read the peptides and their `Score_EL`.
2. Sort by `Score_EL` descending and keep the TOP 15 (`--top-n`, default 15).
3. Register each peptide in a global map `{peptide: {combinations: set, scores: list}}`.
4. Keep peptides that appear in ≥ `MIN_COMBINATIONS` combinations (`--min-combinations`, default 1).
5. Sort by `(num_combinations DESC, avg_score_el DESC)`.

```
python scripts/02_epitope_prediction/mhc_ii/06_epitope_promiscuity.py \
    --input-dir Sasa_prediccion_MHCII_TOP15_antigenicas_2 \
    --top-n 15 \
    --min-combinations 1
```

