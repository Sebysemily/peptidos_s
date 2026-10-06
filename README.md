# PEPTIDOS Predictive Workflow
<!--toc:start-->
- [Requirements](#requirements)
- [Useful Commands](#useful-commands)
- [Current Flow](#current-flow)
- [Main Outputs](#main-outputs)
- [Global Parameters](#global-parameters)
- [Tools and Core Parameters](#tools-and-core-parameters)
  - [Pre-processing](#pre-processing)
  - [Toxicity Prediction](#toxicity-prediction)
  - [Hemolysis Prediction](#hemolysis-prediction)
  - [Immunogenicity Prediction](#immunogenicity-prediction)
  - [Anticancer Prediction (ACP)](#anticancer-prediction-acp)
<!--toc:end-->

Snakemake workflow for peptide preprocessing and multi-model prediction. The pipeline systematically evaluates candidate peptides across 4 predictive modules (Toxicity, Hemolysis, Immunogenicity, and Anticancer properties) to identify safe and effective therapeutic candidates.

Requirements
----------

- `snakemake` with `conda` or `mamba`.
- Internet connection to install the environments.
- Clone the repository with external resources using `git clone --recursive <URL>` or run `git submodule update --init --recursive` if already cloned.

Useful Commands
---------------

```bash
# Complete main flow
snakemake --cores all --use-conda
```

Current Flow
------------

The pipeline executes 4 predictive modules iteratively, filtering and characterizing peptide candidates:

1. **Inputs & Pre-processing:** Curated FASTA files configured in `config/config.yml` are clustered with MMseqs2.
2. **Toxicity Prediction (Module 1):** Representative sequences are evaluated using ToxinPred3, ToxTeller, and CAPTP. A combined toxicity summary is generated, filtering out toxic candidates.
3. **Hemolysis Prediction (Module 2):** Non-toxic peptides are evaluated using HemoPI2, Macrel, HEPAD, and Hemo_DL. Generates a hemolytic summary and a non-hemolytic sequence subset.
4. **Immunogenicity Prediction (Module 3):** Non-hemolytic peptides are screened using Algpred2, AllergenAI, and AllerTrans. Generates an immunogenicity summary.
5. **Anticancer Prediction (Module 4):** Evaluates sequences for anti-cancer properties (ACP) using AntiCP2.
6. **Properties & Final Report:** Final physical/chemical characteristics and overall filtering reports (`metadata/characteristics.csv`, `metadata/filtering.csv`) are produced.

Main Outputs
-------------------

- Pre-processed sequences: `data/curated_md-lais/mmseqs2/`
- Filtered FASTA subsets: `data/derived/non_toxic/` and `data/derived/non_hemo/`
- Combined Summaries:
  - Toxicity: `results/tox_check/toxicity_summary/`
  - Hemolysis: `results/hemo_check/hemolytic_summary/`
  - Immunogenicity: `results/inmuno_check/inmuno_summary/`
- Final physical/chemical properties: `metadata/characteristics.csv`
- Final general filtering report: `metadata/filtering.csv`

Global Parameters
-----------------

Files: `config/config.yml`.

To avoid cluttering the configuration files, not all parameters used across the pipeline are exposed as global variables. The main parameters in the config file are:

- **Inputs:** `curated_fastas` maps the peptide set name (e.g., `"25_50"`) to its file path.
- **Resources:** `max_threads: 18`.
- **Batching:** `batching` controls the number of shared check batches per subset to manage memory limits.
- **Tool Settings:** MMseqs2 clustering parameters and global workflow behaviors (e.g., `process_groups_together`).

Tools and Core Parameters
-------------------------

The pipeline relies on several bioinformatics and machine learning tools. Tools are isolated in specific conda environments (e.g., `envs/tox_check/toxinpred3_captp.yml`, `envs/hem_check/hemopi2.yml`).

### Pre-processing

- **MMseqs2:** Clusters curated peptide FASTA files to extract representative sequences.

### Toxicity Prediction

- **ToxinPred3:** ToxinPred3 and CAPTP use shared batches controlled by the `batching` config. Temporary workflow artifacts (like `seq.aac`) are cleaned up.
- **ToxTeller:** Uses its own batches of up to 9,500 sequences because it refuses inputs above 10,000 sequences. Requires older model-compatible scikit-learn stack.
- **CAPTP:** Only receives sequences up to 49 amino acids because its preprocessing adds a `[CLS]` token and fails on 50-aa peptides. Longer sequences are automatically omitted.

### Hemolysis Prediction

- **HemoPI2:** Classification uses Hybrid1 RF+MERCI (`-m 2`) and regression reports HC50. Installed via pip as the local standalone does not include the large model directory.
- **Macrel:** Runs in a separate environment due to Bioconda dependency conflicts with the modern Python/PyTorch stack.
- **HEPAD:** Hemolytic activity predictor.
- **Hemo_DL:** Deep learning based hemolysis prediction.

### Immunogenicity Prediction

- **Algpred2:** Predicts allergenic peptides.
- **AllergenAI:** Assesses sequence allergenicity.
- **AllerTrans:** Transformer-based prediction model for allergens.

### Anticancer Prediction (ACP)

- **AntiCP2:** Predicts Anticancer Peptides based on sequence features.
