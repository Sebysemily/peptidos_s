# PEPTIDOS Toxicity Workflow

Snakemake workflow for peptide preprocessing and toxicity prediction.

## Workflow

Clone the repository with external resources:

```bash
git clone --recursive https://github.com/<user>/<repo>.git
```

If the repository was already cloned, initialize resources with:

```bash
git submodule update --init --recursive
```

Run the pipeline with:

```bash
snakemake --cores all --use-conda
```

Current stages:

1. Cluster curated peptide FASTA files with MMseqs2.
2. Run ToxinPred3 on representative sequences.
3. Run ToxTeller on representative sequences.
4. Run CAPTP on compatible representative sequences.
5. Run HemoPI2 classification and regression on representative sequences.

## Inputs

Curated FASTA inputs are configured in `config/config.yml`:

```yaml
curated_fastas:
  "25_50": "data/curated_md-lais/25_50/alls_RangeSelected_25_50_curado.fasta"
  # "10_25": "data/curated_md-lais/10_25/alls_RangeSelected_10_25_curado.fasta"
```

Each key under `curated_fastas` is used as the workflow identifier and output
folder name. Add or uncomment entries here to request additional subsets.

These FASTA files were obtained from MD-Lais GUI FASTA curator using automatic
assignment and replacing each selected position with the most frequent residue.

## Configuration

Main configuration file:

```text
config/config.yml
```

Important fields:

- `curated_fastas`: FASTA files to process. The keys define the subset names.
- `batching`: number of shared toxin-check batches per subset. Missing subsets
  default to 1 batch.
- `max_threads`: global thread limit used by all threaded rules.
- `mmseqs`: MMseqs2 clustering parameters.
- `toxteller.program_dir`: local checkout path for ToxTeller.
- `captp.program_dir`: local checkout path for CAPTP.

Example:

```yaml
max_threads: 18

batching:
  "25_50": 1
  # "10_25": 1

toxteller:
  program_dir: "resources/ToxTeller/"

captp:
  program_dir: "resources/CAPTP/"
```

ToxTeller is stored inside the project under `resources/ToxTeller/`. If
`program_resource/toxteller.py` is not found there, initialize the resource
submodules with:

```bash
git submodule update --init --recursive
```

The checkout must contain `program_resource/toxteller.py` and the
ToxTeller model/scaler pickle files.

CAPTP is stored inside the project under `resources/CAPTP/`. The checkout must
contain `main.py`, `data/AAindex.pkl`, and `data/model_saved.pkl`.

ToxTeller and CAPTP are expected to be git submodules pinned by the parent
repository. Snakemake validates that their expected files exist, but it does not
clone them during normal workflow execution.

## Environments

Snakemake creates conda environments from:

- `envs/pre_processing.yml` for MMseqs2.
- `envs/tox_check.yml` for ToxinPred3, ToxTeller, and CAPTP.
- `envs/hemolisis_check.yml` for HemoPI2.

`envs/tox_check.yml` pins `toxinpred3==1.4` and includes PyTorch for CAPTP.
`envs/hemolisis_check.yml` installs HemoPI2 from pip.

## Outputs

Batch FASTA files:

```text
results/batches/tox_check/{peptide_set}/batch_{batch_id}.fasta
results/batches/toxteller/{peptide_set}/batch_{batch_id}.fasta
results/batches/captp/{peptide_set}/batch_{batch_id}.fasta
```

MMseqs2 representative sequences:

```text
data/curated_md-lais/mmseqs2/{peptide_set}/clusters_{peptide_set}_rep_seq.fasta
```

ToxinPred3 report:

```text
results/toxinpred3/{peptide_set}/clusters_{peptide_set}_rep_seq_toxinpred3.csv
```

ToxTeller report:

```text
results/toxteller/{peptide_set}/clusters_{peptide_set}_rep_seq_toxteller.csv
```

CAPTP report:

```text
results/captp/{peptide_set}/clusters_{peptide_set}_rep_seq_captp.csv
```

HemoPI2 reports:

```text
results/hemopi2_classification/{peptide_set}/clusters_{peptide_set}_rep_seq_hemopi2_classification.csv
results/hemopi2_regression/{peptide_set}/clusters_{peptide_set}_rep_seq_hemopi2_regression.csv
```

Known completed output:

```text
results/toxinpred3/25_50/clusters_25_50_rep_seq_toxinpred3.csv
```

## Notes

- ToxinPred3 and CAPTP use shared batches controlled by `batching`, then each
  tool merges its own CSVs.
- ToxTeller uses its own batches of up to 9,500 sequences because it refuses
  inputs above 10,000 sequences.
- CAPTP only receives sequences up to 49 aa because the bundled CAPTP
  preprocessing adds a `[CLS]` token and fails on 50-aa peptides. Empty
  sequences and longer peptides are omitted from CAPTP outputs and can be left
  blank in a combined final report.
- HemoPI2 is installed from pip because the local standalone checkout does not
  include the large model directory. Classification uses Hybrid1 RF+MERCI
  (`-m 2`) and regression reports HC50. HemoPI2 rules live in
  `rules/hem_check.smk`; future HemoPI2 helper code should live under
  `code/hem_check/`.
- CSV batch merging uses the shared helper `code/merge_csv_reports.py`.
- External ToxTeller/CAPTP versions are fixed by the submodule commits recorded
  in the parent repository. To restore them, run
  `git submodule update --init --recursive`.
- `10_25` is much larger than `25_50`; review runtime and storage before it is
  enabled as a default target.
- ToxinPred3 may create temporary files such as `seq.aac`, `seq.dpc`, and
  `Sequence_1`; these are workflow artifacts and should not be committed.
- If a previous Snakemake run was interrupted, unlock the working directory:

```bash
snakemake --unlock
```
