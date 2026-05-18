# PEPTIDOS Toxicity Workflow

Snakemake workflow for peptide preprocessing and toxicity prediction.

## Workflow

Run the pipeline with:

```bash
snakemake --cores all --use-conda
```

Current stages:

1. Cluster curated peptide FASTA files with MMseqs2.
2. Run ToxinPred3 on representative sequences.
3. Run ToxTeller on representative sequences.

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
- `batching`: number of ToxinPred3 batches per subset. Missing subsets default
  to 1 batch.
- `max_threads`: global thread limit used by all threaded rules.
- `mmseqs`: MMseqs2 clustering parameters.
- `toxteller.program_dir`: local checkout path for ToxTeller.

Example:

```yaml
max_threads: 18

batching:
  "25_50": 1
  # "10_25": 1

toxteller:
  program_dir: "resources/ToxTeller/"
```

ToxTeller is stored inside the project under `resources/ToxTeller/`. If
`program_resource/toxteller.py` is not found there, the workflow clones the
public ToxTeller repository.

The cloned repository must contain `program_resource/toxteller.py` and the
ToxTeller model/scaler pickle files.

## Environments

Snakemake creates conda environments from:

- `envs/pre_processing.yml` for MMseqs2.
- `envs/tox_check.yml` for ToxinPred3 and ToxTeller.

`envs/tox_check.yml` pins `toxinpred3==1.4`.

## Outputs

Batch FASTA files:

```text
results/batches/toxinpred3/{peptide_set}/batch_{batch_id}.fasta
results/batches/toxteller/{peptide_set}/batch_{batch_id}.fasta
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

Known completed output:

```text
results/toxinpred3/25_50/clusters_25_50_rep_seq_toxinpred3.csv
```

## Notes

- ToxinPred3 batches are controlled manually by `batching`.
- ToxTeller batches are split automatically into chunks of 10,000 sequences.
- `10_25` is much larger than `25_50` and should use multiple batches before it
  is enabled as a default target.
- ToxTeller has an internal sequence limit, so the default batch size is 10,000
  sequences per call.
- ToxinPred3 may create temporary files such as `seq.aac`, `seq.dpc`, and
  `Sequence_1`; these are workflow artifacts and should not be committed.
- If a previous Snakemake run was interrupted, unlock the working directory:

```bash
snakemake --unlock
```
