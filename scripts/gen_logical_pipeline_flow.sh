#!/usr/bin/env bash
set -euo pipefail

OUTDIR=${1:-analysis}
mkdir -p "$OUTDIR"

DOT="$OUTDIR/pipeline_logical_flow.dot"
SVG="$OUTDIR/pipeline_logical_flow.svg"

python3 - "$DOT" <<'PY'
from pathlib import Path
import sys

out = Path(sys.argv[1])
out.write_text(
    r'''digraph pipeline_logical_flow {
  graph [
    rankdir = LR,
    bgcolor = "white",
    pad = 0.35,
    nodesep = 0.55,
    ranksep = 0.8,
    splines = ortho,
    fontname = "DejaVu Sans"
  ];

  node [
    shape = box,
    style = "rounded,filled",
    fontname = "DejaVu Sans",
    fontsize = 11,
    margin = "0.12,0.08",
    color = "#4b5563",
    fillcolor = "#f8fafc"
  ];

  edge [
    color = "#64748b",
    arrowsize = 0.7,
    fontname = "DejaVu Sans",
    fontsize = 9
  ];

  curated [
    label = "Curated peptide FASTA\n10_25 and 25_50",
    fillcolor = "#eef2ff"
  ];

  mmseqs [
    label = "MMseqs2 clustering\nrepresentative sequences",
    fillcolor = "#eef2ff"
  ];

  external_resources [
    label = "Check external resources\nToxTeller + CAPTP",
    fillcolor = "#f1f5f9"
  ];

  tox_batches [
    label = "Indexed toxicity batches\nshared IDs + mapping CSV",
    fillcolor = "#ecfeff"
  ];

  toxteller_batches [
    label = "ToxTeller batches\nmax 9,500 seq/batch",
    fillcolor = "#ecfeff"
  ];

  toxinpred3 [
    label = "ToxinPred3\nbatch predictions",
    fillcolor = "#fff7ed"
  ];

  toxteller [
    label = "ToxTeller\nLR/SVM/RF/XGBoost",
    fillcolor = "#fff7ed"
  ];

  captp_filter [
    label = "CAPTP input filter\nlength <= 49 aa",
    fillcolor = "#fefce8"
  ];

  captp [
    label = "CAPTP\nbatch predictions",
    fillcolor = "#fff7ed"
  ];

  tox_merge [
    label = "Merge toxicity reports\nper tool + peptide set",
    fillcolor = "#fffbeb"
  ];

  tox_summary [
    label = "Toxicity summary\npass = no tool predicts toxic",
    fillcolor = "#dcfce7"
  ];

  non_toxic_fasta [
    label = "Non-toxic representative FASTA\nkeeps original peptide IDs",
    fillcolor = "#dcfce7"
  ];

  non_toxic_batches [
    label = "Indexed non-toxic batches\nnew internal IDs + mapping CSV",
    fillcolor = "#ccfbf1"
  ];

  hemopi2_class [
    label = "HemoPI2 classification\nHemolytic / Non-Hemolytic",
    fillcolor = "#f3e8ff"
  ];

  hemopi2_reg [
    label = "HemoPI2 regression\nHC50 prediction",
    fillcolor = "#f3e8ff"
  ];

  macrel [
    label = "Macrel peptides\nAMP + hemolysis prediction",
    fillcolor = "#f3e8ff"
  ];

  hemo_merge [
    label = "Merge hemolysis reports\nper tool + peptide set",
    fillcolor = "#ede9fe"
  ];

  final_outputs [
    label = "Final outputs\nToxicity summaries +\nnon-toxic FASTA +\nhemolysis reports",
    fillcolor = "#e0f2fe"
  ];

  curated -> mmseqs;
  mmseqs -> tox_batches;
  mmseqs -> toxteller_batches;

  tox_batches -> toxinpred3;
  toxteller_batches -> toxteller;
  tox_batches -> captp_filter;
  captp_filter -> captp;

  external_resources -> toxteller [style = dashed];
  external_resources -> captp [style = dashed];

  toxinpred3 -> tox_merge;
  toxteller -> tox_merge;
  captp -> tox_merge;

  mmseqs -> tox_summary;
  tox_merge -> tox_summary;

  tox_summary -> non_toxic_fasta;
  mmseqs -> non_toxic_fasta;

  non_toxic_fasta -> non_toxic_batches;

  non_toxic_batches -> hemopi2_class;
  non_toxic_batches -> hemopi2_reg;
  non_toxic_batches -> macrel;

  hemopi2_class -> hemo_merge;
  hemopi2_reg -> hemo_merge;
  macrel -> hemo_merge;

  tox_summary -> final_outputs;
  non_toxic_fasta -> final_outputs;
  hemo_merge -> final_outputs;
}
''',
    encoding="utf-8",
)
PY

if command -v dot >/dev/null 2>&1; then
  dot -Tsvg "$DOT" -o "$SVG"
  echo "Wrote $DOT and $SVG"
else
  echo "Wrote $DOT"
  echo "Warning: 'dot' not found; SVG conversion skipped." >&2
fi
