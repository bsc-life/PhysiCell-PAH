#!/usr/bin/env python3
"""
personalize_intracellular_model.py
====================================

Derive cell-type- and condition-specific MaBoSS ``.cfg`` files for the PASMC
model (``models/ModeloPASMC_base_v2.bnd`` -- the team's final, SIGNOR/NeKo-
curated 34-node model; supersedes the earlier 11-node prototype
``pasmc_calcium_nfat_kv15.bnd``) from snRNA-seq pseudobulk marker-gene
expression, so each PhysiCell/PhysiBoSS cell type/condition in the digital
twin runs the intracellular model with a biologically appropriate starting
point rather than one fixed default.

Self-contained: depends only on ``anndata``, ``pandas``, ``numpy`` and
``pyyaml``.

--------------------------------------------------------------------------
WHAT GETS PERSONALIZED, AND WHY ONLY FOUR NODES
--------------------------------------------------------------------------
``ModeloPASMC_base_v2.bnd`` has nine "free" input nodes (``rate_up =
rate_down = 0`` -- state set only by ``.istate``, never by network logic),
but only FOUR are gene/channel-identity nodes a resting-state transcriptome
can set: ``CACNA1C``, ``STIM1``, ``ORAI1``, ``ATP2A2``. The other five
(``AMLO``, ``ET1``, ``HIPOXIA``, ``K_ALTO``, ``CPA``) are experimental/
pharmacological SCENARIO SWITCHES (drug exposure, stimulation, hypoxia,
depolarization assay, SERCA-inhibitor probe) set by the simulation's
experimental design, not by this script -- see ``gene_node_mapping.yaml``'s
``scenario_inputs`` section. Every other node in the model (Kv1_5, Em_depol,
Bkca, Ltcc, Plc, Serca, Sr_lleno, Stim1_pm, Soce, Trpc6, Pkc, Ncx_rev,
Ca_global, Ca_local, Ca_ltcc, Rock, Mlcp, Mlck, Calcineurina, Nfat, Hif1a,
CONTR_1/2, CREC_1/2) is derived by the network's own logic and is
deliberately NOT personalized independently -- see ``gene_node_mapping.yaml``
for the full rationale and the gene <-> node mapping (editable without
touching this script).

For each (cell_type, condition) group in the input data, this script:
  1. Computes pseudobulk mean expression of each mapped gene (mean across
     all cells in that group; expects ``adata.X`` -- or ``--layer`` -- to
     already be normalized/log-transformed, consistent with standard
     scanpy preprocessing).
  2. Min-max-scales each node's expression to [0, 1] across ALL groups in
     the file (so "SMCs in PAH" and "fibroblasts in Control" sit on the
     same scale), then clips to the ``floor``/``ceiling`` from the mapping
     YAML (keeps every node stochastically reachable -- never hard-forced).
  3. Writes that scaled value as a MaBoSS probabilistic initial state:
     ``[NODE].istate = (1-p) [0], p [1];`` -- i.e. the node starts ON with
     probability p, in every one of MaBoSS's Monte Carlo trajectories. Note
     the square brackets around the node name -- MaBoSS's probabilistic-
     istate grammar requires them (bare ``NODE.istate = p0 [0], p1 [1];``
     is a syntax error); confirmed against a real MaBoSS run, not just a
     string-level check.
  4. Copies all other settings (non-personalized istates, is_internal
     flags, simulation parameters) from ``--base-cfg`` unchanged, and
     writes the result to ``<out-dir>/pasmc_<cell_type>_<condition>.cfg``.

A ``personalization_summary.csv`` is also written with one row per
(cell_type, condition, node): raw pseudobulk expression, the genes
averaged, the scaled [0,1] value, and the final istate probability -- so
every generated ``.cfg`` is traceable back to the expression it came from.

--------------------------------------------------------------------------
EXPECTED INPUT SCHEMA
--------------------------------------------------------------------------
``adata.obs`` must contain a cell-type column (``--celltype-col``, e.g.
``"cell_type"``) and a condition column (``--condition-col``, e.g.
``"condition"``, values like ``"Control"``/``"PAH"``). ``adata.var_names``
must include the gene symbols named in ``gene_node_mapping.yaml``
(CACNA1C, STIM1, ORAI1, ATP2A2 by default) -- symbols are matched
case-sensitively; use ``--list-genes`` to check what's actually in the file
first.

--------------------------------------------------------------------------
EXAMPLE USAGE (once the real file is available)
--------------------------------------------------------------------------
# Confirm CACNA1C/STIM1/ORAI1/ATP2A2 (and your obs columns) are present:
python personalize_intracellular_model.py --input PAH_snRNAseq.h5ad --list-genes

# Generate one personalized .cfg per (cell_type, condition):
python personalize_intracellular_model.py \\
    --input PAH_snRNAseq.h5ad \\
    --celltype-col cell_type --condition-col condition \\
    --mapping gene_node_mapping.yaml \\
    --base-cfg ../models/ModeloPASMC_base_v2.cfg \\
    --out-dir ../results/personalized_cfgs \\
    --cell-types "Smooth muscle cells" "STEAP4hi SMC-like" "EC arterial"

No patient data is bundled with this script -- point --input at your own
snRNA-seq .h5ad. It has not been run against any fabricated/synthetic data.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import yaml

try:
    import anndata
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires the 'anndata' package (pip install anndata). "
        f"Import failed with: {exc}"
    )


def sanitize_name(label: str) -> str:
    s = re.sub(r"[^0-9A-Za-z_]+", "_", str(label).strip())
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unknown"


def load_mapping(path: Path) -> Dict[str, dict]:
    with open(path) as fh:
        spec = yaml.safe_load(fh)
    nodes = spec.get("nodes", {})
    if not nodes:
        raise SystemExit(f"No 'nodes' section found in mapping file {path}")
    return nodes


def list_genes(adata: "anndata.AnnData", mapping: Dict[str, dict]) -> None:
    all_genes = sorted({g for node in mapping.values() for g in node["genes"]})
    present = set(adata.var_names)
    print("Gene availability check:")
    for g in all_genes:
        print(f"  {g:15s} {'FOUND' if g in present else 'MISSING'}")
    print(f"\nadata.obs columns: {list(adata.obs.columns)}")


def pseudobulk_expression(
    adata: "anndata.AnnData",
    celltype_col: str,
    condition_col: str,
    genes: List[str],
    layer: str = None,
) -> pd.DataFrame:
    """Return a long DataFrame: cell_type, condition, gene, mean_expression."""
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise SystemExit(
            f"Gene(s) not found in adata.var_names: {missing}. "
            "Run with --list-genes to see what is available."
        )
    sub = adata[:, genes]
    X = sub.layers[layer] if layer else sub.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)

    df = pd.DataFrame(X, columns=genes, index=adata.obs_names)
    df["cell_type"] = adata.obs[celltype_col].astype(str).to_numpy()
    df["condition"] = adata.obs[condition_col].astype(str).to_numpy()

    long = df.melt(id_vars=["cell_type", "condition"], var_name="gene", value_name="expression")
    pb = long.groupby(["cell_type", "condition", "gene"], as_index=False)["expression"].mean()
    return pb


def build_node_signal(
    pb: pd.DataFrame,
    node: str,
    node_spec: dict,
) -> pd.DataFrame:
    """Average the mapped genes per (cell_type, condition), min-max scale
    across all groups, clip to [floor, ceiling]."""
    genes = node_spec["genes"]
    floor = float(node_spec.get("floor", 0.02))
    ceiling = float(node_spec.get("ceiling", 0.98))

    node_pb = pb[pb["gene"].isin(genes)]
    signal = node_pb.groupby(["cell_type", "condition"], as_index=False)["expression"].mean()
    signal = signal.rename(columns={"expression": "raw_expression"})

    lo, hi = signal["raw_expression"].min(), signal["raw_expression"].max()
    if hi > lo:
        scaled = (signal["raw_expression"] - lo) / (hi - lo)
    else:
        scaled = pd.Series(0.5, index=signal.index)  # no variation -> neutral midpoint
    signal["scaled_expression"] = scaled.clip(0.0, 1.0)
    signal["istate_prob"] = (floor + signal["scaled_expression"] * (ceiling - floor)).clip(floor, ceiling)
    signal["node"] = node
    signal["genes_used"] = ", ".join(genes)
    return signal


def write_personalized_cfg(
    base_cfg_text: str,
    node_probs: Dict[str, float],
    out_path: Path,
) -> None:
    """Replace each `NODE.istate = ...;` line (for personalized nodes only)
    in the base cfg with a probabilistic istate line; leave every other
    line untouched.

    MaBoSS's probabilistic-istate grammar REQUIRES square brackets around
    the node name (`[NODE].istate = p0 [0], p1 [1];`) -- a bare
    `NODE.istate = ...` with probability terms is a configuration syntax
    error. This was caught by actually running MaBoSS on a generated cfg,
    not by the string-level unit test alone -- see
    docs/amlodipine_adaptation_mechanism.md for the real-data analysis that
    surfaced it.
    """
    lines = base_cfg_text.splitlines()
    pattern = re.compile(r"^\s*(\w+)\.istate\s*=")
    out_lines = []
    replaced = set()
    for line in lines:
        m = pattern.match(line)
        if m and m.group(1) in node_probs:
            node = m.group(1)
            p = node_probs[node]
            out_lines.append(f"[{node}].istate = {1 - p:.4f} [0], {p:.4f} [1];")
            replaced.add(node)
        else:
            out_lines.append(line)

    missing = set(node_probs) - replaced
    if missing:
        out_lines.append("")
        out_lines.append("// personalize_intracellular_model.py: nodes not found in base cfg, appended")
        for node in sorted(missing):
            p = node_probs[node]
            out_lines.append(f"[{node}].istate = {1 - p:.4f} [0], {p:.4f} [1];")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out_lines) + "\n")


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Personalize the PASMC intracellular MaBoSS model per cell type/condition from snRNA-seq.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", required=True, type=Path, help="Path to the snRNA-seq .h5ad file.")
    p.add_argument("--mapping", type=Path,
                    default=Path(__file__).resolve().parent.parent / "configs" / "gene_node_mapping.yaml",
                    help="Path to gene_node_mapping.yaml (default: ../configs/gene_node_mapping.yaml next to this script).")
    p.add_argument("--base-cfg", type=Path, default=None, help="Base .cfg file to personalize (non-personalized settings copied verbatim).")
    p.add_argument("--layer", default=None, help="adata.layers key to use instead of adata.X (e.g. 'lognorm').")

    p.add_argument("--celltype-col", default="cell_type", help="adata.obs column with cell-type labels.")
    p.add_argument("--condition-col", default="condition", help="adata.obs column with condition labels (e.g. Control/PAH).")
    p.add_argument("--cell-types", nargs="*", default=None, help="Restrict to these cell-type labels (default: all).")
    p.add_argument("--conditions", nargs="*", default=None, help="Restrict to these condition labels (default: all).")

    p.add_argument("--out-dir", type=Path, default=Path("personalized_cfgs"), help="Output directory for personalized .cfg files.")
    p.add_argument("--summary-csv", type=Path, default=None, help="Output path for the traceability summary CSV (default: <out-dir>/personalization_summary.csv).")

    p.add_argument("--list-genes", action="store_true", help="Print gene availability / obs columns and exit.")
    p.add_argument("--dry-run", action="store_true", help="Compute and print the summary; do not write .cfg files.")

    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    if not args.input.exists():
        raise SystemExit(f"--input file not found: {args.input}")
    if not args.mapping.exists():
        raise SystemExit(f"--mapping file not found: {args.mapping}")

    mapping = load_mapping(args.mapping)
    adata = anndata.read_h5ad(args.input, backed="r" if args.layer is None else None)

    if args.list_genes:
        list_genes(adata, mapping)
        return 0

    for col in (args.celltype_col, args.condition_col):
        if col not in adata.obs.columns:
            raise SystemExit(f"Column {col!r} not found in adata.obs. Available: {list(adata.obs.columns)}")

    all_genes = sorted({g for node in mapping.values() for g in node["genes"]})
    pb = pseudobulk_expression(adata, args.celltype_col, args.condition_col, all_genes, layer=args.layer)

    if args.cell_types:
        pb = pb[pb["cell_type"].isin(args.cell_types)]
    if args.conditions:
        pb = pb[pb["condition"].isin(args.conditions)]
    if pb.empty:
        raise SystemExit("No (cell_type, condition) groups remain after filtering -- check --cell-types/--conditions.")

    signals = [build_node_signal(pb, node, spec) for node, spec in mapping.items()]
    summary = pd.concat(signals, ignore_index=True)
    summary = summary[["cell_type", "condition", "node", "genes_used", "raw_expression", "scaled_expression", "istate_prob"]]

    print(summary.to_string(index=False))

    if args.dry_run:
        print("\n--dry-run set: no files written.")
        return 0

    if args.base_cfg is None:
        raise SystemExit("--base-cfg is required unless --dry-run/--list-genes.")
    base_cfg_text = args.base_cfg.read_text()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    n_written = 0
    for (cell_type, condition), grp in summary.groupby(["cell_type", "condition"]):
        node_probs = dict(zip(grp["node"], grp["istate_prob"]))
        out_name = f"pasmc_{sanitize_name(cell_type)}_{sanitize_name(condition)}.cfg"
        write_personalized_cfg(base_cfg_text, node_probs, args.out_dir / out_name)
        n_written += 1
    print(f"\nWrote {n_written} personalized .cfg file(s) to {args.out_dir}")

    summary_csv = args.summary_csv or (args.out_dir / "personalization_summary.csv")
    summary.to_csv(summary_csv, index=False)
    print(f"Wrote traceability summary to {summary_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
