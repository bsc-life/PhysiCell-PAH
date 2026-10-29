#!/usr/bin/env python3
"""
build_vessel_cross_section.py
==============================

Procedurally build a 2D pulmonary-artery cross-section initial condition
(concentric Intima / Media / Adventitia layers around a lumen) and write it
as a PhysiCell-format ``cells.csv`` -- for use BEFORE a real Xenium slice is
available, or as a clean idealized reference condition to compare a real
slice against.

Self-contained: only depends on ``numpy`` and ``pandas``.

--------------------------------------------------------------------------
ANATOMY ENCODED (2D cross-section, z = 0 -- a Xenium slice is a thin
optical section, so a 2D cross-section is the natural match)
--------------------------------------------------------------------------
                 (LUMEN -- empty, no agents)
                       |
                       v
      Intima   : Endothelial_Cell   monolayer lining the lumen
                       |
                       v
      Media    : Smooth_Muscle_Cell concentric rings (+ a minority
                 Synthetic_SMC subpopulation under --condition pah)
                       |
                       v
      Adventitia: Fibroblast        outer, more loosely packed layer

All cell-type names are plain, consistent English identifiers (no mixed
Spanish/English labels): ``Endothelial_Cell``, ``Smooth_Muscle_Cell``,
``Synthetic_SMC``, ``Fibroblast``, optionally ``Pericyte``.

Each layer is filled with concentric rings of cells (not random packing):
within a layer spanning radius [r_inner, r_outer], N radial sub-rings are
placed (N chosen from the layer thickness and the cell diameter), each
sub-ring holding as many cells as fit its circumference at the requested
packing fraction, with alternating angular offset between sub-rings
(brick-like offset) and a small random jitter for a natural look. This
directly mirrors real small-pulmonary-artery histology -- concentric SMC
rings in the media, a single endothelial monolayer at the lumen.

--------------------------------------------------------------------------
CONTROL vs PAH-REMODELED GEOMETRY
--------------------------------------------------------------------------
``--condition pah`` thickens the media and raises its SMC packing fraction
(hyperplasia/hypertrophy of the remodeled wall), and seeds a minority
Synthetic_SMC subpopulation within the media -- both choices are
simple, documented approximations set from the proportions in this
project's own cell-type inventory report (control vs PAH SMC counts;
STEAP4hi_SMC-like fraction of Xenium SMCs), not fabricated single-cell
data. Override any of these with explicit flags; see --help.

--------------------------------------------------------------------------
EXAMPLE USAGE
--------------------------------------------------------------------------
# Default control-like vessel cross-section:
python build_vessel_cross_section.py --output cells.csv

# PAH-remodeled vessel (thicker media, denser SMCs, STEAP4hi minority):
python build_vessel_cross_section.py --condition pah --output cells_pah.csv

# Fully custom geometry:
python build_vessel_cross_section.py \\
    --lumen-radius 120 --media-thickness 90 --adventitia-thickness 70 \\
    --smc-diameter 15 --media-packing-fraction 0.8 \\
    --include-pericytes --output cells.csv --dry-run
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

_TAB20_HEX = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5",
]


def sanitize_name(label: str) -> str:
    """Turn an arbitrary cell-type label into a PhysiCell-safe identifier."""
    s = re.sub(r"[^0-9A-Za-z_]+", "_", str(label).strip())
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        s = "unknown_type"
    if s[0].isdigit():
        s = f"type_{s}"
    return s


def fill_annulus_with_rings(
    r_inner: float,
    r_outer: float,
    cell_diameter: float,
    packing_fraction: float,
    rng: np.random.Generator,
    jitter_frac: float = 0.15,
) -> np.ndarray:
    """Fill the annulus [r_inner, r_outer) with concentric rings of cells.

    Returns an (N, 2) array of (x, y) coordinates centered at the origin.
    ``packing_fraction`` in (0, 1]: 1.0 packs cells edge-to-edge along each
    ring's circumference; smaller values space them further apart (lower
    density, e.g. for loosely-packed adventitial connective tissue).
    """
    thickness = r_outer - r_inner
    if thickness <= 0:
        return np.empty((0, 2))

    n_rings = max(1, round(thickness / cell_diameter))
    points: List[Tuple[float, float]] = []
    for i in range(n_rings):
        r = r_inner + (i + 0.5) * thickness / n_rings
        circumference = 2 * np.pi * r
        spacing = cell_diameter / max(packing_fraction, 1e-6)
        n_in_ring = max(1, int(np.floor(circumference / spacing)))
        angle_offset = (np.pi / n_in_ring) if (i % 2 == 1) else 0.0
        thetas = np.linspace(0, 2 * np.pi, n_in_ring, endpoint=False) + angle_offset
        thetas = thetas + rng.uniform(-jitter_frac, jitter_frac, size=n_in_ring) * (2 * np.pi / n_in_ring)
        radii = r + rng.uniform(-jitter_frac, jitter_frac, size=n_in_ring) * (thickness / n_rings)
        xs = radii * np.cos(thetas)
        ys = radii * np.sin(thetas)
        points.extend(zip(xs, ys))
    return np.array(points) if points else np.empty((0, 2))


def build_vessel(args: argparse.Namespace) -> pd.DataFrame:
    rng = np.random.default_rng(args.seed)
    rows = []

    r_lumen = args.lumen_radius
    r_intima_outer = r_lumen + args.ec_diameter
    r_media_outer = r_intima_outer + args.media_thickness
    r_adventitia_outer = r_media_outer + args.adventitia_thickness

    # Intima: a single-cell-thick endothelial monolayer lining the lumen.
    ec_xy = fill_annulus_with_rings(
        r_lumen, r_intima_outer, args.ec_diameter, args.ec_packing_fraction, rng
    )
    rows.append(pd.DataFrame({"x": ec_xy[:, 0], "y": ec_xy[:, 1], "cell_type": "Endothelial_Cell"}))

    media_inner = r_intima_outer
    if args.include_pericytes:
        pericyte_thickness = args.pericyte_diameter * 1.0
        pc_xy = fill_annulus_with_rings(
            media_inner, media_inner + pericyte_thickness, args.pericyte_diameter,
            args.pericyte_packing_fraction, rng,
        )
        rows.append(pd.DataFrame({"x": pc_xy[:, 0], "y": pc_xy[:, 1], "cell_type": "Pericyte"}))
        media_inner = media_inner + pericyte_thickness

    # Media: concentric SMC rings, with an optional Synthetic_SMC minority.
    smc_xy = fill_annulus_with_rings(
        media_inner, media_inner + args.media_thickness, args.smc_diameter,
        args.media_packing_fraction, rng,
    )
    n_smc = len(smc_xy)
    is_steap4 = rng.uniform(size=n_smc) < args.steap4_fraction
    smc_types = np.where(is_steap4, "Synthetic_SMC", "Smooth_Muscle_Cell")
    rows.append(pd.DataFrame({"x": smc_xy[:, 0], "y": smc_xy[:, 1], "cell_type": smc_types}))

    # Adventitia: fibroblasts, more loosely packed.
    adv_inner = media_inner + args.media_thickness
    fib_xy = fill_annulus_with_rings(
        adv_inner, adv_inner + args.adventitia_thickness, args.fibroblast_diameter,
        args.adventitia_packing_fraction, rng,
    )
    rows.append(pd.DataFrame({"x": fib_xy[:, 0], "y": fib_xy[:, 1], "cell_type": "Fibroblast"}))

    df = pd.concat(rows, ignore_index=True)
    df["z"] = args.z_value
    return df


def build_legend(df: pd.DataFrame) -> pd.DataFrame:
    counts = df["cell_type"].value_counts()
    rows = []
    for type_id, (label, n) in enumerate(counts.items()):
        rows.append(
            {
                "type_id": type_id,
                "cell_type_name": sanitize_name(label),
                "original_label": label,
                "n_cells": int(n),
                "suggested_color_hex": _TAB20_HEX[type_id % len(_TAB20_HEX)],
            }
        )
    return pd.DataFrame(rows)


def write_outputs(df: pd.DataFrame, legend: pd.DataFrame, output: Path, legend_output: Path, type_format: str) -> None:
    label_to_id = dict(zip(legend["original_label"], legend["type_id"]))
    label_to_name = dict(zip(legend["original_label"], legend["cell_type_name"]))
    out = pd.DataFrame({"x": df["x"].round(4), "y": df["y"].round(4), "z": df["z"].round(4)})
    out["type"] = df["cell_type"].map(label_to_name if type_format == "name" else label_to_id)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    legend_output.parent.mkdir(parents=True, exist_ok=True)
    legend.to_csv(legend_output, index=False)


def apply_condition_defaults(args: argparse.Namespace) -> None:
    """Fill in --condition-dependent defaults for any flag the user did not
    explicitly set (argparse default sentinel: None)."""
    if args.condition == "control":
        defaults = dict(media_thickness=60.0, media_packing_fraction=0.72, steap4_fraction=0.0)
    else:  # "pah"
        # Geometry scaling approximated from this project's own inventory
        # report (control -> PAH): SMC counts rose ~3x in snRNA-seq and the
        # Xenium SMC population expanded further in the remodeled wall;
        # STEAP4hi_SMC-like cells were ~794/13,023 (~6%) of Xenium PAH SMCs.
        defaults = dict(media_thickness=170.0, media_packing_fraction=0.85, steap4_fraction=0.06)
    for key, val in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, val)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a 2D pulmonary-artery cross-section initial condition as a PhysiCell cells.csv.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--condition", choices=["control", "pah"], default="control",
                    help="Preset geometry: healthy vs PAH-remodeled wall (thicker/denser media + STEAP4hi minority).")

    p.add_argument("--lumen-radius", type=float, default=100.0, help="Lumen radius, µm.")
    p.add_argument("--media-thickness", type=float, default=None, help="Media thickness, µm (default depends on --condition).")
    p.add_argument("--adventitia-thickness", type=float, default=90.0, help="Adventitia thickness, µm.")

    p.add_argument("--ec-diameter", type=float, default=12.0, help="Endothelial cell diameter, µm.")
    p.add_argument("--smc-diameter", type=float, default=15.0, help="Smooth muscle cell diameter, µm.")
    p.add_argument("--fibroblast-diameter", type=float, default=15.0, help="Fibroblast diameter, µm.")
    p.add_argument("--pericyte-diameter", type=float, default=10.0, help="Pericyte diameter, µm (if --include-pericytes).")

    p.add_argument("--ec-packing-fraction", type=float, default=0.95, help="Intima packing fraction (0-1].")
    p.add_argument("--media-packing-fraction", type=float, default=None, help="Media packing fraction (default depends on --condition).")
    p.add_argument("--adventitia-packing-fraction", type=float, default=0.45, help="Adventitia packing fraction (0-1].")
    p.add_argument("--pericyte-packing-fraction", type=float, default=0.5, help="Pericyte-layer packing fraction (if --include-pericytes).")

    p.add_argument("--include-pericytes", action="store_true", help="Insert a thin Pericyte layer between intima and media.")
    p.add_argument("--steap4-fraction", type=float, default=None,
                    help="Fraction of media cells labeled Synthetic_SMC instead of Smooth_Muscle_Cell (default depends on --condition).")

    p.add_argument("--z-value", type=float, default=0.0, help="Constant z coordinate to write (default: 0.0).")
    p.add_argument("--type-format", choices=["name", "id"], default="name", help="Write 'type' as sanitized name or integer ID.")
    p.add_argument("--seed", type=int, default=0, help="Random seed for jitter/offsets (reproducible layout).")

    p.add_argument("--output", type=Path, default=Path("cells.csv"), help="Output cells.csv path.")
    p.add_argument("--legend-output", type=Path, default=None, help="Output legend CSV path (default: <output>_legend.csv).")
    p.add_argument("--dry-run", action="store_true", help="Print summary only; do not write files.")

    args = p.parse_args(argv)
    apply_condition_defaults(args)
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    df = build_vessel(args)
    legend = build_legend(df)

    r_media_outer = args.lumen_radius + args.ec_diameter + (args.pericyte_diameter if args.include_pericytes else 0) + args.media_thickness
    r_adventitia_outer = r_media_outer + args.adventitia_thickness

    print(f"Condition: {args.condition}")
    print(f"Lumen radius: {args.lumen_radius:.1f} µm")
    print(f"Media outer radius: {r_media_outer:.1f} µm  (thickness {args.media_thickness:.1f} µm, packing {args.media_packing_fraction:.2f})")
    print(f"Adventitia outer radius: {r_adventitia_outer:.1f} µm")
    print(f"\nTotal cells: {len(df):,}")
    for _, row in legend.iterrows():
        print(f"  {row.type_id:>2d}  {row.cell_type_name!r:25s} n={row.n_cells:,}")
    print(
        f"\nSuggested PhysiCell domain (copy into PhysiCell_settings.xml):\n"
        f"  x: [{-r_adventitia_outer:.1f}, {r_adventitia_outer:.1f}]\n"
        f"  y: [{-r_adventitia_outer:.1f}, {r_adventitia_outer:.1f}]\n"
        f"  2D simulation (z held at {args.z_value})"
    )

    if args.dry_run:
        print("\n--dry-run set: no files written.")
        return 0

    legend_output = args.legend_output or args.output.with_name(args.output.stem + "_legend.csv")
    write_outputs(df, legend, args.output, legend_output, args.type_format)
    print(f"\nWrote {len(df):,} cells to {args.output}")
    print(f"Wrote legend to {legend_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
