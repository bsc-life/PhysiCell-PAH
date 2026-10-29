#!/usr/bin/env python3
"""
generate_physiboss_vessel_gradient_config.py
================================================
Generates a buildable PhysiBoSS project that replaces the well-mixed MaBoSS
gradient sweep (docs/amlodipine_adaptation_mechanism.md section 9) with a
SPATIAL one: a real vessel cross-section (Intima/Media/Adventitia, from
build_vessel_cross_section.py) with 4 diffusible substrates (AMLO, ET1, O2,
K_ALTO), each SMC agent running its own instance of
ModeloPASMC_v8.bnd/.cfg, sensing its LOCAL substrate
concentration rather than one researcher-chosen combination shared by the
whole well-mixed population.

ANATOMICAL RATIONALE FOR EACH SUBSTRATE'S SOURCE
----------------------------------------------------
Sources are modeled as ENDOTHELIAL (Intima) secretion, not custom
Dirichlet-node C++ code -- simpler, and more biologically literal (the
endothelium genuinely is the source/gatekeeper for all three signals in
vivo):
  - AMLO:   EC secretes toward a target concentration representing the
            blood-side drug level; diffuses inward through Media/Adventitia
            as the wall takes it up -- a real, if simplified, proxy for
            systemic drug crossing the endothelium into the wall.
  - ET1:    EC secretes it (endothelium is ET1's classic source) AND, in
            the feedback-ON arm, so do activated SMCs themselves (Nfat ->
            "ET1 secretion" -- PASMC autocrine ET-1 production under basal
            conditions is independently well-documented: Fagan et al.,
            "Regulation of endothelin-1 synthesis in human pulmonary
            arterial smooth muscle cells," 1999; MacLean & Dempsie,
            "The role of endothelin-1 in the pathogenesis of pulmonary
            arterial hypertension," 2010 -- NOT the contested Johnson et
            al. 2020 claim; Nfat is used here as the model's general
            activation/stress-integration proxy, not a literal claim that
            NFAT specifically transcribes ET1).
  - O2:     EC secretes toward a physiological target (~38 mmHg, the
            PhysiCell default), diffusing inward and consumed by all cells
            -- the classic diffusion-limited-hypoxia pattern. In the PAH
            geometry, media thickness is ~170 micron (vs ~60 in control),
            so diffusion-limited hypoxia in the thickened media is itself
            a real, well-established PAH pathophysiology concept (medial
            thickening -> diffusion-limited medial hypoxia), independent
            of the feedback question. In the feedback-ON arm, SMCs ALSO
            consume MORE O2 when active (Ca_global -> "O2 uptake") -- a
            generic, uncontested cell-biology principle (metabolically
            active cells consume more oxygen), not a PASMC-specific claim.
  - K_ALTO: uniform initial condition only, no secretion/diffusion source
            -- the 40mM KCl depolarization test is a bath-applied
            experimental maneuver, not a physiological spatial gradient;
            modeling it as a true gradient would misrepresent how the
            experiment is actually done.

INPUT MAPPING ACTIONS (PhysiBoSS native, addons/PhysiBoSS/src/utils.h)
---------------------------------------------------------------------------
  AMLO, ET1, K_ALTO -> their Boolean nodes: action="activation" (node
    turns on above threshold) -- same mechanism validated in the
    well-mixed ecology study (section 7).
  O2 -> HIPOXIA: action="inhibition" -- confirmed in utils.h's
    updateNode(): with action="inhibition", the node activates when the
    substrate value drops BELOW threshold. This is exactly the low-O2-is-
    hypoxic relationship; no custom C++ needed.

OUTPUT MAPPING (feedback arm only)
---------------------------------------
  CREC_1 -> "cycle entry" (division rate, both arms, same as section 7).
  Nfat -> "ET1 secretion", Ca_global -> "O2 uptake" (feedback arm only) --
  both confirmed as native PhysiCell behavior-vocabulary strings
  ("<substrate> secretion"/"<substrate> uptake", core/
  PhysiCell_signal_behavior.cpp) addressable via the same native <output>
  mapping mechanism as "cycle entry" -- no custom C++ needed.

CAVEATS
-------
- Division/mechanics are minimized (near-zero division rate, no active
  motility) for this first test, to isolate the spatial-substrate-response
  question from population dynamics -- the explicit reason being that
  section 7's ecology study already covers growth/selection dynamics; this
  test is about whether POSITION changes a cell's intracellular attractor.
- The feedback rate parameters (secretion_rate under Nfat=1 / uptake_rate
  under Ca_global=1) are new, uncalibrated placeholders, following this
  project's established convention of flagging new rate constants
  explicitly rather than hiding them as if validated.
- Section 9's documented HIPOXIA-dependent sign-flip caveat for the
  TRPC6_feedback model (a possible ceiling effect, not independently
  confirmed) still applies here and may complicate interpretation of the
  O2 gradient's effect near the adventitia (relatively hypoxic zone).

USAGE
-----
    python generate_physiboss_vessel_gradient_config.py \
        --out-dir ../physicell_runs/vessel_gradient_feedback_on \
        --feedback on --max-time 4320
"""
import argparse
import math
import re
import shutil
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# Real snRNA-seq-fitted Beta_minmax distributions (PAH condition, n=3510 SMCs),
# from smc_heterogeneity_summary.csv (analyze_smc_expression_heterogeneity.py,
# AIC-selected best-fit family). CACNA1C is deliberately excluded: ModeloPASMC_v8
# made it Myocd-derived rather than a frozen free input, so personalizing it
# directly would fight the new contractile/synthetic toggle's own dynamics.
REAL_GENE_BETA_PARAMS_PAH = {
    "STIM1":  dict(alpha=0.09969838718917508, beta=3.033976979540983,  vmin=0.0, vmax=12.0),
    "ATP2A2": dict(alpha=0.11293322211791486, beta=2.7785016987965463, vmin=0.0, vmax=16.0),
    # ORAI1: NOT the raw snRNA-seq Beta_minmax fit (alpha=0.0959, beta=16.95,
    # vmax=9 -> ~1-3% ON) -- that real-data fit is contradicted by (a) the PAH
    # literature, which treats STIM1/Orai1-mediated SOCE as a broadly-active,
    # often-upregulated pathway in PASMC, not a rare-subpopulation phenomenon,
    # and (b) ORAI1's gene length (16.6 kb) being 14-44x shorter than
    # STIM1/ATP2A2/CACNA1C -- single-NUCLEUS RNA-seq relies heavily on
    # intronic reads for capture, so a short-intron gene is structurally
    # disadvantaged for detection regardless of true expression (see docs/
    # amlodipine_adaptation_mechanism.md for the full writeup). This keeps
    # genuine per-agent heterogeneity (NOT deterministic ORAI1=1 for everyone
    # either) but with a literature-informed, dropout-corrected baseline: a
    # smooth unimodal Beta(3,2) shape (no boundary pile-up) solved to give
    # ~65-69% of agents ON at the shared 0.5 raw-count threshold -- a
    # deliberately uncalibrated placeholder (no paper reports "fraction of
    # PASMC with functional SOCE"), chosen only to sit in the literature's
    # "broadly active" regime instead of the dropout-suppressed ~1-3%.
    "ORAI1":  dict(alpha=3.0, beta=2.0, vmin=0.0, vmax=1.0),
    # TRPC6: raw snRNA-seq Beta_minmax fit, used WITHOUT a dropout correction
    # (unlike ORAI1) -- TRPC6's gene span (425 kb) is LONGER than STIM1's
    # (238 kb), arguing against the short-gene/intronic-read dropout
    # mechanism, and the raw fit's implied ON-fraction at threshold=0.5
    # (~13.6%) is already close to the measured snRNA-seq pct_expressing
    # (9.9%) -- internally consistent, no correction applied. Gates the
    # dynamic Trpc6 channel node (ModeloPASMC_v8.bnd v8.1 edit): a cell
    # without enough TRPC6 cannot open the channel regardless of upstream
    # Plc/Nfat/Klf4 signaling (see escape_route_knockout_with_et1.py /
    # amlodipine_adaptation_mechanism.md for why this matters -- TRPC6 was
    # found to be the dominant chronic escape route under chronic ET-1 tone,
    # consistent with the TRPC6 promoter SNP literature, PMID 19380626).
    "TRPC6":  dict(alpha=0.09272146957651581, beta=4.655985304504226, vmin=0.0, vmax=15.0),
}


def augment_cells_csv_with_heterogeneity(cells_csv_path: Path, out_path: Path, seed: int) -> None:
    """Add custom:<gene>_level columns to a v2 cells.csv, sampling each
    Smooth_Muscle_Cell / Synthetic_SMC row's STIM1/ORAI1/ATP2A2 level
    independently from the real PAH-fitted Beta distribution above --
    giving every agent its own value instead of one population-level
    average. Fibroblast/Endothelial_Cell rows get blank entries (skipped by
    PhysiCell's CSV loader, which leaves their custom_data at its default).
    """
    import csv
    rng = np.random.default_rng(seed)
    with open(cells_csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    assert header[:4] == ["x", "y", "z", "type"], f"unexpected cells.csv header: {header}"
    gene_cols = [f"custom:{g.lower()}_level" for g in REAL_GENE_BETA_PARAMS_PAH] + ["custom:padding_unused"]
    new_header = header + gene_cols
    smc_types = {"Smooth_Muscle_Cell", "Synthetic_SMC"}
    out_rows = []
    for row in rows:
        row = list(row)
        if row[3] in smc_types:
            for gene, p in REAL_GENE_BETA_PARAMS_PAH.items():
                draw = rng.beta(p["alpha"], p["beta"])
                raw = p["vmin"] + draw * (p["vmax"] - p["vmin"])
                row.append(f"{raw:.6f}")
            row.append("0.0")  # padding_unused -- see heterogeneity_inputs_xml comment
        else:
            row.extend([""] * len(gene_cols))
        out_rows.append(row)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(new_header)
        writer.writerows(out_rows)

def fill_disk_with_rings(r_outer: float, cell_diameter: float, packing_fraction: float,
                          rng: np.random.Generator, jitter_frac: float = 0.15,
                          theta_min: float = 0.0, theta_max: float = 2 * np.pi,
                          theta_phase: float = 0.0) -> np.ndarray:
    """Fill the angular wedge [theta_min, theta_max) of the disk [0, r_outer) with concentric rings
    of cells -- same brick-offset/jitter scheme as build_vessel_wedge.py's fill_annulus_with_rings,
    specialized to r_inner=0 (a solid disk/wedge, not an annulus) for the Lumen_Core rigid plug.
    Defaults (full circle) are for --geometry=annulus; pass the SAME theta_min/theta_max used for
    the surrounding tissue's --geometry=wedge, or Lumen_Core will fill angles with no EC/SMC/
    Fibroblast around them. theta_phase is a whole-layer angular shift in units of one ring's
    angular cell spacing (0.5 = half-cell brick offset) used when replicating Lumen_Core across
    z layers. Returns an (N, 2) array of (x, y)."""
    wedge_span = theta_max - theta_min
    n_rings = max(1, round(r_outer / cell_diameter))
    points = []
    for i in range(n_rings):
        r = (i + 0.5) * r_outer / n_rings
        if r < 1e-6:
            if wedge_span >= 2 * np.pi - 1e-9:
                points.append((0.0, 0.0))
            continue
        arc_length = wedge_span * r
        spacing = cell_diameter / max(packing_fraction, 1e-6)
        n_in_ring = max(1, int(np.floor(arc_length / spacing)))
        angle_offset = (wedge_span / n_in_ring / 2.0) if (i % 2 == 1) else 0.0
        thetas = np.linspace(theta_min, theta_max, n_in_ring, endpoint=False) + angle_offset
        thetas = thetas + theta_phase * (wedge_span / n_in_ring)
        thetas = thetas + rng.uniform(-jitter_frac, jitter_frac, size=n_in_ring) * (wedge_span / n_in_ring)
        radii = r + rng.uniform(-jitter_frac, jitter_frac, size=n_in_ring) * (r_outer / n_rings)
        points.extend(zip(radii * np.cos(thetas), radii * np.sin(thetas)))
    return np.array(points) if points else np.empty((0, 2))


def detect_z_layers(cells_csv_path: Path):
    """Return the distinct z values present in a v2 cells.csv, as an ordered list of
    (float_value, original_token) pairs. A single-layer (2D) file yields one entry --
    the trigger for emitting the original use_2D=true one-voxel z slab; a z-extruded
    file (build_vessel_wedge.py --z-layers N) yields N exact layer planes, used to
    size the 3D z domain and replicate Lumen_Core axially."""
    import csv
    with open(cells_csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header[:4] == ["x", "y", "z", "type"], f"unexpected cells.csv header: {header}"
        seen = {}
        for row in reader:
            if row and row[2] != "":
                seen.setdefault(round(float(row[2]), 4), row[2])
    return sorted(seen.items())


def append_lumen_core_rows(cells_csv_path: Path, lumen_radius: float, diameter: float,
                            packing_fraction: float, seed: int,
                            theta_min: float = 0.0, theta_max: float = 2 * np.pi) -> None:
    """Append Lumen_Core rows filling the [theta_min, theta_max) wedge of r < lumen_radius to an
    already-written cells.csv (run AFTER augment_cells_csv_with_heterogeneity / the plain copy, so
    it works whether or not per-agent heterogeneity columns are present -- extra columns, if any,
    are padded blank for these rows, matching the blank-padding augment_cells_csv_with_heterogeneity
    already does for Fibroblast/EC).

    The disk is replicated at EVERY distinct z plane already present in the tissue rows -- a 2D
    (single-plane) cells.csv gets the original single Lumen_Core disk (byte-identical z token from
    the tissue rows themselves); a z-extruded cells.csv gets one disk per axial layer, brick-offset
    (half-cell angular phase) on alternate layers, forming a 3D cylindrical core plug."""
    import csv
    with open(cells_csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    assert header[:4] == ["x", "y", "z", "type"], f"unexpected cells.csv header: {header}"
    z_planes = detect_z_layers(cells_csv_path)
    rng = np.random.default_rng(seed + 1)  # offset from the heterogeneity seed, not reused
    n_extra_cols = len(header) - 4
    for k, (_z_val, z_token) in enumerate(z_planes):
        theta_phase = 0.5 if (k % 2 == 1) else 0.0
        xy = fill_disk_with_rings(lumen_radius, diameter, packing_fraction, rng,
                                  theta_min=theta_min, theta_max=theta_max, theta_phase=theta_phase)
        for x, y in xy:
            rows.append([f"{x:.4f}", f"{y:.4f}", z_token, "Lumen_Core"] + [""] * n_extra_cols)
    with open(cells_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


# Minimal phenotype block (no PhysiBoSS) for Fibroblast / Endothelial_Cell,
# with per-cell-definition secretion/uptake overridden by the caller.
def bystander_phenotype_block(secretion_xml: str, total_volume: float = 1000.0,
                               adhesion_strength: float = 0.4, repulsion_strength: float = 10.0) -> str:
    return f"""\
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0.0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">0</death_rate>
                        <phase_durations units="min"><duration index="0" fixed_duration="true">516</duration></phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min"><duration index="0" fixed_duration="true">0</duration><duration index="1" fixed_duration="true">86400</duration></phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <!-- See smc_phenotype_block's <volume> comment -- same fix (regrowth rates were
                         0.0, a confirmed bug: cells never regrow toward target volume after division,
                         shrinking by half every generation with no floor). Fibroblast/EC/Lumen_Core
                         divide rarely-to-never, so this mattered less for them, but fixing it here too
                         for consistency with the standard PhysiCell regrowth rates. -->
                    <total units="micron^3">{total_volume:.2f}</total>
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">1.666667e-02</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">5.833333e-03</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">{adhesion_strength}</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">{repulsion_strength}</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <cell_BM_adhesion_strength units="micron/min">4.0</cell_BM_adhesion_strength>
                    <cell_BM_repulsion_strength units="micron/min">10.0</cell_BM_repulsion_strength>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.0</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">0.0</migration_bias>
                    <options><enabled>false</enabled><use_2D>true</use_2D>
                        <chemotaxis><enabled>false</enabled><substrate>O2</substrate><direction>1</direction></chemotaxis>
                        <advanced_chemotaxis><enabled>false</enabled><normalize_each_gradient>false</normalize_each_gradient>
                            <chemotactic_sensitivities><chemotactic_sensitivity substrate="O2">0.0</chemotactic_sensitivity></chemotactic_sensitivities>
                        </advanced_chemotaxis>
                    </options>
                </motility>
                <secretion>
{secretion_xml}
                </secretion>
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0.0</dead_phagocytosis_rate>
                    <damage_rate units="1/min">1.0</damage_rate>
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates></transformation_rates>
                </cell_transformations>
"""


def secretion_substrate(name, rate, target, uptake):
    return (f'                    <substrate name="{name}">\n'
            f'                        <secretion_rate units="1/min">{rate}</secretion_rate>\n'
            f'                        <secretion_target units="substrate density">{target}</secretion_target>\n'
            f'                        <uptake_rate units="1/min">{uptake}</uptake_rate>\n'
            f'                        <net_export_rate units="total substrate/min">0</net_export_rate>\n'
            f'                    </substrate>\n')


def smc_phenotype_block(amlo_threshold, et1_threshold, o2_threshold, o2_inact_threshold,
                         division_value, base_division_value, feedback_on,
                         et1_secretion_rate, o2_uptake_rate,
                         is_founder=False, founder_migration_speed=0.0, migration_bias=0.0,
                         pathway_node="CREC_1", apoptosis_rate=0.0,
                         transformation_target=None, transformation_rate_value=0.0,
                         transformation_node="Klf4", heterogeneity_on=False,
                         o2_baseline_uptake_rate=10.0,
                         smc_adhesion_strength=0.4, smc_repulsion_strength=10.0):
    secretion_xml = (
        secretion_substrate("AMLO", 0, 1, 0) +
        secretion_substrate("ET1", 0, 1, 0) +
        secretion_substrate("O2", 0, 38, o2_baseline_uptake_rate) +  # SMCs consume O2 at a baseline rate (tunable)
        secretion_substrate("K_ALTO", 0, 1, 0)
    )
    feedback_outputs = ""
    if feedback_on:
        feedback_outputs = f"""
                        <output physicell_name="ET1 secretion" intracellular_name="Nfat">
                            <settings>
                                <action>activation</action>
                                <value>{et1_secretion_rate:.6f}</value>
                                <base_value>0.0</base_value>
                                <smoothing>0</smoothing>
                            </settings>
                        </output>
                        <output physicell_name="O2 uptake" intracellular_name="Ca_global">
                            <settings>
                                <action>activation</action>
                                <value>{o2_uptake_rate:.6f}</value>
                                <base_value>10.0</base_value>
                                <smoothing>0</smoothing>
                            </settings>
                        </output>"""
    founder_migration_output = ""
    if is_founder:
        # CREC_1 -> "migration speed": only the Synthetic_SMC founder
        # subpopulation gains directed motility when its Boolean network is
        # in the CREC_1=1 state -- paralleling the Notch3+ founder-cell
        # literature's "neointimal founder cells" that leave the media and
        # invade the (here, empty) lumen. Direction is set by a custom
        # update_migration_bias function (toward the lumen center), since
        # PhysiCell's native XML mapping has no "direction toward a point"
        # primitive.
        founder_migration_output = f"""
                        <output physicell_name="migration speed" intracellular_name="{pathway_node}">
                            <settings>
                                <action>activation</action>
                                <value>{founder_migration_speed:.6f}</value>
                                <base_value>0.0</base_value>
                                <smoothing>0</smoothing>
                            </settings>
                        </output>"""
    # v8 -- Boolean-driven phenotype transformation (Myocd/Klf4 toggle),
    # replacing the old fixed-fraction "clone" seeding: cells transform
    # between Smooth_Muscle_Cell and Synthetic_SMC at runtime based on their
    # OWN current Boolean state (Klf4 drives SMC->Synthetic_SMC; Myocd drives
    # the reverse), not a population split decided at t=0.
    transformation_rates_xml = ""
    transformation_output_xml = ""
    if transformation_target is not None:
        transformation_rates_xml = (
            f'<transformation_rate name="{transformation_target}" units="1/min">0.0</transformation_rate>'
        )
        transformation_output_xml = f"""
                        <output physicell_name="transform to {transformation_target}" intracellular_name="{transformation_node}">
                            <settings>
                                <action>activation</action>
                                <value>{transformation_rate_value:.8g}</value>
                                <base_value>0.0</base_value>
                                <smoothing>0</smoothing>
                            </settings>
                        </output>"""
    # v8 -- Resist_apop (Nfat OR Klf4 driven) now sets the agent's OWN
    # apoptosis rate continuously, superseding the old --founder-apoptosis-rate
    # CLI constant: apoptosis_rate (passed in) becomes the BASE_VALUE (what a
    # cell with Resist_apop=0 experiences), and Resist_apop=1 drives it toward
    # ~0 (apoptosis-resistant), for BOTH Smooth_Muscle_Cell and Synthetic_SMC.
    resist_apop_output_xml = f"""
                        <output physicell_name="apoptosis" intracellular_name="Resist_apop">
                            <settings>
                                <action>activation</action>
                                <value>0.0</value>
                                <base_value>{apoptosis_rate:.8g}</base_value>
                                <smoothing>0</smoothing>
                            </settings>
                        </output>"""
    # Per-agent continuous heterogeneity (real snRNA-seq-fitted STIM1/ORAI1/
    # ATP2A2 levels, section "SMC agent representation"): each agent's own
    # custom:<gene>_level (set per-row in cells.csv, see
    # augment_cells_csv_with_heterogeneity) drives that SAME agent's free
    # input node via the native signal/behavior dictionary -- the identical
    # mechanism already used for AMLO/ET1/HIPOXIA/K_ALTO, just sourced from
    # custom_data instead of a microenvironment substrate. No new C++.
    heterogeneity_inputs_xml = ""
    heterogeneity_custom_data_xml = ""
    if heterogeneity_on:
        heterogeneity_inputs_xml = "".join(f"""
                        <input physicell_name="custom:{gene.lower()}_level" intracellular_name="{gene}">
                            <settings><action>activation</action><threshold>0.5</threshold><smoothing>0</smoothing></settings>
                        </input>""" for gene in REAL_GENE_BETA_PARAMS_PAH)
        heterogeneity_custom_data_xml = "".join(
            f'\n                <{gene.lower()}_level conserved="false" units="dimensionless" description="">0.0</{gene.lower()}_level>'
            for gene in REAL_GENE_BETA_PARAMS_PAH
        )
        # WORKAROUND for a confirmed PhysiBoSS bug: the LAST custom_data-sourced
        # <input> declared in a cell_definition's intracellular mapping never
        # receives its per-agent signal value correctly (reads as if always
        # below threshold, regardless of the real custom_data value) --
        # reproduced with every one of the 3 real genes in that slot, and with
        # a renamed dummy variable in that slot, ruling out a name-specific
        # cause. Confirmed fixed by appending a disposable 4th custom-data
        # input that absorbs the bug: mapped to CPA (an existing, otherwise
        # always-off scenario node in this model) with a threshold far above
        # any value this dummy signal will ever take, so it never actually
        # changes CPA's behavior even though ITS OWN read is unreliable.
        heterogeneity_inputs_xml += """
                        <input physicell_name="custom:padding_unused" intracellular_name="CPA">
                            <settings><action>activation</action><threshold>999999</threshold><smoothing>0</smoothing></settings>
                        </input>"""
        heterogeneity_custom_data_xml += (
            '\n                <padding_unused conserved="false" units="dimensionless" description="">0.0</padding_unused>'
        )
    return f"""\
                <cycle code="5" name="live">
                    <phase_transition_rates units="1/min">
                        <rate start_index="0" end_index="0" fixed_duration="false">0.0</rate>
                    </phase_transition_rates>
                </cycle>
                <death>
                    <model code="100" name="apoptosis">
                        <death_rate units="1/min">{apoptosis_rate:.8g}</death_rate>
                        <phase_durations units="min"><duration index="0" fixed_duration="true">516</duration></phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                    <model code="101" name="necrosis">
                        <death_rate units="1/min">0.0</death_rate>
                        <phase_durations units="min"><duration index="0" fixed_duration="true">0</duration><duration index="1" fixed_duration="true">86400</duration></phase_durations>
                        <parameters>
                            <unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
                            <lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
                            <cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
                            <nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
                            <calcification_rate units="1/min">0</calcification_rate>
                            <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                        </parameters>
                    </model>
                </death>
                <volume>
                    <!-- Regrowth rates (fluid/cytoplasmic/nuclear_biomass_change_rate) were 0.0 here --
                         a confirmed bug: standard_volume_update_function's three relaxation terms are
                         each dt*RATE*(target-current), so RATE=0 means volume NEVER changes except at
                         division (which halves it). Over repeated generations this geometrically
                         shrinks cells with no floor -- confirmed in a 45-day run: median SMC volume
                         fell to 220.9 um^3 (target/8, 3rd generation) and the minimum to 6.9 um^3
                         (target/256, 8th generation) by day 45. Since mechanics radius/repulsion derive
                         from volume, this silently undersized cells and UNDERSTATED how much the
                         repulsion/dt fixes earlier this session needed to do -- restored to PhysiCell's
                         standard regrowth rates (same values already used, correctly, for the
                         apoptosis/necrosis SHRINKAGE-while-dying rates just above) so daughters regrow
                         toward the target volume between divisions instead of halving forever. -->
                    <total units="micron^3">1767.15</total>  <!-- 15 um diameter, matches build_vessel_cross_section.py --smc-diameter -->
                    <fluid_fraction units="dimensionless">0.75</fluid_fraction>
                    <nuclear units="micron^3">540</nuclear>
                    <fluid_change_rate units="1/min">0.05</fluid_change_rate>
                    <cytoplasmic_biomass_change_rate units="1/min">1.666667e-02</cytoplasmic_biomass_change_rate>
                    <nuclear_biomass_change_rate units="1/min">5.833333e-03</nuclear_biomass_change_rate>
                    <calcified_fraction units="dimensionless">0</calcified_fraction>
                    <calcification_rate units="1/min">0</calcification_rate>
                    <relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
                </volume>
                <mechanics>
                    <cell_cell_adhesion_strength units="micron/min">{smc_adhesion_strength}</cell_cell_adhesion_strength>
                    <cell_cell_repulsion_strength units="micron/min">{smc_repulsion_strength}</cell_cell_repulsion_strength>
                    <relative_maximum_adhesion_distance units="dimensionless">1.25</relative_maximum_adhesion_distance>
                    <options>
                        <set_relative_equilibrium_distance enabled="false" units="dimensionless">1.8</set_relative_equilibrium_distance>
                        <set_absolute_equilibrium_distance enabled="false" units="micron">15.12</set_absolute_equilibrium_distance>
                    </options>
                    <cell_BM_adhesion_strength units="micron/min">4.0</cell_BM_adhesion_strength>
                    <cell_BM_repulsion_strength units="micron/min">10.0</cell_BM_repulsion_strength>
                    <attachment_elastic_constant units="1/min">0.01</attachment_elastic_constant>
                    <attachment_rate units="1/min">0.0</attachment_rate>
                    <detachment_rate units="1/min">0.0</detachment_rate>
                </mechanics>
                <motility>
                    <speed units="micron/min">0.0</speed>
                    <persistence_time units="min">5</persistence_time>
                    <migration_bias units="dimensionless">{migration_bias if is_founder else 0.0}</migration_bias>
                    <options><enabled>{str(is_founder).lower()}</enabled><use_2D>true</use_2D>
                        <chemotaxis><enabled>{str(is_founder).lower()}</enabled><substrate>ET1</substrate><direction>1</direction></chemotaxis>
                        <advanced_chemotaxis><enabled>false</enabled><normalize_each_gradient>false</normalize_each_gradient>
                            <chemotactic_sensitivities><chemotactic_sensitivity substrate="O2">0.0</chemotactic_sensitivity></chemotactic_sensitivities>
                        </advanced_chemotaxis>
                    </options>
                </motility>
                <secretion>
{secretion_xml}
                </secretion>
                <cell_interactions>
                    <dead_phagocytosis_rate units="1/min">0.0</dead_phagocytosis_rate>
                    <damage_rate units="1/min">1.0</damage_rate>
                </cell_interactions>
                <cell_transformations>
                    <transformation_rates>{transformation_rates_xml}</transformation_rates>
                </cell_transformations>
                <intracellular type="maboss">
                    <bnd_filename>./config/ModeloPASMC_v8.bnd</bnd_filename>
                    <cfg_filename>./config/ModeloPASMC_v8.cfg</cfg_filename>
                    <settings>
                        <intracellular_dt>60</intracellular_dt>
                        <time_stochasticity>0.5</time_stochasticity>
                        <scaling>60</scaling>
                        <start_time>0.0</start_time>
                        <inheritance global="True" />
                    </settings>
                    <mapping>{heterogeneity_inputs_xml}
                        <input physicell_name="AMLO" intracellular_name="AMLO">
                            <settings><action>activation</action><threshold>{amlo_threshold}</threshold><smoothing>0</smoothing></settings>
                        </input>
                        <input physicell_name="ET1" intracellular_name="ET1">
                            <settings><action>activation</action><threshold>{et1_threshold}</threshold><smoothing>0</smoothing></settings>
                        </input>
                        <input physicell_name="O2" intracellular_name="HIPOXIA">
                            <settings><action>inhibition</action><threshold>{o2_threshold}</threshold><inact_threshold>{o2_inact_threshold}</inact_threshold><smoothing>0</smoothing></settings>
                        </input>
                        <input physicell_name="K_ALTO" intracellular_name="K_ALTO">
                            <settings><action>activation</action><threshold>0.5</threshold><smoothing>0</smoothing></settings>
                        </input>
                        <output physicell_name="cycle entry" intracellular_name="{pathway_node}">
                            <settings><action>activation</action><value>{division_value:.8f}</value><base_value>{base_division_value:.8f}</base_value><smoothing>0</smoothing></settings>
                        </output>{feedback_outputs}{founder_migration_output}{transformation_output_xml}{resist_apop_output_xml}
                    </mapping>
                </intracellular>
            </phenotype>
            <custom_data>{heterogeneity_custom_data_xml}
                <pressure_debug conserved="false" units="dimensionless" description="diagnostic: pCell-&gt;state.simple_pressure, written every post_update_intracellular call">0.0</pressure_debug>
                <placeholder conserved="false" units="dimensionless" description="">1.0</placeholder>
            </custom_data>
        </cell_definition>
"""


def build_settings_xml(args) -> str:
    hw = args.domain_half_width
    feedback_on = (args.feedback == "on")

    smc_block = ('        <cell_definition name="Smooth_Muscle_Cell" ID="0">\n            <phenotype>\n'
                 + smc_phenotype_block(args.amlo_threshold, args.et1_threshold, args.o2_threshold, args.o2_inact_threshold,
                                        args.division_rate, args.base_division_rate, feedback_on,
                                        args.et1_secretion_rate, args.o2_uptake_rate,
                                        apoptosis_rate=args.smc_apoptosis_rate,
                                        transformation_target="Synthetic_SMC",
                                        transformation_rate_value=args.transformation_rate,
                                        transformation_node="Klf4",
                                        heterogeneity_on=args.per_agent_heterogeneity,
                                        o2_baseline_uptake_rate=args.o2_baseline_uptake_rate,
                                        smc_adhesion_strength=args.smc_adhesion_strength,
                                        smc_repulsion_strength=args.smc_repulsion_strength))
    synthetic_block = ('        <cell_definition name="Synthetic_SMC" ID="3">\n            <phenotype>\n'
                     + smc_phenotype_block(args.amlo_threshold, args.et1_threshold, args.o2_threshold, args.o2_inact_threshold,
                                            args.founder_division_rate, args.base_division_rate, feedback_on,
                                            args.et1_secretion_rate, args.o2_uptake_rate,
                                            is_founder=True, founder_migration_speed=args.founder_migration_speed,
                                            migration_bias=args.migration_bias,
                                            pathway_node=args.founder_pathway_node,
                                            apoptosis_rate=args.smc_apoptosis_rate,
                                            transformation_target="Smooth_Muscle_Cell",
                                            transformation_rate_value=args.reverse_transformation_rate,
                                            transformation_node="Myocd",
                                            heterogeneity_on=args.per_agent_heterogeneity,
                                            o2_baseline_uptake_rate=args.o2_baseline_uptake_rate,
                                        smc_adhesion_strength=args.smc_adhesion_strength,
                                        smc_repulsion_strength=args.smc_repulsion_strength))

    # Vasa-vasorum proxy: real medium/large arteries (especially a remodeled, thickened PAH wall)
    # are NOT oxygenated by luminal diffusion alone beyond the first ~100-150um (a typical tissue
    # O2 diffusion/Krogh length) -- they carry a dedicated adventitial microvascular supply (the vasa
    # vasorum) specifically because the wall is too thick for that. This model's 272um
    # (EC+media+adventitia) PAH wall is 2-3x a representative diffusion length, so without some
    # secondary source, nearly the entire wall is predicted hypoxic by pure single-source diffusion
    # physics -- confirmed directly (Hif1a ON in 99.6% of SMC-lineage cells at the old O2 params).
    # Representing vasa vasorum as a small secretion_rate on Fibroblast (which populate exactly the
    # adventitial layer where vasa vasorum actually sit) is more robust than a single domain-edge
    # Dirichlet value: it is a source DISTRIBUTED throughout the outer wall (scales with the
    # adventitia's fixed ~367-cell count, which --unlike the SMC population -- never grows/shrinks
    # with transformation or death), rather than a single far-field boundary value whose effect on
    # the hypoxic fraction drifted substantially as SMC density changed during tuning.
    fibro_secretion = secretion_substrate("AMLO", 0, 1, 0) + secretion_substrate("ET1", 0, 1, 0) + \
        secretion_substrate("O2", args.o2_vasa_vasorum_secretion_rate, args.o2_vasa_vasorum_target,
                             args.o2_fibro_uptake_rate) + secretion_substrate("K_ALTO", 0, 1, 0)
    fibro_block = ('        <cell_definition name="Fibroblast" ID="1">\n            <phenotype>\n'
                    + bystander_phenotype_block(fibro_secretion, total_volume=1767.15)  # 15 um diameter, matches build_vessel_cross_section.py's --fibroblast-diameter
                    + '            </phenotype>\n            <custom_data>\n'
                      '                <placeholder conserved="false" units="dimensionless" description="">1.0</placeholder>\n'
                      '            </custom_data>\n        </cell_definition>\n')

    ec_secretion = (secretion_substrate("AMLO", args.amlo_secretion_rate, 1.0, 0) +
                     secretion_substrate("ET1", args.et1_secretion_rate_ec, 1.0, 0) +
                     secretion_substrate("O2", args.o2_secretion_rate, 38, 0) +
                     secretion_substrate("K_ALTO", 0, 1, 0))
    ec_block = ('        <cell_definition name="Endothelial_Cell" ID="2">\n            <phenotype>\n'
                + bystander_phenotype_block(ec_secretion, total_volume=904.78,  # 12 um diameter, matches build_vessel_cross_section.py's --ec-diameter
                                            adhesion_strength=args.ec_barrier_adhesion, repulsion_strength=args.ec_barrier_repulsion)
                # Endothelium as a tight, non-traversable barrier: neointima
                # formation happens in the SUBENDOTHELIAL space (between the
                # endothelium and the internal elastic lamina), not by SMCs
                # physically crossing past an intact endothelial monolayer
                # into the open lumen -- strong EC-EC adhesion/repulsion
                # keeps the ring cohesive and mechanically resistant to being
                # pushed through by invading founder cells.
                + '            </phenotype>\n            <custom_data>\n'
                  '                <placeholder conserved="false" units="dimensionless" description="">1.0</placeholder>\n'
                  '            </custom_data>\n        </cell_definition>\n')

    lumen_core_secretion = (secretion_substrate("AMLO", 0, 1, 0) + secretion_substrate("ET1", 0, 1, 0) +
                            secretion_substrate("O2", 0, 38, 0) + secretion_substrate("K_ALTO", 0, 1, 0))
    lumen_core_block = ('        <cell_definition name="Lumen_Core" ID="4">\n            <phenotype>\n'
                         + bystander_phenotype_block(lumen_core_secretion, total_volume=904.78,
                                                      adhesion_strength=args.lumen_core_adhesion_strength,
                                                      repulsion_strength=args.lumen_core_repulsion_strength)
                         # Chemically inert (zero secretion/uptake on every substrate) so it doesn't
                         # perturb the AMLO/ET1/O2/K_ALTO fields -- purely a mechanical rigid plug.
                         # is_movable=false is set in custom.cpp (not exposed in this XML schema).
                         + '            </phenotype>\n            <custom_data>\n'
                           '                <placeholder conserved="false" units="dimensionless" description="">1.0</placeholder>\n'
                           '            </custom_data>\n        </cell_definition>\n')

    k_alto_init = 1.0 if args.k_alto_on else 0.0

    # z domain: single-plane (2D) cells.csv -> the original one-voxel z slab, use_2D=true
    # (byte-identical to all prior runs). A z-extruded cells.csv (build_vessel_wedge.py
    # --z-layers N, detected here from the CSV's distinct exact z values) -> use_2D=false
    # with a z range covering all layer planes plus a half-layer and a 10 um margin on
    # each side; virtual_wall_at_domain_edge then makes zmin/zmax reflective walls -- the
    # correct zero-flux symmetry-plane treatment for an axial segment of an axisymmetric
    # vessel (mirroring the adjacent axial segment reproduces the same boundary forces/
    # fluxes a longer segment would produce), exactly analogous to the wedge's two angular
    # cut edges. Substrate Dirichlet stays OFF on zmin/zmax in either mode for the same
    # reason (see the O2 block's wedge comment).
    z_planes = getattr(args, "_z_planes", None) or [(0.0, "0.0")]
    if len(z_planes) > 1:
        z_vals = [v for v, _t in z_planes]
        gaps = [b - a for a, b in zip(z_vals, z_vals[1:])]
        z_spacing = float(np.median(gaps)) if gaps else 15.0
        z_lo = min(z_vals) - z_spacing / 2.0 - 10.0
        z_hi = max(z_vals) + z_spacing / 2.0 + 10.0
        z_block = (f"<z_min>{z_lo:.1f}</z_min><z_max>{z_hi:.1f}</z_max>\n"
                   f"        <dx>20</dx><dy>20</dy><dz>20</dz>\n"
                   f"        <!-- 3D z-extrusion: {len(z_planes)} axial layer planes detected in cells.csv "
                   f"(spacing ~{z_spacing:.1f} um) -->\n"
                   f"        <use_2D>false</use_2D>")
    else:
        z_block = ("<z_min>-10</z_min><z_max>10</z_max>\n"
                   "        <dx>20</dx><dy>20</dy><dz>20</dz>\n"
                   "        <use_2D>true</use_2D>")

    return f"""<PhysiCell_settings version="devel-version">

    <domain>
        <x_min>{0 if args.geometry == "wedge" else -hw}</x_min><x_max>{hw}</x_max>
        <y_min>{0 if args.geometry == "wedge" else -hw}</y_min><y_max>{hw}</y_max>
        {z_block}
    </domain>

    <overall>
        <max_time units="min">{args.max_time}</max_time>
        <time_units>min</time_units>
        <space_units>micron</space_units>
        <dt_diffusion units="min">0.05</dt_diffusion>
        <dt_mechanics units="min">{args.dt_mechanics}</dt_mechanics>
        <dt_phenotype units="min">1</dt_phenotype>
    </overall>

    <parallel><omp_num_threads>{args.n_threads}</omp_num_threads></parallel>

    <save>
        <folder>output</folder>
        <full_data><interval units="min">{args.save_interval}</interval><enable>true</enable></full_data>
        <SVG><interval units="min">{args.save_interval}</interval><enable>{str(args.save_svg).lower()}</enable></SVG>
        <legacy_data><enable>true</enable></legacy_data>
    </save>

    <options>
        <legacy_random_points_on_sphere_in_divide>false</legacy_random_points_on_sphere_in_divide>
        <virtual_wall_at_domain_edge>true</virtual_wall_at_domain_edge>
        <disable_automated_spring_adhesions>true</disable_automated_spring_adhesions>
    </options>

    <microenvironment_setup>
        <variable name="AMLO" units="dimensionless" ID="0">
            <physical_parameter_set>
                <diffusion_coefficient units="micron^2/min">{args.amlo_diffusion_coefficient}</diffusion_coefficient>
                <decay_rate units="1/min">{args.amlo_decay_rate}</decay_rate>
            </physical_parameter_set>
            <initial_condition units="dimensionless">0.0</initial_condition>
            <Dirichlet_boundary_condition units="dimensionless" enabled="false">0.0</Dirichlet_boundary_condition>
            <Dirichlet_options>
                <boundary_value ID="xmin" enabled="false">0</boundary_value><boundary_value ID="xmax" enabled="false">0</boundary_value>
                <boundary_value ID="ymin" enabled="false">0</boundary_value><boundary_value ID="ymax" enabled="false">0</boundary_value>
                <boundary_value ID="zmin" enabled="false">0</boundary_value><boundary_value ID="zmax" enabled="false">0</boundary_value>
            </Dirichlet_options>
        </variable>
        <variable name="ET1" units="dimensionless" ID="1">
            <physical_parameter_set>
                <diffusion_coefficient units="micron^2/min">1500</diffusion_coefficient>
                <decay_rate units="1/min">0.01</decay_rate>
            </physical_parameter_set>
            <initial_condition units="dimensionless">0.0</initial_condition>
            <Dirichlet_boundary_condition units="dimensionless" enabled="false">0.0</Dirichlet_boundary_condition>
            <Dirichlet_options>
                <boundary_value ID="xmin" enabled="false">0</boundary_value><boundary_value ID="xmax" enabled="false">0</boundary_value>
                <boundary_value ID="ymin" enabled="false">0</boundary_value><boundary_value ID="ymax" enabled="false">0</boundary_value>
                <boundary_value ID="zmin" enabled="false">0</boundary_value><boundary_value ID="zmax" enabled="false">0</boundary_value>
            </Dirichlet_options>
        </variable>
        <variable name="O2" units="mmHg" ID="2">
            <physical_parameter_set>
                <diffusion_coefficient units="micron^2/min">{args.o2_diffusion_coefficient}</diffusion_coefficient>
                <decay_rate units="1/min">{args.o2_decay_rate}</decay_rate>
            </physical_parameter_set>
            <initial_condition units="mmHg">38.0</initial_condition>
            <Dirichlet_boundary_condition units="mmHg" enabled="{str(args.o2_dirichlet_wall).lower()}">{args.o2_dirichlet_value}</Dirichlet_boundary_condition>
            <Dirichlet_options>
                <!-- Wedge geometry: xmin/ymin are the angular-cut (symmetry) edges, not the true
                     outer boundary -- Dirichlet must stay OFF there (zero-flux is the physically
                     correct symmetry condition; see fill_annulus_with_rings docstring in
                     build_vessel_wedge.py). Only xmax/ymax (the true far/adventitial boundary) get
                     the Dirichlet value. For the full annulus (geometry=="annulus"), all four
                     faces are the same "outer" boundary and all follow --o2-dirichlet-wall as before. -->
                <boundary_value ID="xmin" enabled="{"false" if args.geometry == "wedge" else str(args.o2_dirichlet_wall).lower()}">0</boundary_value><boundary_value ID="xmax" enabled="{str(args.o2_dirichlet_wall).lower()}">{args.o2_dirichlet_value}</boundary_value>
                <boundary_value ID="ymin" enabled="{"false" if args.geometry == "wedge" else str(args.o2_dirichlet_wall).lower()}">0</boundary_value><boundary_value ID="ymax" enabled="{str(args.o2_dirichlet_wall).lower()}">{args.o2_dirichlet_value}</boundary_value>
                <boundary_value ID="zmin" enabled="false">0</boundary_value><boundary_value ID="zmax" enabled="false">0</boundary_value>
            </Dirichlet_options>
        </variable>
        <variable name="K_ALTO" units="dimensionless" ID="3">
            <physical_parameter_set>
                <diffusion_coefficient units="micron^2/min">4000</diffusion_coefficient>
                <decay_rate units="1/min">0.0</decay_rate>
            </physical_parameter_set>
            <initial_condition units="dimensionless">{k_alto_init}</initial_condition>
            <Dirichlet_boundary_condition units="dimensionless" enabled="false">0.0</Dirichlet_boundary_condition>
            <Dirichlet_options>
                <boundary_value ID="xmin" enabled="false">0</boundary_value><boundary_value ID="xmax" enabled="false">0</boundary_value>
                <boundary_value ID="ymin" enabled="false">0</boundary_value><boundary_value ID="ymax" enabled="false">0</boundary_value>
                <boundary_value ID="zmin" enabled="false">0</boundary_value><boundary_value ID="zmax" enabled="false">0</boundary_value>
            </Dirichlet_options>
        </variable>
        <options>
            <calculate_gradients>false</calculate_gradients>
            <track_internalized_substrates_in_each_agent>false</track_internalized_substrates_in_each_agent>
            <initial_condition enabled="false" type="matlab"><filename>./config/initial.mat</filename></initial_condition>
            <dirichlet_nodes enabled="false" type="matlab"><filename>./config/dirichlet.mat</filename></dirichlet_nodes>
        </options>
    </microenvironment_setup>

    <cell_definitions>
{smc_block}
{fibro_block}
{ec_block}
{synthetic_block}
{lumen_core_block if args.fill_lumen_core else ""}
    </cell_definitions>

    <initial_conditions>
        <cell_positions enabled="true" type="csv">
            <folder>./config</folder>
            <filename>cells.csv</filename>
        </cell_positions>
    </initial_conditions>

    <user_parameters>
        <random_seed type="int" units="dimensionless" description="">{args.seed}</random_seed>
        <node_to_visualize type="string" units="" description="">CREC_1</node_to_visualize>
        <pressure_half_max type="double" units="dimensionless" description="PhysiCell simple_pressure value at which cycle-entry rate is suppressed 2x (pressure-gated contact inhibition)">{args.pressure_half_max}</pressure_half_max>
        <pressure_hill_power type="double" units="dimensionless" description="Hill coefficient for the pressure-gating suppression of cycle entry">{args.pressure_hill_power}</pressure_hill_power>
        <population_half_max type="double" units="dimensionless" description="total agent count at which cycle-entry rate is halved (global backstop cap)">{args.population_half_max}</population_half_max>
        <population_hill_power type="double" units="dimensionless" description="Hill coefficient for the global population cap">{args.population_hill_power}</population_hill_power>
        <lumen_radius type="double" units="micron" description="hard radial clamp for founder cells -- basement membrane/IEL">{args.lumen_radius}</lumen_radius>
        <contraction_speed type="double" units="micron/min" description="CONTR_1=1 -> inward radial velocity added to Smooth_Muscle_Cell/Synthetic_SMC (reversible vasoconstriction proxy, uncalibrated placeholder)">{args.contraction_speed}</contraction_speed>
        <amlo_periodic_dosing type="bool" units="dimensionless" description="enable once-daily oral PK AMLO dosing instead of constant EC secretion_rate">{str(args.amlo_periodic_dosing).lower()}</amlo_periodic_dosing>
        <amlo_pk_ka type="double" units="1/min" description="oral absorption rate constant, solved from --amlo-tmax-hours">{args.amlo_pk_ka:.8g}</amlo_pk_ka>
        <amlo_pk_ke type="double" units="1/min" description="elimination rate constant, from --amlo-half-life-hours">{args.amlo_pk_ke:.8g}</amlo_pk_ke>
        <amlo_dose_interval_min type="double" units="min" description="redosing interval (1440 = once daily)">{args.amlo_dose_interval_min}</amlo_dose_interval_min>
        <amlo_pk_dose_amount type="double" units="dimensionless" description="dose scale = amlo_dose_mg/10 (linear PK assumption)">{args.amlo_pk_dose_amount:.8g}</amlo_pk_dose_amount>
        <amlo_pk_ss_peak type="double" units="dimensionless" description="steady-state peak of the unit-dose (10mg reference) PK curve, used to normalize EC secretion_rate">{args.amlo_pk_ss_peak:.8g}</amlo_pk_ss_peak>
        <amlo_pk_base_secretion_rate type="double" units="1/min" description="EC AMLO secretion_rate at steady-state peak for the 10mg reference dose">{args.amlo_secretion_rate}</amlo_pk_base_secretion_rate>
    </user_parameters>
</PhysiCell_settings>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--physicell-root", type=Path,
                     default=Path("/home/oth/BSC/TEBN/pah_digital_twin/physicell_runs/PhysiCell-master"))
    ap.add_argument("--cells-csv", type=Path, required=True)
    ap.add_argument("--feedback", choices=["on", "off"], default="off")
    ap.add_argument("--domain-half-width", type=float, default=420.0)
    ap.add_argument("--geometry", choices=["annulus", "wedge"], default="annulus",
                     help="'annulus' (default): full 360-degree circular cross-section, domain "
                          "centered at the origin. 'wedge': an angular-symmetry-reduction -- see "
                          "--wedge-angle-degrees and build_vessel_wedge.py -- domain occupies the "
                          "first quadrant only (x,y >= 0), with the two straight edges (x=0, y=0) "
                          "treated as reflective/zero-flux symmetry cuts rather than the true outer "
                          "boundary. Must be paired with a --cells-csv built by build_vessel_wedge.py "
                          "(NOT build_vessel_cross_section.py), or the cell layout and the domain/"
                          "Dirichlet geometry will disagree.")
    ap.add_argument("--wedge-angle-degrees", type=float, default=90.0,
                     help="Angular span of the wedge, degrees -- only used to auto-scale "
                          "--population-half-max (by wedge_angle_degrees/360) when --geometry=wedge "
                          "and --population-half-max is left at its default; must match the value "
                          "passed to build_vessel_wedge.py when generating --cells-csv.")
    ap.add_argument("--z-layers", type=int, default=None,
                     help="EXPECTED number of axial z layers in --cells-csv (build_vessel_wedge.py "
                          "--z-layers N). Purely a cross-check: the layer count is auto-detected from "
                          "the CSV's distinct z values, which also drives the 3D switch (N>1 -> "
                          "use_2D=false, z-sized domain, per-layer Lumen_Core fill, and "
                          "--population-half-max auto-scaled by N on top of the wedge fraction). "
                          "Leave unset to trust the CSV.")
    ap.add_argument("--max-time", type=float, default=4320.0)
    ap.add_argument("--save-interval", type=float, default=360.0)
    ap.add_argument("--save-svg", action="store_true")
    ap.add_argument("--n-threads", type=int, default=8,
                     help="OpenMP thread count. Raised from an earlier 4 -- this sandbox has 8 cores and was "
                          "leaving half of them idle.")
    ap.add_argument("--amlo-threshold", type=float, default=0.3)
    ap.add_argument("--et1-threshold", type=float, default=0.3)
    ap.add_argument("--o2-threshold", type=float, default=15.0, help="O2 (mmHg) below which HIPOXIA turns ON")
    ap.add_argument("--o2-inact-threshold", type=float, default=20.0, help="O2 (mmHg) above which HIPOXIA turns back OFF")
    ap.add_argument("--k-alto-on", action="store_true")
    ap.add_argument("--amlo-secretion-rate", type=float, default=0.02)
    ap.add_argument("--amlo-diffusion-coefficient", type=float, default=600.0,
                     help="AMLO diffusion coefficient, um^2/min. Default lowered from the original 4000 so the "
                          "decay length sqrt(D/decay) sits well within the vessel-wall thickness (~150-250 um) "
                          "instead of vastly exceeding the whole ~820 um domain -- giving a real radial gradient "
                          "across the SMC layer instead of near-uniform saturation (see section 10.17/10.20).")
    ap.add_argument("--amlo-decay-rate", type=float, default=0.03,
                     help="AMLO decay rate, 1/min. Default raised from the original 0.002; paired with "
                          "--amlo-diffusion-coefficient=600 this gives decay length sqrt(600/0.03)~=141 um.")
    ap.add_argument("--o2-diffusion-coefficient", type=float, default=100000.0,
                     help="O2 diffusion coefficient, um^2/min -- the standard PhysiCell tissue-oxygen value. "
                          "(Earlier this session, D was lowered to 20000 to try to force a steeper gradient via "
                          "bulk decay alone -- but the actual gradient-shaping mechanism is now the Dirichlet-0 "
                          "outer wall below, combined with a corrected, much lower per-cell uptake rate -- see "
                          "--o2-baseline-uptake-rate -- so D is reverted to the physiologically standard value.)")
    ap.add_argument("--o2-decay-rate", type=float, default=0.1,
                     help="O2 bulk decay rate, 1/min -- the standard PhysiCell tissue-oxygen value (reverted from "
                          "an earlier 0.4/0.05 tuning attempt). The real sink shaping the radial gradient is the "
                          "Dirichlet-0 outer wall below, combined with --o2-baseline-uptake-rate; bulk decay plays "
                          "only a minor supporting role at this value.")
    ap.add_argument("--o2-dirichlet-wall", type=lambda s: s.lower() not in ("0", "false", "no"), default=True,
                     help="If true (default), clamp O2 to --o2-dirichlet-value at the domain's outer boundary "
                          "(all 6 box faces) -- a Dirichlet 'wall' representing the adventitia/interstitium "
                          "boundary. Combined with the EC source at the lumen (secretion target 38 mmHg) this "
                          "creates a real source-to-sink diffusion gradient across the vessel wall, independent "
                          "of (and more robust than) tuning the bulk decay rate alone.")
    ap.add_argument("--o2-dirichlet-value", type=float, default=28.0,
                     help="O2 value (mmHg) imposed at the domain's outer boundary when --o2-dirichlet-wall is "
                          "set. History: originally 0.0 (a total-sink 'Dirichlet-0' wall), which drove the whole "
                          "wall thickness beyond ~r=100um below the 15mmHg HIPOXIA threshold by day 180 of the "
                          "long campaign (Hif1a ON in 99.6% of SMC-lineage cells; 56% of all CREC_1 activations "
                          "came via Hif1a alone, Nfat OFF -- the hypoxia/HIF arm, not the intended Ca2+/"
                          "calcineurin/NFAT/amlodipine pathway, was doing most of the work of driving "
                          "transformation). Raised to 28 mmHg (paired with --o2-baseline-uptake-rate=8, see that "
                          "argument) to treat the adventitial edge as modestly-perfused tissue instead of a "
                          "total sink -- validated (short equilibration run, initial cells.csv layout) to give a "
                          "genuine U-shaped radial profile (~30 mmHg near lumen, dipping to ~15 mmHg in the "
                          "mid-media -- a real but partial hypoxic core, consistent with diffusion-limited medial "
                          "hypoxia in a thickened PAH wall -- and recovering to ~23 mmHg at the outer wall), "
                          "putting ~50% rather than ~100% of the SMC-lineage population below the HIPOXIA "
                          "threshold. NOTE: this was tuned at the campaign's INITIAL cell density/extent; the "
                          "wall thickens further as the tissue grows toward the ~2000-2500-agent population cap, "
                          "which will likely shift the hypoxic fraction somewhat -- re-check at steady-state "
                          "density before trusting the 50% figure for the full 180-day campaign.")
    ap.add_argument("--et1-secretion-rate-ec", type=float, default=0.05)
    ap.add_argument("--o2-secretion-rate", type=float, default=300.0,
                     help="Endothelial_Cell O2 secretion rate, 1/min. Default raised sharply from an earlier 1.0: "
                          "with only 52 EC agents supplying O2 against a MUCH larger SMC/Fibroblast-consuming "
                          "population, a modest secretion_rate could not offset total tissue uptake even after "
                          "lowering --o2-baseline-uptake-rate and adding the Dirichlet-0 wall -- the real-world "
                          "equivalent is that blood flowing through the lumen is an effectively unlimited O2 "
                          "reservoir, which a small secreting-cell population only approximates at a high rate. "
                          "300 pins the EC ring near its secretion_target (38 mmHg) and lets a real radial "
                          "gradient form outward (~30-35 mmHg near the lumen -> ~4-5 mmHg at the outer wall), "
                          "instead of either near-uniform saturation (the original default) or near-total "
                          "depletion everywhere (an intermediate, now-superseded tuning attempt).")
    ap.add_argument("--division-rate", type=float, default=2e-5,
                     help="Smooth_Muscle_Cell's own CREC_1-gated division rate (1/min), ~24-day doubling. "
                          "History: originally 1e-5 (~48-day doubling); raised to 2e-4 (~2.4-day doubling) to fix "
                          "a tissue-emptying bug -- but that bug's real cause was --founder-division-rate "
                          "defaulting to a MISMATCHED, much slower value (now fixed, defaults to match this "
                          "argument), not this rate itself. 2.4-day doubling is unrealistically fast for "
                          "differentiated vascular SMC turnover: the PAH literature describes medial hyperplasia "
                          "as a MINOR contributor to wall thickening relative to cell hypertrophy/ECM deposition, "
                          "with proliferation (labeling) indices in hypoxic/MCT rat models typically single-digit "
                          "percent/day even in actively remodeling vessels -- consistent with a weeks-scale "
                          "population doubling time, not days. Lowered 10x to 2e-5 (~24-day doubling) accordingly; "
                          "re-validate population stability (no emptying) given the founder-rate bug is fixed.")
    ap.add_argument("--base-division-rate", type=float, default=1e-7)
    ap.add_argument("--pressure-half-max", type=float, default=0.1,
                     help="simple_pressure value at which cycle-entry rate is halved (pressure-gated contact "
                          "inhibition). Lowered from an earlier 0.25: once --founder-division-rate was corrected "
                          "to match --division-rate (see above) and --transformation-rate was slowed, there was "
                          "no longer an (accidental) population brake from Synthetic_SMC dividing far slower than "
                          "Smooth_Muscle_Cell -- at 0.25 the population exploded to ~24,000 agents (from ~1,200) "
                          "within 36 days, since pressure-gating alone was too permissive to compensate. 0.1 gives "
                          "gentle, bounded growth (~1,200->~1,500 over 30 days) without tipping into the opposite "
                          "failure mode (0.05/0.02 both caused net tissue shrinkage in testing).")
    ap.add_argument("--pressure-hill-power", type=float, default=4.0)
    ap.add_argument("--population-half-max", type=float, default=None,
                     help="Total agent count at which cycle-entry rate is halved -- a GLOBAL Hill-type backstop "
                          "on top of the local pressure gate above. Default: 2000 for --geometry=annulus, or "
                          "2000*(--wedge-angle-degrees/360) for --geometry=wedge (same areal density target, "
                          "scaled to the wedge's smaller cell count) -- pass an explicit value to override "
                          "either default. (History: added after local pressure-gating alone was found unable "
                          "to reliably bound growth over multi-week timescales -- a continuously expanding "
                          "low-pressure tissue edge kept outrunning local suppression, reaching 10,000-24,000 "
                          "agents within 25-65 days across several pressure_half_max/pressure_hill_power "
                          "combinations tested. 2000 sits comfortably above the validated ~1,200-1,700 "
                          "steady-state range from this project's earlier 30-45 day runs.)")
    ap.add_argument("--population-hill-power", type=float, default=8.0,
                     help="Hill coefficient for the global population cap -- sharp by default so the cap acts "
                          "as a genuine ceiling rather than a slow brake.")
    ap.add_argument("--founder-division-rate", type=float, default=None,
                     help="Synthetic_SMC (neointimal founder) CREC_1=1 division rate, 1/min. Defaults to the SAME "
                          "value as --division-rate (Smooth_Muscle_Cell's own) unless explicitly overridden -- an "
                          "earlier mismatched default (1e-5, ~48-day doubling, vs --division-rate's 2e-4, ~2.4-day) "
                          "meant that once most agents transformed to Synthetic_SMC, the tissue's effective "
                          "division capacity collapsed ~20x, emptying out the tissue over months even without "
                          "AMLO. Both phenotypes now divide at the same rate by default, consistent with there "
                          "being no literature basis for synthetic SMC proliferating slower than contractile SMC "
                          "in situ (if anything the opposite is reported) -- only the BOOLEAN-gated cycle-entry "
                          "logic (which pathway_node drives it) differs between them.")
    ap.add_argument("--smc-apoptosis-rate", type=float, default=5.31667e-06,
                     help="BASE apoptosis (death model 100) rate, 1/min, for BOTH Smooth_Muscle_Cell and "
                          "Synthetic_SMC -- what a cell with Resist_apop=0 experiences. History: originally "
                          "5.31667e-05 (PhysiCell's standard template default -- a cancer-cell-line apoptosis "
                          "rate, ~13-day mean lifetime; biologically fast for differentiated vascular SMC, whose "
                          "healthy in-vivo turnover is on the order of months). That rate was tolerable only "
                          "because Hif1a was ON in ~99.6% of SMC-lineage cells (near-universal Resist_apop "
                          "protection masked it); once the O2 gradient fix (see --o2-dirichlet-value) made "
                          "hypoxia partial (~50%) instead of near-total, the ~50% of cells left unprotected "
                          "died faster than the population-cap-gated division rate could replace them, AND -- "
                          "because fewer live cells consume less O2 -- that death raised ambient O2 further, "
                          "turning off MORE cells' Hif1a and accelerating the decline (SMC count fell from 1765 "
                          "to 1090, i.e. -20 to -34/day, over days 20-45 of a 45-day test at the old rate). "
                          "Lowered 10x to 5.31667e-06 (~130-day mean lifetime, closer to real differentiated-SMC "
                          "turnover) so unprotected cells no longer die faster than gated division can replace "
                          "them. Re-validate population trajectory after any future change to the pressure/"
                          "population gate or O2 parameters, since this balance is sensitive to both. v8: "
                          "apoptosis resistance is driven continuously per-agent by the Resist_apop Boolean node "
                          "(Nfat OR Klf4), via the 'apoptosis' output mapping.")
    ap.add_argument("--transformation-rate", type=float, default=1.0667e-5,
                     help="Smooth_Muscle_Cell -> Synthetic_SMC transformation rate, 1/min, applied when Klf4=1 "
                          "(v8's contractile->synthetic toggle). Default slowed ~15x from an earlier 1.6e-4 "
                          "(~3-day half-life once Klf4=1): with Klf4 itself reaching ~90%+ quickly near the "
                          "ET1-rich lumen (chronic ET1 secretion from Endothelial_Cell), the original rate "
                          "transformed ~91% of the SMC lineage within 21 days even with NO drug at all -- far too "
                          "fast for a model meant to show disease progression over months. 1.0667e-5 gives a "
                          "~45-day half-life once Klf4 is ON, spreading transformation over the 3-6 month "
                          "campaign timescale instead of saturating in 3 weeks -- still UNCALIBRATED (no "
                          "published rate exists), but chosen to make the long-run time axis meaningful.")
    ap.add_argument("--reverse-transformation-rate", type=float, default=1.0667e-6,
                     help="Synthetic_SMC -> Smooth_Muscle_Cell transformation rate, 1/min, applied when Myocd=1 "
                          "(reversion to contractile). Kept at 10x slower than --transformation-rate (same ratio "
                          "as before, both now slowed ~15x) -- UNCALIBRATED.")
    ap.add_argument("--founder-migration-speed", type=float, default=0.05,
                     help="Synthetic_SMC CREC_1=1 migration speed (micron/min), directed up the ET1 gradient "
                          "(native PhysiCell <chemotaxis>, already enabled for this cell type -- see the "
                          "is_founder branch). Default raised from 0.0 (which left chemotaxis structurally "
                          "enabled but functionally inert, since direction without speed is a no-op) to a "
                          "modest 0.05 um/min -- deliberately much smaller than the invasive migration speeds "
                          "explored in this project's dedicated neointima-formation experiments (section 10.9-10.10 "
                          "of the docs), since this long-run campaign is not specifically testing lumen invasion.")
    ap.add_argument("--migration-bias", type=float, default=0.8,
                     help="directedness (0=random walk, 1=fully directed) of founder-cell migration toward the lumen")
    ap.add_argument("--ec-barrier-adhesion", type=float, default=4.0,
                     help="Endothelial_Cell cell-cell adhesion strength. Raised from an earlier 0.8: at that value "
                          "the EC ring deformed into a chaotic scatter under sustained pressure from the "
                          "proliferating media (individual EC agents end up at wildly different radii -- e.g. "
                          "41.9 to 174.3 um by day 180 in one long run -- rather than compressing together as a "
                          "cohesive sheet, which is what real endothelial cell-cell junctions would do). 4.0 keeps "
                          "neighboring EC agents glued together so the ring can still be compressed inward as a "
                          "whole (the intended lumen-narrowing readout) without scattering individual cells "
                          "through gaps in it.")
    ap.add_argument("--ec-barrier-repulsion", type=float, default=100.0,
                     help="Endothelial_Cell cell-cell repulsion strength. Raised from an earlier 20.0 alongside "
                          "--smc-repulsion-strength -- see that argument's help for the severe-overlap finding "
                          "that motivated both increases (--lumen-radius clamp still provides the actual lumen "
                          "impenetrability; this is about the ring's own mechanical integrity, not the barrier).")
    ap.add_argument("--smc-adhesion-strength", type=float, default=0.4,
                     help="Smooth_Muscle_Cell/Synthetic_SMC cell-cell adhesion strength -- kept at the original "
                          "PhysiCell-template default (unlike repulsion below, adhesion wasn't implicated in the "
                          "overlap finding).")
    ap.add_argument("--dt-mechanics", type=float, default=0.1,
                     help="Mechanics timestep, min. History: the original PhysiCell-template default of 1 min "
                          "crashed (core dump within 2 frames) once --smc-repulsion-strength/--ec-barrier-"
                          "repulsion were raised to fix the overlap problem (see that argument's help). "
                          "0.1 was tried first; a dt=0.5 attempt (whole-run crash check only) looked OK but "
                          "turned out to silently break EC-ring cohesion specifically (no crash, but the ring's "
                          "radial std grew from ~1 to ~128 over 45 days -- a numerical-accuracy failure, not an "
                          "instability the solver flags). Reducing --ec-barrier-repulsion from 100 to 40 (even "
                          "with --fill-lumen-core's rigid backstop) only delayed the onset (day 1 to day 4), "
                          "not prevented it. 0.1 is restored as the default -- the EC-EC repulsion constant "
                          "(the largest force constant in the system) needs this fine a timestep for ACCURATE "
                          "integration, not just a non-crashing one. The --n-threads default (8, up from 4) and "
                          "--fill-lumen-core (an additional, no-tradeoff cohesion improvement at this same dt) "
                          "are the validated speed/quality wins instead -- see their help text. If you need to "
                          "push past the 0.1 default again, validate with the EC-ring r_std metric specifically "
                          "over a multi-week run, not just a crash check, and not just the cheap sandbox (its "
                          "cleaner initial layout masks instabilities the real founder geometry exposes).")
    ap.add_argument("--smc-repulsion-strength", type=float, default=60.0,
                     help="Smooth_Muscle_Cell/Synthetic_SMC cell-cell repulsion strength. Raised sharply from an "
                          "earlier 10.0 (the original PhysiCell-template default): measured directly in a "
                          "completed 180-day run, essentially every nearby cell pair was overlapping by day 180 "
                          "(80,849 of 80,856 pairs checked within 15 um, mean overlap 42% of the pair's combined "
                          "radii, 30,821 pairs over 50% overlap) -- the repulsion force was far too weak to keep "
                          "up with how fast new cells were being added by division across a dense, ~2,500-agent "
                          "tissue. 60.0 (6x) is a first correction, not independently re-validated at the full "
                          "180-day duration -- re-check overlap statistics on any new long run before trusting "
                          "this value fully.")
    ap.add_argument("--f-snx29", type=float, default=1.0,
                     help="SNX29 brake level on Stim1_pm trafficking (0.5 low/1.0 normal/2.0 high, Karnes 2026 GWAS, section 8)")
    ap.add_argument("--escape-rate-scale", type=float, default=1.0,
                     help="Multiplier on $u_stim/$d_stim/$d_stim_off (default 1.0 = original, fast: the AMLO-escape "
                          "compensatory Stim1_pm pathway reaches quasi-steady-state within ~1 simulated day regardless "
                          "of run length). A value <1 (e.g. 0.02-0.1) stretches this pathway's own characteristic "
                          "timescale without touching electrical-speed ($k_fast) or protein-turnover ($u_cav/$d_cav) "
                          "rates elsewhere in the network -- pair with a longer --max-time so the slowed transient "
                          "actually unfolds within the observed window instead of plateauing immediately. Uncalibrated "
                          "-- no published rate exists for this specific compensatory pathway's real-time kinetics.")
    ap.add_argument("--lumen-radius", type=float, default=100.0,
                     help="Lumen radius, um -- founder cells are hard-clamped to stay >= this (basement membrane/IEL, cannot be mechanically crossed)")
    ap.add_argument("--fill-lumen-core", type=lambda s: s.lower() not in ("0", "false", "no"), default=True,
                     help="If true (default), fill the lumen interior (r < --lumen-radius) with a dense, "
                          "chemically-inert 'Lumen_Core' cell population -- a mechanical anchor for the "
                          "Endothelial_Cell ring's inner face, so it isn't trying to stay cohesive against empty "
                          "space. Soft/movable by default (see --lumen-core-rigid) so sustained SMC-lineage "
                          "pressure can still compress the EC ring inward (functional/structural lumen "
                          "narrowing, a real readout this model already uses) instead of hitting a wall that "
                          "can never yield.")
    ap.add_argument("--lumen-core-rigid", type=lambda s: s.lower() not in ("0", "false", "no"), default=False,
                     help="If true, make Lumen_Core cells perfectly immovable (is_movable=false) -- maximum EC "
                          "ring cohesion, but blocks lumen narrowing entirely (the ring can never push an "
                          "immovable floor out of its way). Default false: Lumen_Core stays movable, relying on "
                          "--lumen-core-repulsion-strength being weak enough to yield under sustained pressure "
                          "while still being a presence the EC ring can lean on day-to-day.")
    ap.add_argument("--lumen-core-diameter", type=float, default=12.0, help="Lumen_Core cell diameter, um (matches EC).")
    ap.add_argument("--lumen-core-packing-fraction", type=float, default=0.95)
    ap.add_argument("--lumen-core-adhesion-strength", type=float, default=0.4)
    ap.add_argument("--lumen-core-repulsion-strength", type=float, default=10.0,
                     help="Lumen_Core cell-cell/cell-EC repulsion strength. History: 60.0 (matched to "
                          "--smc-repulsion-strength) when the core was rigid (--lumen-core-rigid=true) -- fixed "
                          "EC scattering but also blocked lumen narrowing, since a rigid floor with strong "
                          "repulsion is doubly immovable. Lowered to 10.0 (the project's original, pre-fix SMC "
                          "repulsion value) now that the core is soft/movable by default: a WEAK repulsion "
                          "lets sustained SMC-lineage pressure actually displace/compress Lumen_Core cells "
                          "toward the center rather than having them resist like a wall.")
    ap.add_argument("--contraction-speed", type=float, default=0.5,
                     help="um/min inward radial velocity added to Smooth_Muscle_Cell/Synthetic_SMC agents while CONTR_1 (myosin phosphorylation/MLCK-MLCP "
                          "balance) is ON -- a reversible, functional vasoconstriction proxy distinct from the irreversible structural narrowing driven by "
                          "Synthetic_SMC accumulation. Removed (not reversed) when CONTR_1 is OFF, letting existing repulsion/pressure mechanics relax the "
                          "ring back outward. Uncalibrated placeholder -- no published PhysiCell-scale SMC contraction velocity to fit against.")
    ap.add_argument("--founder-pathway-node", type=str, default="CREC_1",
                     help="Boolean node gating founder cell division + migration (\"cycle entry\"/\"migration speed\" outputs). "
                          "Default CREC_1 (= Nfat OR Hif1a) lets self-generated hypoxia in the packed subendothelial mass dominate "
                          "and swamp any LCCB/SNX29/NFAT signal -- set to \"Nfat\" to isolate the calcium/calcineurin/NFAT pathway "
                          "from the hypoxia confound and cleanly test the amlodipine/SNX29 hypothesis.")
    ap.add_argument("--et1-secretion-rate", type=float, default=0.02, help="feedback arm: Nfat=1 -> ET1 secretion rate")
    ap.add_argument("--o2-uptake-rate", type=float, default=20.0, help="feedback arm: Ca_global=1 -> O2 uptake rate")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--per-agent-heterogeneity", action="store_true",
                     help="Sample each Smooth_Muscle_Cell/Synthetic_SMC agent's own STIM1/ORAI1/ATP2A2/TRPC6 "
                          "level from the real snRNA-seq-fitted PAH Beta distributions (REAL_GENE_BETA_PARAMS_PAH), "
                          "instead of every agent sharing one population-level istate. Off by default so existing "
                          "results remain exactly reproducible. CACNA1C is excluded (v8 made it Myocd-derived).")
    ap.add_argument("--amlo-periodic-dosing", action="store_true",
                     help="Replace AMLO's constant EC secretion_rate with a once-daily, literature-grounded oral PK "
                          "profile (one-compartment absorption/elimination, redosed every --amlo-dose-interval-min). "
                          "Off by default -- constant secretion_rate (the original behavior) is unchanged.")
    ap.add_argument("--amlo-half-life-hours", type=float, default=40.0,
                     help="Amlodipine terminal elimination half-life, hours (real range 30-50h; PMID 16296697, "
                          "FDA label) -- sets the PK elimination rate constant ke=ln(2)/t_half.")
    ap.add_argument("--amlo-tmax-hours", type=float, default=8.0,
                     help="Time to peak plasma concentration after a dose, hours (real range 6-12h) -- sets the "
                          "PK absorption rate constant ka (solved numerically given ke).")
    ap.add_argument("--amlo-dose-interval-min", type=float, default=1440.0,
                     help="Redosing interval, minutes (1440 = once daily).")
    ap.add_argument("--o2-baseline-uptake-rate", type=float, default=8.0,
                     help="Smooth_Muscle_Cell/Synthetic_SMC baseline O2 uptake rate, 1/min. History: originally "
                          "10/min, crashed O2 near-zero everywhere under the old Dirichlet-0 outer wall (see "
                          "--o2-dirichlet-value) and was lowered to 1.5/min as a workaround -- but that produced a "
                          "monotonic radial decline with NO recovery, driving Hif1a ON in 99.6% of SMC-lineage "
                          "cells by day 180 of the long campaign (56% of all CREC_1 activations came via Hif1a "
                          "alone, swamping the intended Ca2+/calcineurin/NFAT/amlodipine pathway). Restored to "
                          "8/min now that --o2-dirichlet-value is a real nonzero floor (not 0): combined, these "
                          "give a genuine U-shaped profile (~30 mmHg near lumen -> ~15 mmHg dip mid-media -> ~23 "
                          "mmHg recovering at the outer wall), validated to put ~50% of the SMC-lineage population "
                          "below the HIPOXIA threshold -- a real but partial, not total, hypoxic core.")
    ap.add_argument("--o2-fibro-uptake-rate", type=float, default=4.0,
                     help="Fibroblast baseline O2 uptake rate, 1/min. Kept at half the SMC rate, as before; see "
                          "--o2-baseline-uptake-rate for the rationale behind the current value.")
    ap.add_argument("--o2-vasa-vasorum-secretion-rate", type=float, default=4.0,
                     help="Fibroblast O2 secretion rate, 1/min -- a distributed vasa-vasorum proxy (see the "
                          "code comment at fibro_secretion). 0.0 disables it, recovering the old single-source "
                          "(EC-only) diffusion picture. Tune jointly with --o2-vasa-vasorum-target and "
                          "--o2-baseline-uptake-rate; re-check the achieved hypoxic fraction at STEADY-STATE "
                          "density (day 40+ of a real run), not just the low-density initial condition -- the "
                          "whole point of a distributed source tied to the fixed-count Fibroblast population is "
                          "to be more robust to SMC density drift than the old single Dirichlet-wall value was, "
                          "but it still needs validating at the density that actually matters.")
    ap.add_argument("--o2-vasa-vasorum-target", type=float, default=25.0,
                     help="Fibroblast O2 secretion_target, mmHg -- the local O2 level the vasa-vasorum proxy "
                          "pulls the adventitia toward. Set below the EC's 38 mmHg secretion_target (a "
                          "real vasa vasorum is a smaller, secondary supply, not full arterial-blood "
                          "oxygenation).")
    ap.add_argument("--amlo-dose-mg", type=float, default=10.0,
                     help="Simulated oral dose, mg (clinically 5/10/20mg once daily). Assumes linear PK (amlodipine "
                          "is a low-clearance, dose-proportional drug at these doses): secretion rate at steady-"
                          "state peak scales as dose_mg/10 relative to --amlo-secretion-rate, which calibrates the "
                          "10mg reference dose.")
    args = ap.parse_args()

    # Axial (z) layer detection from the actual cells.csv -- single source of truth
    # for the 2D-vs-3D switch; --z-layers, if given, is only a consistency check.
    args._z_planes = detect_z_layers(args.cells_csv)
    n_z = len(args._z_planes)
    if args.z_layers is not None and args.z_layers != n_z:
        raise SystemExit(f"--z-layers={args.z_layers} disagrees with the cells.csv "
                         f"({n_z} distinct z planes detected: {[t for _v, t in args._z_planes]})")
    if n_z > 1:
        print(f"3D z-extrusion detected: {n_z} axial layer planes "
              f"(z = {[t for _v, t in args._z_planes]}) -> use_2D=false, z-sized domain, "
              f"per-layer Lumen_Core fill")

    if args.population_half_max is None:
        frac = (args.wedge_angle_degrees / 360.0) if args.geometry == "wedge" else 1.0
        args.population_half_max = 2000.0 * frac * n_z
        print(f"--population-half-max not set: resolved to {args.population_half_max:.1f} "
              f"(geometry={args.geometry}, frac={frac:.4g}, z_layers={n_z})")

    # --- Amlodipine oral PK parameter solving (pure Python, no scipy dependency) ---
    def _solve_amlo_pk_ka(t_half_hours: float, tmax_hours: float):
        ke = math.log(2) / (t_half_hours * 60.0)
        tmax_min = tmax_hours * 60.0
        def tmax_pred(ka):
            return math.log(ka / ke) / (ka - ke)
        lo, hi = ke * 1.001, 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if tmax_pred(mid) > tmax_min:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi), ke

    def _amlo_pk_steady_state_peak(ka: float, ke: float, dose_interval: float, n_doses: int = 30) -> float:
        """Steady-state PEAK of the superposed one-compartment PK curve for
        dose_amount=1 (the --amlo-dose-mg=10 reference), by brute-force
        superposition over n_doses -- used to normalize EC secretion_rate
        so the 10mg reference dose exactly recovers --amlo-secretion-rate
        at its steady-state peak, with other doses scaling proportionally."""
        best = 0.0
        t0 = (n_doses - 1) * dose_interval
        for i in range(2000):
            t = t0 + i * dose_interval / 2000.0
            total = 0.0
            for n in range(n_doses):
                td = t - n * dose_interval
                if td < 0:
                    continue
                total += ka / (ka - ke) * (math.exp(-ke * td) - math.exp(-ka * td))
            best = max(best, total)
        return best

    if args.founder_division_rate is None:
        args.founder_division_rate = args.division_rate

    args.amlo_pk_ka, args.amlo_pk_ke = _solve_amlo_pk_ka(args.amlo_half_life_hours, args.amlo_tmax_hours)
    args.amlo_pk_ss_peak = _amlo_pk_steady_state_peak(args.amlo_pk_ka, args.amlo_pk_ke, args.amlo_dose_interval_min)
    args.amlo_pk_dose_amount = args.amlo_dose_mg / 10.0

    out = args.out_dir
    (out / "config").mkdir(parents=True, exist_ok=True)
    (out / "custom_modules").mkdir(parents=True, exist_ok=True)

    (out / "config" / "PhysiCell_settings.xml").write_text(build_settings_xml(args))
    if args.per_agent_heterogeneity:
        augment_cells_csv_with_heterogeneity(args.cells_csv, out / "config" / "cells.csv", args.seed)
    else:
        shutil.copy(args.cells_csv, out / "config" / "cells.csv")
    if args.fill_lumen_core:
        lc_theta_max = math.radians(args.wedge_angle_degrees) if args.geometry == "wedge" else 2 * math.pi
        append_lumen_core_rows(out / "config" / "cells.csv", args.lumen_radius, args.lumen_core_diameter,
                                args.lumen_core_packing_fraction, args.seed, theta_min=0.0, theta_max=lc_theta_max)

    models_dir = HERE.parent / "models"
    shutil.copy(models_dir / "ModeloPASMC_v8.bnd", out / "config")
    cfg_text = (models_dir / "ModeloPASMC_v8.cfg").read_text()
    cfg_text = re.sub(
        r"\$f_snx29\s*=\s*[0-9.]+\s*;",
        f"$f_snx29   = {args.f_snx29};",
        cfg_text, count=1,
    )
    assert f"$f_snx29   = {args.f_snx29};" in cfg_text, "f_snx29 override patch failed"
    if args.escape_rate_scale != 1.0:
        # Stretch the SNX29-braked Stim1_pm compensatory ("escape") pathway's
        # own characteristic timescale, leaving every other node's rates
        # (electrical/fast-process $k_fast, CACNA1C protein turnover $u_cav/
        # $d_cav, Soce's own gating) untouched. At escape_rate_scale=1 (the
        # original $u_stim=$d_stim=1, $d_stim_off=5), the whole network
        # relaxes to a quasi-steady-state within ~1 simulated day regardless
        # of simulated duration -- confirmed empirically (pct_stenosis/
        # CONTR_1_frac/Nfat_frac all plateau by ~1440 min and stay flat for
        # the remaining 6 days of a 7-day run), which makes an "acute day-1
        # vs chronic day-7" comparison just two points on the same plateau,
        # not a genuine slow transient. Scaling these three rates down
        # stretches Stim1_pm's (and everything causally downstream of it:
        # Soce, CONTR_1, Nfat) rise/fall time proportionally, matched with a
        # correspondingly longer --max-time so the transient actually
        # unfolds within the observation window instead of snapping to
        # plateau immediately.
        for param, base in [("u_stim", 1.0), ("d_stim", 1.0), ("d_stim_off", 5.0)]:
            new_val = base * args.escape_rate_scale
            cfg_text, n = re.subn(
                rf"\${param}\s*=\s*[0-9.]+\s*;",
                f"${param}    = {new_val};",
                cfg_text, count=1,
            )
            assert n == 1, f"{param} override patch failed"
    if args.per_agent_heterogeneity:
        # Make the personalized nodes visible in boolean_intracellular.csv so
        # the per-agent heterogeneity effect can actually be observed/plotted
        # (they're is_internal=TRUE, i.e. hidden, in the base v8 cfg).
        for gene in REAL_GENE_BETA_PARAMS_PAH:
            cfg_text = re.sub(rf"{gene}\.is_internal\s*=\s*TRUE", f"{gene}.is_internal = FALSE", cfg_text)
    (out / "config" / "ModeloPASMC_v8.cfg").write_text(cfg_text)

    sample_dir = args.physicell_root / "sample_projects_intracellular" / "boolean" / "physiboss_cell_lines"
    custom_h_text = (sample_dir / "custom_modules" / "custom.h").read_text()
    custom_h_text = custom_h_text.replace(
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt );",
        "void founder_barrier_velocity( Cell* pCell, Phenotype& phenotype, double dt );\n"
        "void apply_contraction_pull( Cell* pCell, double dt );\n"
        "void apply_lumen_barrier( Cell* pCell, double dt );\n"
        "void contraction_velocity( Cell* pCell, Phenotype& phenotype, double dt );\n"
        "void update_amlodipine_pk( double current_time, double dt );\n"
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt );",
        1,
    )
    assert "founder_barrier_velocity" in custom_h_text, "custom.h barrier forward-decl patch failed"
    assert "contraction_velocity" in custom_h_text, "custom.h contraction forward-decl patch failed"
    (out / "custom_modules" / "custom.h").write_text(custom_h_text)
    custom_cpp_text = (sample_dir / "custom_modules" / "custom.cpp").read_text()
    custom_cpp_text = re.sub(
        r"void pre_update_intracellular\( Cell\* pCell, Phenotype& phenotype, double dt \)\r?\n\{.*?\r?\n\}",
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )\n"
        "{\n\t// no-op: sample's $time_scale desync logic doesn't apply to this model\n}",
        custom_cpp_text, count=1, flags=re.DOTALL,
    )
    # The physiboss_cell_lines sample template this project is built on
    # disables continuous mechanics/motility velocity updates entirely
    # (update_velocity = NULL) -- fine for a non-dividing, non-moving
    # sample population, but it silently breaks BOTH adhesion/repulsion
    # relaxation after division AND any motility (migration_bias/speed)
    # for this project's growth/migration experiments, since neither
    # pathway runs without it. Re-enable the standard PhysiCell mechanics
    # update (already linked in via PhysiCell_standard_models.h).
    custom_cpp_text = custom_cpp_text.replace(
        "cell_defaults.functions.update_velocity = NULL;",
        "cell_defaults.functions.update_velocity = standard_update_cell_velocity;",
        1,
    )
    assert "standard_update_cell_velocity;" in custom_cpp_text, "update_velocity patch failed"
    # Lumen_Core (see --fill-lumen-core): a dense, inert disk filling the lumen interior
    # (r < lumen_radius) to give the Endothelial_Cell ring a mechanical anchor from the INSIDE,
    # instead of empty space. Without ANY core, the EC ring's only resistance to being dented/
    # scattered by SMC-lineage pressure from the outside is its own self-adhesion/repulsion -- which
    # needed extreme values (and a correspondingly fine dt_mechanics, costing ~5x wall-clock time) to
    # stay cohesive.
    #
    # --lumen-core-rigid (default FALSE, i.e. soft/movable): a perfectly rigid core (is_movable=false)
    # fixed EC scattering but has a worse failure mode -- it also makes functional/structural LUMEN
    # NARROWING impossible, since the EC ring can never push a rigid, immovable floor out of the way.
    # Lumen narrowing (the EC ring compressing inward under sustained SMC-lineage pressure) is a real,
    # already-used readout in this model (the 5th-percentile EC radial position). Default to a SOFT
    # core instead: is_movable stays at PhysiCell's normal default (true), and
    # --lumen-core-repulsion-strength is set much weaker than the EC ring's own (see that argument),
    # so Lumen_Core cells get pushed toward the center and compressed rather than acting as a floor --
    # giving the EC ring a cohesion-supporting neighbor on its inner face that still YIELDS to sustained
    # pressure, instead of blocking the biology the mechanics fix was never meant to touch.
    if args.lumen_core_rigid:
        custom_cpp_text = re.sub(
            r"display_cell_definitions\( std::cout \); ?\r?\n",
            "display_cell_definitions( std::cout ); \n"
            "\n"
            "\tfor ( int i=0; i < cell_definitions_by_index.size(); i++ )\n"
            "\t{\n"
            "\t\tif ( cell_definitions_by_index[i]->name == \"Lumen_Core\" )\n"
            "\t\t{ cell_definitions_by_index[i]->is_movable = false; }\n"
            "\t}\n",
            custom_cpp_text, count=1,
        )
        assert 'cell_definitions_by_index[i]->name == "Lumen_Core"' in custom_cpp_text, "Lumen_Core is_movable patch failed"
    # Fibroblast/Endothelial_Cell have NO <intracellular> block in this
    # project (phenotype.intracellular == NULL for them) -- the sample's
    # coloring functions call intracellular methods unconditionally on
    # every cell and segfault on the bystander types. Guard both.
    custom_cpp_text = re.sub(
        r"if \( !pCell->phenotype\.intracellular->get_boolean_variable_value\( parameters\.strings\(\"node_to_visualize\"\) \) \)",
        'if ( pCell->phenotype.intracellular == NULL ) { return output; }\n'
        '\tif ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )',
        custom_cpp_text, count=1,
    )
    custom_cpp_text = re.sub(
        r"void color_node\(Cell\* pCell\)\{\r?\n\tstd::string node_name = parameters\.strings\(\"node_to_visualize\"\);\r?\n\tpCell->custom_data\[node_name\] = pCell->phenotype\.intracellular->get_boolean_variable_value\(node_name\);\r?\n\}",
        "void color_node(Cell* pCell){\n"
        "\tif ( pCell->phenotype.intracellular == NULL ) { return; }\n"
        "\tstd::string node_name = parameters.strings(\"node_to_visualize\");\n"
        "\tpCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);\n"
        "}",
        custom_cpp_text, count=1,
    )
    assert "pCell->phenotype.intracellular == NULL ) { return output; }" in custom_cpp_text, "my_coloring_function patch failed"
    assert "pCell->phenotype.intracellular == NULL ) { return; }" in custom_cpp_text, "color_node patch failed"

    # Pressure-gated contact inhibition on "cycle entry": PhysiBoSS's own
    # CREC_1 -> "cycle entry" <output> mapping sets the division rate every
    # intracellular_dt from the Boolean state alone, with no awareness of
    # local crowding (the "live" cycle model has no built-in contact
    # inhibition). Left unchecked, exponential Boolean-driven growth
    # outruns what cell-cell repulsion mechanics can equilibrate between
    # divisions, producing non-physical overlapping clumps instead of
    # outward tissue growth. This adds a Hill-type suppression of
    # "cycle entry" by PhysiCell's native simple_pressure signal (the
    # standard contact-inhibition idiom for growth models in PhysiCell),
    # applied once per intracellular_dt right after PhysiBoSS sets the
    # Boolean-driven rate, so growth concentrates at the low-density
    # growing edge instead of an already-packed interior.
    custom_cpp_text = custom_cpp_text.replace(
        "void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )\n{\n\tcolor_node(pCell);\n}",
        "void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )\n"
        "{\n"
        "\tcolor_node(pCell);\n"
        "\tif ( pCell->phenotype.intracellular == NULL ) { return; }\n"
        "\tstatic double p_half = parameters.doubles(\"pressure_half_max\");\n"
        "\tstatic double p_hill = parameters.doubles(\"pressure_hill_power\");\n"
        "\tdouble pressure = pCell->state.simple_pressure;\n"
        "\tpCell->custom_data[\"pressure_debug\"] = pressure;\n"
        "\tdouble ratio_n = pow( pressure / p_half, p_hill );\n"
        "\tdouble gate = 1.0 / (1.0 + ratio_n); // Hill-type suppression: 1 at zero pressure, 0.5 at p_half, ->0 when crowded\n"
        # GLOBAL population cap: local simple_pressure alone proved unable to
        # reliably bound growth over multi-week timescales (tested at several
        # pressure_half_max/pressure_hill_power combinations -- a continuously
        # expanding low-pressure tissue edge kept outrunning local suppression,
        # reaching ~10,000-24,000 agents within 25-65 days regardless). This
        # adds a second, GLOBAL Hill-type gate on total agent count (the same
        # functional form, just reading *all_cells->size() instead of local
        # pressure) as a robust backstop -- multiplicative with the local
        # pressure gate, so local contact inhibition still shapes WHERE growth
        # concentrates (the low-density edge), while this caps HOW MANY cells
        # can exist in total, independent of whatever is happening with local
        # mechanics.
        "\tstatic double n_half = parameters.doubles(\"population_half_max\");\n"
        "\tstatic double n_hill = parameters.doubles(\"population_hill_power\");\n"
        "\tdouble n_total = (double) (*all_cells).size();\n"
        "\tdouble n_ratio_n = pow( n_total / n_half, n_hill );\n"
        "\tdouble n_gate = 1.0 / (1.0 + n_ratio_n);\n"
        "\tdouble current_rate = PhysiCell::get_single_behavior(pCell, \"cycle entry\");\n"
        "\tPhysiCell::set_single_behavior(pCell, \"cycle entry\", current_rate * gate * n_gate);\n"
        "}"
    )
    assert "pressure_half_max" in custom_cpp_text, "pressure-gating patch failed"
    assert "population_half_max" in custom_cpp_text, "population-cap patch failed"

    # (Neointimal founder-cell DIRECTION is handled by PhysiCell's NATIVE
    # <chemotaxis> mechanism -- up the ET1 gradient toward the endothelial
    # source. But soft mechanics (EC-EC adhesion/repulsion alone) was
    # tested and found insufficient to keep the endothelium intact under
    # sustained invasive pressure -- the EC ring itself got dragged apart
    # and founder cells still ended up deep in the lumen. Neointima
    # formation happens in the SUBENDOTHELIAL space, below an intact
    # endothelium/internal elastic lamina that normal mechanical pressure
    # cannot breach (crossing it in vivo requires active proteolytic
    # degradation, not simple crowding) -- so this is modeled as a hard
    # positional constraint: founder cells are clamped to r >= lumen_radius
    # every mechanics step, regardless of how much inward force/chemotaxis
    # is pushing on them.
    custom_cpp_text = custom_cpp_text.replace(
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )",
        "void apply_contraction_pull( Cell* pCell, double dt )\n"
        "{\n"
        "\t// Reversible functional vasoconstriction: while the Boolean node\n"
        "\t// CONTR_1 (myosin phosphorylation / MLCK-MLCP balance, requires\n"
        "\t// Myocd -- i.e. only contractile-phenotype SMCs can do this) is ON,\n"
        "\t// add a small inward radial velocity pulling the agent toward the\n"
        "\t// vessel centerline. This is a VELOCITY term, not a position clamp:\n"
        "\t// when CONTR_1 turns back OFF (e.g. amlodipine blocking Ltcc), no\n"
        "\t// inward term is added and the existing repulsion/pressure\n"
        "\t// mechanics relax the ring back outward on their own -- genuine\n"
        "\t// vasorelaxation, distinct from the irreversible structural\n"
        "\t// narrowing driven by Synthetic_SMC accumulation.\n"
        "\tif ( pCell->phenotype.intracellular == NULL ) { return; }\n"
        "\tstatic double contraction_speed = parameters.doubles(\"contraction_speed\");\n"
        "\tif ( !pCell->phenotype.intracellular->get_boolean_variable_value(\"CONTR_1\") ) { return; }\n"
        "\tdouble r = sqrt( pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1] );\n"
        "\tif (r > 1e-6) {\n"
        "\t\tdouble nx = pCell->position[0]/r, ny = pCell->position[1]/r;\n"
        "\t\tpCell->velocity[0] -= contraction_speed*nx;\n"
        "\t\tpCell->velocity[1] -= contraction_speed*ny;\n"
        "\t}\n"
        "}\n"
        "\n"
        "void apply_lumen_barrier( Cell* pCell, double dt )\n"
        "{\n"
        "\t// Basement membrane / internal elastic lamina as an impenetrable\n"
        "\t// barrier: no SMC-lineage agent (contractile OR synthetic/founder)\n"
        "\t// may ever occupy the flowing-blood lumen, whether pushed there by\n"
        "\t// neointimal migration or pulled there by its own CONTR_1-driven\n"
        "\t// contraction -- crossing an intact endothelium/IEL in vivo needs\n"
        "\t// active proteolytic degradation, not simple crowding or myosin\n"
        "\t// contraction. The Endothelial_Cell ring itself is NOT subject to\n"
        "\t// this clamp -- it is free to be mechanically compressed inward\n"
        "\t// (below the nominal lumen_radius) by pressure from the SMC-lineage\n"
        "\t// mass outside it, which is exactly how structural/functional\n"
        "\t// lumen narrowing is read out (5th-percentile EC radial position).\n"
        "\t// Implemented as a VELOCITY correction, never a position teleport --\n"
        "\t// an earlier version called Cell::assign_position() every mechanics\n"
        "\t// step to clamp position directly, which calls register_agent()\n"
        "\t// (just a list push_back) without removing the cell from its\n"
        "\t// previous voxel, silently duplicating its neighbor-list\n"
        "\t// registration on every call and causing a runaway force/position\n"
        "\t// blowup within 1-2 simulated days. Cancelling the inward velocity\n"
        "\t// component here instead lets position/voxel bookkeeping go through\n"
        "\t// PhysiCell\'s own (correct) update_position() pathway, same as\n"
        "\t// every other cell.\n"
        "\tstatic double lumen_r = parameters.doubles(\"lumen_radius\");\n"
        "\tdouble r = sqrt( pCell->position[0]*pCell->position[0] + pCell->position[1]*pCell->position[1] );\n"
        "\tif (r < lumen_r && r > 1e-6) {\n"
        "\t\tdouble nx = pCell->position[0]/r, ny = pCell->position[1]/r;\n"
        "\t\tdouble radial_v = pCell->velocity[0]*nx + pCell->velocity[1]*ny;\n"
        "\t\tif (radial_v < 0) { pCell->velocity[0] -= radial_v*nx; pCell->velocity[1] -= radial_v*ny; }\n"
        "\t\tdouble penetration = lumen_r - r;\n"
        "\t\tdouble push = penetration * 0.5;\n"
        "\t\tpCell->velocity[0] += push*nx; pCell->velocity[1] += push*ny;\n"
        "\t}\n"
        "}\n"
        "\n"
        "void contraction_velocity( Cell* pCell, Phenotype& phenotype, double dt )\n"
        "{\n"
        "\t// Wired to Smooth_Muscle_Cell. Same three-stage composition as\n"
        "\t// founder_barrier_velocity below: standard mechanics, then the\n"
        "\t// CONTR_1-driven inward pull, then the impenetrable-lumen clamp\n"
        "\t// (which always has the final say -- contraction cannot push a\n"
        "\t// contractile SMC through the endothelium any more than migration\n"
        "\t// can push a founder cell through it).\n"
        "\tstandard_update_cell_velocity( pCell, phenotype, dt );\n"
        "\tapply_contraction_pull( pCell, dt );\n"
        "\tapply_lumen_barrier( pCell, dt );\n"
        "}\n"
        "\n"
        "void founder_barrier_velocity( Cell* pCell, Phenotype& phenotype, double dt )\n"
        "{\n"
        "\t// Normal mechanics first (adhesion/repulsion/EC contact, etc.) --\n"
        "\t// this is the SAME function cell_defaults.functions.update_velocity\n"
        "\t// already uses, just called explicitly so we can add a hard\n"
        "\t// constraint on top for this one cell type.\n"
        "\tstandard_update_cell_velocity( pCell, phenotype, dt );\n"
        "\n"
        "\t// Reversible functional vasoconstriction (CONTR_1), same as\n"
        "\t// Smooth_Muscle_Cell's contraction_velocity -- applied BEFORE the\n"
        "\t// lumen clamp below, so the clamp still has the final say if\n"
        "\t// contraction would otherwise pull this cell past the boundary.\n"
        "\tapply_contraction_pull( pCell, dt );\n"
        "\tapply_lumen_barrier( pCell, dt );\n"
        "}\n\n"
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )",
        1,
    )
    custom_cpp_text = custom_cpp_text.replace(
        "\tbuild_cell_definitions_maps(); ",
        "\tCell_Definition* pFounderCD2 = find_cell_definition(\"Synthetic_SMC\");\n"
        "\tif (pFounderCD2 != NULL) { pFounderCD2->functions.update_velocity = founder_barrier_velocity; }\n\n"
        "\tCell_Definition* pSMC_CD = find_cell_definition(\"Smooth_Muscle_Cell\");\n"
        "\tif (pSMC_CD != NULL) { pSMC_CD->functions.update_velocity = contraction_velocity; }\n\n"
        "\tbuild_cell_definitions_maps(); ",
        1,
    )
    assert "founder_barrier_velocity" in custom_cpp_text, "barrier function patch failed"
    assert "pFounderCD2" in custom_cpp_text, "barrier wiring failed"
    assert "pSMC_CD" in custom_cpp_text, "contraction wiring failed"

    # One-compartment oral-absorption PK model for amlodipine, redosed every
    # --amlo-dose-interval-min (default 1440 min = 24h). ka/ke solved offline
    # in Python (see build_settings_xml's PK block below) from real
    # amlodipine pharmacokinetics (terminal half-life 30-50h -> ke; Tmax
    # 6-12h -> ka; PMID 16296697 and FDA label) and passed in as
    # user_parameters, along with ss_peak (the theoretical steady-state PEAK
    # of the C_level ODE at dose_amount=1, also solved offline) so that
    # EC secretion_rate = base_rate * (C_level / ss_peak) recovers the
    # original constant-dosing calibration at peak. Disabled (falls back to
    # the original constant XML secretion_rate) unless amlo_periodic_dosing
    # is true. Only touches Endothelial_Cell agents (the AMLO source,
    # cached once since ECs neither divide nor die in this model).
    custom_cpp_text = custom_cpp_text.replace(
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )",
        "void update_amlodipine_pk( double current_time, double dt )\n"
        "{\n"
        "\tstatic bool periodic_dosing = parameters.bools(\"amlo_periodic_dosing\");\n"
        "\tif ( !periodic_dosing ) { return; }\n"
        "\n"
        "\tstatic double ka = parameters.doubles(\"amlo_pk_ka\");\n"
        "\tstatic double ke = parameters.doubles(\"amlo_pk_ke\");\n"
        "\tstatic double dose_interval = parameters.doubles(\"amlo_dose_interval_min\");\n"
        "\tstatic double dose_amount = parameters.doubles(\"amlo_pk_dose_amount\");\n"
        "\tstatic double ss_peak = parameters.doubles(\"amlo_pk_ss_peak\");\n"
        "\tstatic double base_rate = parameters.doubles(\"amlo_pk_base_secretion_rate\");\n"
        "\n"
        "\tstatic double A_gut = 0.0;\n"
        "\tstatic double C_level = 0.0;\n"
        "\tstatic double next_dose_time = 0.0;\n"
        "\tstatic std::vector<Cell*> ec_cells;\n"
        "\tstatic bool ec_cells_found = false;\n"
        "\tstatic int amlo_index = microenvironment.find_density_index(\"AMLO\");\n"
        "\n"
        "\tif ( !ec_cells_found )\n"
        "\t{\n"
        "\t\tfor ( auto pC : *all_cells ) { if ( pC->type_name == \"Endothelial_Cell\" ) { ec_cells.push_back(pC); } }\n"
        "\t\tec_cells_found = true;\n"
        "\t}\n"
        "\n"
        "\tif ( current_time >= next_dose_time - 1e-6 )\n"
        "\t{\n"
        "\t\tA_gut += dose_amount;\n"
        "\t\tnext_dose_time += dose_interval;\n"
        "\t}\n"
        "\n"
        "\tdouble dA = -ka * A_gut;\n"
        "\tdouble dC = ka * A_gut - ke * C_level;\n"
        "\tA_gut += dA * dt;\n"
        "\tC_level += dC * dt;\n"
        "\tif ( A_gut < 0 ) { A_gut = 0; }\n"
        "\tif ( C_level < 0 ) { C_level = 0; }\n"
        "\n"
        "\tdouble secretion_rate = base_rate * (C_level / ss_peak);\n"
        "\tfor ( auto pC : ec_cells ) { pC->phenotype.secretion.secretion_rates[amlo_index] = secretion_rate; }\n"
        "}\n"
        "\n"
        "void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )",
        1,
    )
    assert "update_amlodipine_pk" in custom_cpp_text, "amlodipine PK function patch failed"
    (out / "custom_modules" / "custom.cpp").write_text(custom_cpp_text)

    main_cpp_text = (sample_dir / "main.cpp").read_text()
    main_cpp_text, n_main_patches = re.subn(
        r"(\t+/\*\r?\n\t+ *Custom add-ons could potentially go here\. *\r?\n\t+\*/\r?\n)",
        r"\1\r\n\t\t\tupdate_amlodipine_pk( PhysiCell_globals.current_time, diffusion_dt );\r\n",
        main_cpp_text, count=1,
    )
    assert n_main_patches == 1, "main.cpp custom-addon insertion point not found"
    (out / "main.cpp").write_text(main_cpp_text)
    makefile_text = (sample_dir / "Makefile").read_text()
    makefile_text = makefile_text.replace("PhysiBoSS_Cell_Lines", "PhysiBoSS_Vessel_Gradient")
    (out / "Makefile").write_text(makefile_text)

    print(f"Wrote PhysiBoSS vessel-gradient project to {out} (feedback={args.feedback})")


if __name__ == "__main__":
    main()
