# PAH_vessel — a PhysiBoSS vessel-wall model of pulmonary arterial hypertension

**HAP** = *Hipertensión Arterial Pulmonar* (Spanish for Pulmonary Arterial
Hypertension, PAH), the name used in this project's own source cell-type
inventory report. This sample project is a spatially-resolved PhysiCell +
PhysiBoSS agent-based model of the pulmonary arterial vessel wall, built
around a Ca2+ / Calcineurin / NFAT / Kv1.5 positive-feedback loop in
pulmonary arterial smooth muscle cells (PASMC), extended with a
Myocd/Klf4 contractile<->synthetic phenotype toggle and an apoptosis-
resistance node. It runs each PhysiCell agent's own embedded MaBoSS
Boolean network via PhysiBoSS's native `<input>`/`<output>` mapping — no
custom per-timestep C++ Boolean-network glue beyond PhysiCell's standard
mechanics/transformation/apoptosis hooks.

This project folder is adapted, with the raw multi-terabyte run outputs
and R&D history stripped out, from a larger research repository
(`pah_digital_twin`); see `docs/amlodipine_adaptation_mechanism.md` for
the full methods/results narrative this model was developed against
(sections 10–22 cover the spatial vessel-wall build specifically).

## What's modeled

* **Geometry**: a vessel cross-section (lumen / media / adventitia rings),
  PAH-thickened media with a minority `Synthetic_SMC` founder subpopulation
  that can invade the lumen (neointima formation) via ET1-gradient
  chemotaxis, held back by an impenetrable endothelial barrier.
* **Intracellular model**: `config/ModeloPASMC_v8.bnd`/`.cfg`, a
  SIGNOR/NeKo-curated Boolean network covering the electrical branch
  (Kv1_5, Em_depol, Bkca, Ltcc/CaV1.2), the ER/SOCE branch (Plc, Serca,
  Stim1_pm, Soce), receptor-operated Ca2+ entry (Trpc6, Pkc), a hypoxia
  sensor (Hif1a), Ca2+-sensitization/contraction (Rock, Mlcp, Mlck,
  CONTR_1/2), a Myocd/Klf4 contractile<->synthetic toggle, and a
  `Resist_apop` node.
* **Substrates**: four diffusing fields (AMLO, ET1, O2, K_ALTO) feed
  Boolean inputs via PhysiBoSS's native `<input>` mapping; `Nfat` ->
  ET1 secretion and `Ca_global` -> O2 uptake close two feedback loops onto
  the microenvironment via the native `<output>` mapping. AMLO uses a
  literature-grounded once-daily oral PK dosing model (one-compartment,
  `amlo_periodic_dosing` flag in `PhysiCell_settings.xml`).
* **Per-agent heterogeneity**: each SMC agent's STIM1/ORAI1/ATP2A2/TRPC6
  level is drawn from real snRNA-seq-fitted distributions via
  `custom:*` columns in `config/cells.csv`, rather than one deterministic
  istate shared by the whole population.
* **Mechanics**: `CONTR_1` drives a reversible inward radial contraction
  pull (vasorelaxation/vasoconstriction), and cell volumes/mechanics are
  set to literature diameters (SMC/Fibroblast 15 um, EC 12 um).

## Build & run

```bash
# from the repo root (after `make PhysiBoSS_Vessel_Gradient-setup` style
# libMaBoSS fetch happens automatically on first build):
make pah-vessel
make -j4
./PhysiBoSS_Vessel_Gradient
```

The default `config/` is the **amlo_snx29_normal** arm (chronic AMLO
dosing on, SNX29 level normal). Two alternative scenario arms are provided
in `config/variants/` — copy their `PhysiCell_settings.xml`/
`ModeloPASMC_v8.cfg` over the defaults before building to reproduce them:

| Variant | `amlo_periodic_dosing` | `$f_snx29` |
|---|---|---|
| `config/` (default) | true | 1.0 (normal) |
| `config/variants/no_amlo` | false | 1.0 (normal) |
| `config/variants/amlo_snx29_low` | true | 0.5 (low) |

`custom_modules/custom.cpp`/`.h` and `main.cpp` are identical across all
three arms — only these two config values differ.

## Reproducing / regenerating the project

`scripts/generate_physiboss_vessel_gradient_config.py` is the generator
that produced this project's `config/`+`custom_modules/`+`main.cpp`+
`Makefile` from the vessel geometry + Boolean model + CLI flags (AMLO
dosing, SNX29 level, per-agent heterogeneity, apoptosis rates, etc.).
`scripts/build_vessel_cross_section.py` builds the idealized vessel
cross-section `cells.csv` this project's geometry is based on.
`scripts/personalize_intracellular_model.py` personalizes Boolean istates
from snRNA-seq pseudobulk expression (not run against real patient data in
this project — see `docs/amlodipine_adaptation_mechanism.md`).

## Figures

Headline results from the full research history are in `docs/figures/`:
`v8_transformation_validation.png` (Boolean-driven SMC->Synthetic_SMC
phenotype transformation), `longrun_campaign_summary.png` /
`longrun_6month_phenotype.gif` (6-month radial transformation front),
`snx29_amlo_clinical_metrics.png` (stenosis/PVR-proxy clinical metrics),
`spatial_vessel_gradient_feedback.png` (radial substrate/Boolean-state
gradients), `physiboss_ecology_clonal_sweep.gif` (clonal-selection
PhysiBoSS ecology arm).

## Attribution

Built on [PhysiCell](https://github.com/MathCancer/PhysiCell) (Ghaffarizadeh
et al., PLoS Comput Biol 2018) and its official
[PhysiBoSS](https://github.com/PhysiBoSS/PhysiBoSS) addon (Letort et al.,
Bioinformatics 2019; Ponce-de-Leon et al., 2022), using
[MaBoSS](https://maboss.curie.fr/) for the embedded Boolean simulation.
Model edge provenance is curated against [SIGNOR](https://signor.uniroma2.it/)
via NeKo. See `docs/amlodipine_adaptation_mechanism.md` for literature
citations behind individual mechanisms (AMLO escape, SNX29, TRPC6/NFAT
feedback, etc.).
