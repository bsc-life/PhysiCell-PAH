# The PAH vessel-wall model, in this fork

**HAP** = *Hipertensión Arterial Pulmonar* (Spanish for Pulmonary Arterial
Hypertension, PAH) — the condition this fork's added model simulates. This
page is a short map of what's new; the full model description, build
instructions, and provenance live in the project folder itself.

## What the model does

A spatially-resolved PhysiCell + PhysiBoSS agent-based model of the
pulmonary arterial vessel wall, centered on a Ca2+ / Calcineurin / NFAT /
Kv1.5 positive-feedback loop in pulmonary arterial smooth muscle cells
(PASMC). Each agent runs its own embedded MaBoSS Boolean network — a
SIGNOR/NeKo-curated model extended with a Myocd/Klf4
contractile<->synthetic phenotype toggle and an apoptosis-resistance node —
wired to four diffusing substrates (AMLO, ET1, O2, K_ALTO) through
PhysiBoSS's native `<input>`/`<output>` mapping. The model reproduces
PAH-characteristic neointima formation (a Boolean-gated founder-cell
subpopulation invading the lumen via ET1 chemotaxis) and lets you compare
chronic amlodipine dosing against two SNX29 genetic-brake levels.

## Where everything is

```
sample_projects_intracellular/boolean/PAH_vessel/
├── README.md               # full model description, build & run instructions
├── config/                  # default build: AMLO dosing on, SNX29 normal
│   └── variants/             # two alternative scenario arms (no_amlo, amlo_snx29_low)
├── custom_modules/          # custom.cpp/.h (mechanics, transformations, feedback)
├── main.cpp, Makefile
├── docs/
│   ├── amlodipine_adaptation_mechanism.md   # methods/results writeup
│   ├── gene_node_mapping.yaml               # snRNA-seq gene -> Boolean node mapping
│   └── figures/              # headline result figures/GIFs
└── scripts/                  # the 3 generator scripts used to build this project
```

## Building it

```bash
make pah-vessel
make -j4
./PhysiBoSS_Vessel_Gradient
```

See [`sample_projects_intracellular/boolean/PAH_vessel/README.md`](sample_projects_intracellular/boolean/PAH_vessel/README.md)
for the full writeup: the Boolean network's node-by-node description, the
three scenario variants and what each changes, per-agent heterogeneity,
the mechanics/contraction model, and attribution (PhysiCell, PhysiBoSS,
MaBoSS, SIGNOR).
