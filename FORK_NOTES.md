# This is a fork of PhysiCell

This repository is a fork of [MathCancer/PhysiCell](https://github.com/MathCancer/PhysiCell)
(cloned at tag/version `1.14.2`, upstream history preserved as the
`upstream` remote). Everything under the usual PhysiCell top-level
layout (`core/`, `BioFVM/`, `addons/`, `sample_projects*/`, the root
`Makefile`, etc.) is unmodified upstream PhysiCell.

## What this fork adds

One new sample project:
[`sample_projects_intracellular/boolean/PAH_vessel/`](sample_projects_intracellular/boolean/PAH_vessel/README.md) —
a PhysiBoSS agent-based model of the pulmonary arterial vessel wall in PAH
(**HAP** = *Hipertensión Arterial Pulmonar*), built on a curated PASMC
Boolean network (Ca2+/Calcineurin/NFAT/Kv1.5 loop + Myocd/Klf4 phenotype
toggle). See that folder's own README for the model description, build
instructions, and provenance.
