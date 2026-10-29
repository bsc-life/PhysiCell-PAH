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

The root `Makefile`'s `list-projects` target and project-loading rules
gained one new target, `pah-vessel`, following the exact
convention already used by the other `sample_projects_intracellular/boolean/*`
projects (`physiboss-cell-lines-sample`, `physiboss-tutorial`, etc.) — copy
the project's `config/`/`custom_modules/`/`main.cpp`/`Makefile` over the
root versions, then `make`.

## Pushing this to your own GitHub fork

```bash
git remote -v                                  # shows "upstream" -> MathCancer/PhysiCell
git remote add origin git@github.com:<you>/PhysiCell.git   # your fork, created via GitHub's Fork button
git push -u origin master
```
