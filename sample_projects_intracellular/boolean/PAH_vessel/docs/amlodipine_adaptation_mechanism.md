# Amlodipine adaptation: a store-independent STIM1/ORAI1 escape route

**Question:** can PASMC adapt to chronic amlodipine (L-type Ca²⁺ channel
blocker, LCCB) exposure over time, and if so, which molecular mechanism
allows it?

## 1. The literature mechanism

Amlodipine blocks CaV1.2 (the `Ltcc` node), which is the model's intended
and only direct drug target. Johnson MT et al., *"L-type Ca2+ channel
blockers promote vascular remodeling through activation of STIM proteins"*,
PNAS 2020 (DOI: [10.1073/pnas.2007598117](https://doi.org/10.1073/pnas.2007598117))
reported that chronic LCCB exposure also has an *unintended* effect: it
directly promotes STIM1/Orai1 puncta assembly and store-operated Ca²⁺ entry
(SOCE) **independent of CaV1.2 channel activity and independent of ER/SR
Ca²⁺ store content** -- i.e. not the canonical store-depletion trigger, a
separate route. Their proposed biochemical mechanism has two parts: (1)
LCCBs (all three major chemical classes -- dihydropyridines like
amlodipine, phenylalkylamines, benzothiazepines) directly bind a 10-amino-acid
N-terminal region of STIM1 located in the ER lumen -- the same region
involved in STIM1's native Ca²⁺-sensing -- pharmacologically mimicking the
"store depleted" signal; and (2) STIM1 is independently known to physically
bind CaV1.2's C-terminus at rest, partially sequestering a pool of STIM1 in
an inactive state (Park et al., Science 2010) -- Johnson et al. propose
LCCBs also disrupt this resting STIM1-CaV1.2 complex, "unlocking" that
sequestered STIM1 pool into its active, ORAI-activating conformation. Even
the original authors note isolated-membrane experiments failed to reproduce
the effect, implying cells need some unidentified "intermediary factor" --
i.e. the chain is not fully resolved even in the source paper. This sits
alongside the separately well-documented clinical observation that only a
minority of PAH patients (reported ~26% in acute vasoreactivity testing)
are CCB "responders" at all.

**IMPORTANT, project-critical caveat: this literature mechanism is
scientifically contested, not settled.** A specialist CRAC-channel
laboratory (Bird, Lin, Tucker, Mueller, Shi, Padmanabhan, Parekh --
*"Scrutinizing science to save lives: uncovering flaws in the data linking
L-type calcium channel blockers to CRAC channels and heart failure,"*
bioRxiv 2024, doi:10.1101/2024.02.06.579229; preceded by *"A reappraisal of
the effects of L-type Ca2+ channel blockers on store-operated Ca2+ entry
and heart failure,"* Function 2023, doi:10.1093/function/zqad047; followed
by a 2025 companion paper) directly disputes Johnson et al. 2020's central
claims:

- They argue amlodipine is **intrinsically fluorescent inside cells and
  overwhelms the fura-2 calcium-indicator dye signal** that Johnson et
  al.'s calcium-imaging conclusions were largely built on -- a
  methodological confound, not a real SOCE signal.
- They cite a substantial prior literature showing **micromolar LCCB
  concentrations *block* store-operated Ca²⁺ entry rather than activate
  it** -- the literal opposite direction from Johnson et al.'s claim.
- Their own meta-analysis of published clinical trials and a prospective
  real-world cohort found dihydropyridines (amlodipine's class) were **not
  associated with increased heart failure or other cardiovascular
  events** -- directly contradicting Johnson et al.'s clinical-risk
  framing (which warned against LCCB use in elderly/advanced-hypertension
  patients).
- Their stated conclusion: *"We therefore believe the study by Johnson et
  al. is flawed at many levels."*

**What this means for everything built on this document:** sections 2-8
below treat the Johnson et al. 2020 mechanism as a literature-grounded
hypothesis worth modeling -- which is a defensible choice for a *mechanistic
exploration* (and the model's own caveats throughout already flag
`$u_amlo_escape` as an uncalibrated placeholder, not a validated rate) --
but the premise itself (whether chronic LCCB exposure activates STIM1/SOCE
at all, let alone at the magnitude modeled here) is an open, actively
disputed scientific question, not an established fact. Every downstream
result in this document (escape magnitude, real-heterogeneity bottleneck,
the ORAI1+ clonal-selection hypothesis, the SNX29-bypass structural
argument) should be read as *"if the Johnson et al. mechanism is real, then
this model predicts..."*, not as a claim about real patients. This caveat
was not in earlier versions of this document and should have been --
corrected here after the controversy surfaced while answering a follow-up
question about the mechanism's biochemical plausibility.

## 2. Model extension

`models/ModeloPASMC_base_v2.bnd` (the base model) is left **unmodified**.
`models/ModeloPASMC_base_v2_AMLO_escape.bnd` adds exactly one edge: the
`Stim1_pm` node gains a second, store-independent activating route gated by
`AMLO`:

```
logic     = (STIM1 AND (NOT Sr_lleno)) OR (STIM1 AND AMLO)
rate_up   = @logic ? ( (NOT Sr_lleno) ? <original SOCE rate> : $u_amlo_escape ) : 0.0
```

When the SR store is depleted (`NOT Sr_lleno`), behavior is identical to the
base model regardless of AMLO. When the store is full but `AMLO` is on, the
new branch alone can still drive `Stim1_pm` assembly, at a new,
**uncalibrated** rate `$u_amlo_escape = 0.5` (swept, not fit — there is no
published Boolean transition rate to calibrate against; 0.5 was placed
between `$k_transc` (0.1) and `$u_stim` (1) as a "moderate" placeholder).

## 3. Simulation design

16 MaBoSS runs (`scripts/run_amlo_adaptation_sweep.py`): {base, escape} model
variant × AMLO {off, on} × ET1 {off, on} × HIPOXIA {off, on}, `max_time=600`,
`time_tick=4`, `sample_count=1000`, all other istates/parameters at the
shipped resting-cell defaults. Each AMLO-input node has `rate_up =
rate_down = 0` (fixed for the whole run), so "acute" vs "chronic" response
is read off one continuous trajectory's early vs late time points rather
than from separate pre/post-dose runs.

## 4. Result

As a correctness check: with `AMLO = 0`, the base and escape models are
**numerically identical** at every time point in every (ET1, HIPOXIA)
combination (the new edge is inert unless AMLO is on) — confirming the
extension doesn't silently change baseline behavior.

With `AMLO = 1` and no other stimulus ("amlodipine monotherapy" — the
cleanest test of the mechanism in isolation), steady-state (t≈600)
marginal probabilities:

| Node | Base model | + AMLO escape edge | Interpretation |
|---|---|---|---|
| `Ltcc` (CaV1.2 open) | 0.008 | 0.026 | Drug target equally blocked in both — amlodipine "works" either way |
| `Soce` | 0.042 | 0.340 (**8.1×**) | Store-operated entry rises sharply despite Ltcc blockade |
| `Stim1_pm` | 0.043 | 0.360 (**8.4×**) | STIM1 puncta assemble via the new store-independent route |
| `Nfat` | 0.165 | 0.555 (**3.4×**) | Growth-transcription factor stays substantially active |
| `CONTR_1` (contraction) | 0.108 | 0.176 (**1.6×**) | Some residual contractile signal, smaller effect than on growth |
| `CREC_1` (growth/remodeling) | 0.090 | 0.292 (**3.3×**) | Remodeling signal rebounds despite CaV1.2 blockade |

With `AMLO = 1` **and** `ET1 = 1` (amlodipine prescribed against persistent
vasoconstrictor/disease tone), the canonical store-depletion SOCE route is
already partly engaged by ET1 in the BASE model (`Soce` 0.293, `Nfat`
0.580) — so the escape edge's *additional* contribution is smaller in
relative terms (`Soce` → 0.449, **+53%**) than in monotherapy, though still
present. **The escape mechanism's impact is largest when it is the only
available Ca²⁺-entry route (monotherapy) and proportionally smaller — but
not absent — when disease-associated tone already engages the canonical
pathway.** Full 16-condition table: `amlo_adaptation_summary_early_late.csv`;
full time series: `amlo_adaptation_timecourse.csv`; figure:
`amlo_adaptation_timecourse.png`.

## 5. Caveats

- `$u_amlo_escape` is a placeholder, not a calibrated value — the sweep
  demonstrates the *qualitative* mechanism (escape is possible and
  mechanistically plausible per the cited literature) and its *direction*
  (growth/remodeling escapes more than contraction), not a quantitative
  prediction of its magnitude in a real patient.
- This is a single, long MaBoSS trajectory per condition, not a true
  time-varying dosing protocol (start drug-naive, add drug partway through,
  continue). A sequential two-phase run (feeding phase-1's final marginal
  probabilities back in as phase-2's istates) would be a closer model of an
  actual treatment timeline and is a natural next step.
- Not yet connected to PhysiBoSS/PhysiCell — this is a MaBoSS-only,
  well-mixed (non-spatial) analysis of the intracellular model's own
  dynamics. Wiring a diffusing amlodipine substrate to the `AMLO` node and
  `CONTR_1/2`/`CREC_1/2` to agent phenotype in PhysiCell is scoped but not
  built (see README §5's "Next step").

## 6. Does real SMC expression heterogeneity change the escape magnitude?

The §2–4 sweep above assumed all four personalizable nodes (`CACNA1C`,
`STIM1`, `ORAI1`, `ATP2A2`) start constitutively ON (`istate = 1`) — an
idealized, fully-expressing cell. `smc_heterogeneity_summary.csv` (real
snRNA-seq, Control n=1,205 SMCs / PAH n=3,510 SMCs) lets us instead
personalize each node's istate with the real fraction of SMC nuclei that
detectably express that gene, and ask whether Control-vs-PAH expression
differences change the predicted escape effect.
(`run_amlo_escape_real_data_personalized.py`; full result table:
`amlo_escape_real_data_personalized.csv`; figure:
`amlo_escape_real_heterogeneity.png`.)

**Istate mapping.** `personalize_intracellular_model.py`'s default pairwise
min-max scaling saturates to the floor/ceiling (0.02/0.98) for *any*
two-group comparison, regardless of true fold-change — a 1.06× STIM1
difference would be rendered identically to a 1.31× CACNA1C difference. We
used `pct_expressing / 100` directly (clipped to [0.02, 0.98]) instead,
which preserves real effect-size magnitude:

| Gene | Control istate P(1) | PAH istate P(1) | Fold-change |
|---|---|---|---|
| CACNA1C | 0.418 | 0.466 | 1.11× |
| STIM1 | 0.180 | 0.205 | 1.14× |
| ORAI1 | 0.032 | 0.020 (floor-clipped) | 0.64× (**down** in PAH) |
| ATP2A2 | 0.281 | 0.322 | 1.15× |

**Finding: ORAI1 is a hard bottleneck that the real data make visible, not
one the escape mechanism can route around.** `Soce { logic = Stim1_pm AND
ORAI1 }` is a strict AND-gate in *both* the base and AMLO-escape model
variants — the literature edge added for §1 only changes `Stim1_pm`'s
dynamics, never `Soce`'s dependency on `ORAI1`. Because only ~2–3% of real
SMC nuclei detectably express *ORAI1* mRNA (in either condition), the
population-level `Soce`, `CONTR_2`, `CREC_1` marginals stay below 1% under
AMLO monotherapy regardless of Control/PAH personalization — roughly two
orders of magnitude smaller than the idealized §2–4 sweep (which showed an
8× escape effect assuming every cell carries full ORAI1). `Stim1_pm` itself
(upstream of the gate) does show the expected, modest personalization
effect: the escape edge raises it by +13% in Control vs +18% in PAH under
AMLO monotherapy (consistent with PAH's higher STIM1 fraction), and by only
+3%/+7% respectively once ET1 is also present (confirming §4's "smaller but
not absent when canonical SOCE is already engaged").

**Interpretation.** Taking the real per-gene detection fraction literally
as "this fraction of the whole SMC population can execute the pathway"
implies the amlodipine-escape route is not a uniform SMC property but one
gated by a minority *ORAI1*-expressing subpopulation — and, by this
transcript-level metric, that subpopulation is if anything *smaller* in PAH
than in Control, which runs counter to a naive "PAH SMCs are uniformly more
escape-prone" story. Two caveats temper this: (1) *ORAI1* is a small,
lowly-expressed channel subunit, and snRNA-seq dropout is known to
under-detect such transcripts in both conditions — the true
protein-level-expressing fraction is almost certainly higher than 2–3% in
both groups, so the absolute magnitude here should not be over-read; and
(2) treating pct_expressing as a *frozen* per-trajectory istate (this node
is either fully on or fully off for the whole simulated run, since
`CACNA1C`/`STIM1`/`ORAI1`/`ATP2A2` all have `rate_up = rate_down = 0` in the
model) is itself a binary-population-split modeling choice, not a graded
dose-response — it tests "what fraction of cells can ever engage this arm"
rather than "how strongly does an average cell engage it." A model that let
expression grade the node's own transition *rate* (rather than only its
initial condition) would be a more faithful next step once mapped against
real single-cell RNA counts rather than a condition-level summary.

A second, unrelated bug was caught while building this comparison and is
now fixed in `personalize_intracellular_model.py`: MaBoSS's probabilistic
istate grammar requires square brackets around the node name
(`[NODE].istate = p0 [0], p1 [1];`); the bare `NODE.istate = ...` form used
previously is a configuration syntax error that the original string-level
unit test did not catch (it checked the generated text, but no cfg was
ever actually run through MaBoSS). All cfgs in this document and in
`run_amlo_escape_real_data_personalized.py` use the corrected form.

## 7. Ecological hypothesis: does a rare ORAI1+ subclone outcompete wild-type under chronic AMLO?

Section 6 showed that the AMLO-escape mechanism, taken literally as a
population-averaged istate probability, is bottlenecked by *ORAI1*'s low
(~2-3%) real detection fraction. That framing treats the whole SMC
population as uniformly, partially escape-capable. The alternative,
testable hypothesis: what if *ORAI1* expression is not a population-wide
dial but a stable, binary trait carried by a small, pre-existing
subclone -- and chronic amlodipine exposure gives that subclone a
proliferative advantage large enough to let it take over the tissue over
weeks, purely through selection (no new mutation, no transcriptional
switching)? If so, "acquired" amlodipine resistance would be a *clonal
composition shift*, not a per-cell adaptation -- and the real heterogeneity
data (ORAI1 detected in only ~2-3% of SMCs) is exactly the kind of standing
variation such a sweep could start from.

### 7.1 Model: two genotypes, one Boolean network, one PhysiBoSS build

Built and tested the project's first full PhysiBoSS integration (PhysiCell
agents each running their own MaBoSS instance of
`ModeloPASMC_base_v2_AMLO_escape.bnd/.cfg`), using the official
`addons/PhysiBoSS` module (libMaBoSS v2.5.6, statically linked) shipped
inside the fetched PhysiCell 1.14.2 source tree.

Two cell definitions, `SMC_ORAI1pos` and `SMC_ORAI1neg`, run the *identical*
Boolean network and differ in exactly one istate override: `ORAI1`=1 vs 0.
Because `ORAI1` has `rate_up = rate_down = 0` in the .bnd (a frozen,
input-like node, like `AMLO`/`ET1`/`HIPOXIA`), this is a stable, heritable
trait -- daughters are created as the same cell definition as their mother
(PhysiCell's standard division behavior), so the genotype clones with no
extra bookkeeping. `AMLO` is driven by a diffusible microenvironment
substrate (PhysiBoSS native `<input>` mapping, action=activation,
threshold=0.5); `CREC_1` (the model's own growth/remodeling readout) drives
PhysiCell's native "cycle entry" behavior (native `<output>` mapping) --
a cell's division rate is read directly off its own Boolean network's
`CREC_1` duty cycle, not a hand-tuned phenotype parameter.
`scripts/generate_physiboss_ecology_config.py` generates the full,
buildable project (XML + cells.csv + custom_modules + Makefile) from these
two parameters plus a population seed fraction.

**A real build/runtime bug found and fixed along the way:** the officially
shipped `pre_update_intracellular()` template (copied from
`sample_projects_intracellular/boolean/physiboss_cell_lines`) looks up a
`$time_scale` MaBoSS parameter that only exists in *that sample's own*
`.cfg` file. Calling `get_parameter_value()` on an undefined symbol
dereferences a null pointer inside the embedded engine
(`Symbol::getIndex()`, `BooleanNetwork.h`) and segfaults -- confirmed via
an AddressSanitizer backtrace (`SEGV ... Symbol::getIndex() ... called from
custom.cpp:164`), isolated by writing a minimal standalone harness linking
directly against `libMaBoSS-static.a` (bypassing PhysiCell entirely) to
rule out the Boolean model itself before suspecting the C++ wiring.
Fixed by stripping the sample-specific body; `generate_physiboss_ecology_config.py`
patches this automatically in every project it generates.

### 7.2 Three arms, two run to completion

Seeded 500 agents in a 1x1 mm 2D patch, 3% `SMC_ORAI1pos` / 97%
`SMC_ORAI1neg` (matching the real snRNA-seq detection fraction), no cell
death modeled (pure growth-competition; composition, not absolute count,
is the readout):

- **chronic_amlo**: AMLO substrate fixed ON (Dirichlet boundary = 1) from
  t=0, 30 simulated days.
- **no_drug**: AMLO substrate OFF throughout, 30 simulated days -- the
  neutral-drift counterfactual.
- **delayed_amlo** (designed, not run): drug-naive growth phase then AMLO
  switched on partway through. The CLI flags exist
  (`--arm delayed_amlo` / `--dirichlet-override-on`) but carrying phase-1's
  per-agent Boolean states into a phase-2 run needs a checkpoint-restart
  mechanism not implemented in this PhysiCell version's `main.cpp` without
  extra custom read-in code -- deferred rather than rushed, since the two
  completed arms already answer the core question unambiguously.

### 7.3 Result: the sweep happens in both arms -- AMLO makes it faster and more complete

| | chronic_amlo | no_drug |
|---|---|---|
| Final agent count (day 30) | 26,921 | 2,218 |
| ORAI1+ fraction, day 30 | 98.1% | 77.0% |
| Day to reach 50% ORAI1+ | 17 | 25 |
| Population-mean CREC_1 duty cycle, day 30 | 40.9% | 22.3% |

Both arms show a genuine clonal sweep from the 3% starting fraction --
*even without any drug*, ORAI1+ cells have an intrinsic proliferative edge
in this model and parameterization, because `ORAI1-negative` cells sit at
essentially 0% CREC_1 duty cycle throughout *both* conditions (the model's
resting attractor is proliferation-locked without the channel), while
ORAI1+ cells show a non-trivial duty cycle (30-49%) from day 1 onward.
Chronic AMLO exposure does not single-handedly create this advantage, but
substantially **amplifies and accelerates** it: the 50% crossover arrives
8 days sooner, and the final sweep is both faster-saturating and more
complete (98.1% vs 77.0% ORAI1+ by day 30) -- consistent with the
literature mechanism's own framing (Johnson et al., PNAS 2020) that chronic
LCCB exposure drives a maladaptive, CaV1.2-independent growth signal
specifically through the STIM1/ORAI1 route. Practically: by this model, a
tissue initially 97% "amlodipine-responsive" (ORAI1-negative) SMCs would, by
day 30 of continuous treatment, be >98% composed of cells that never
needed to adapt at all -- they already carried the escape trait, and
selection did the rest.

### 7.4 Caveats

- **No cell death is modeled.** Final agent counts (up to 26,921 for a
  500-cell, 1x1mm 2D seed) reflect unconstrained growth, not a biologically
  realistic carrying capacity -- read the *fraction* and *duty-cycle*
  columns, not the absolute counts.
- **The cycle-entry value/base_value mapping (division rate when CREC_1=1
  vs 0) is a hand-chosen, uncalibrated design choice** (24h doubling when
  "on", ~quiescent when "off") -- there is no published rate linking this
  Boolean readout to an actual PASMC doubling time. The qualitative result
  (which genotype wins, and that AMLO accelerates it) is far more robust
  than the specific day-17-vs-day-25 numbers.
- **The "even without drug, ORAI1+ already has an edge" finding is itself
  parameter-sensitive** -- it reflects where this particular model's
  resting attractor sits relative to the CREC_1 threshold, not an
  independently validated fact about real PASMC. Treat it as "this model,
  as built, predicts..." rather than a settled biological claim.
- **2D, well-mixed seeding**, not the vessel-wall ring geometry from
  `build_vessel_cross_section.py` -- a tractability simplification for a
  compositional question, not a spatial one.
- Builds on section 1's model and section 6's real-data personalization;
  all of section 6's caveats about `$u_amlo_escape` being swept-not-fit
  still apply here, since the same escape edge underlies `Stim1_pm`'s
  AMLO-driven activation in this per-agent model.

Scripts: `generate_physiboss_ecology_config.py` (project generator).
Figure: `physiboss_ecology_clonal_sweep.png`. Data:
`physiboss_ecology_timecourse.csv`, `physiboss_ecology_summary.csv`. Full
project bundles + build/run logs: `ecology_bundles.tar.gz`.

## 8. The SNX29 variant: a germline brake on the canonical route that AMLO-escape bypasses

`ModeloPASMC_base_v2(.bnd)`'s own comments cite "SNX29 -| translocacion de
STIM1 (manual, Karnes 2026)" -- this is a real, identified citation:
**Karnes JH, Bao C, Liang S, et al. "Variation in SNX29 and Acute
Vasodilator Response in Pulmonary Arterial Hypertension." Circ Res.
2026;139(2):e327117.** (doi:10.1161/CIRCRESAHA.125.327117, PMC13186119;
companion editorial: Yu PB, PMC13336244). The study is a GWAS of **acute
vasodilator response** in PAH -- the clinical test that identifies the
minority "vasodilator-responsive" subgroup eligible for long-term CCB
therapy (the same ~26%-responder population motivating section 1's whole
premise). A SNX29 variant, rs8057488, reached genome-wide significance for
this response and was nominally replicated. Functionally: SNX29 localizes
to endosomes in PASMC; *silencing* it relocates STIM1/Orai1 to the plasma
membrane and **enhances** SOCE, while *overexpressing* it **attenuates**
hypoxic vasoconstriction in perfused lung -- i.e. SNX29 is a brake on
STIM1/Orai1-mediated SOCE, and that brake's strength plausibly determines
whether a patient's vasculature can relax enough, acutely, to count as a
"responder."

**This is exactly the `$f_snx29` parameter already in the model, not a new
mechanism.** `Stim1_pm`'s canonical rate is `$u_stim * (...) / $f_snx29`
(the `.cfg` ships three illustrative levels: `0.5` low / `1.0` normal /
`2.0` high SNX29) -- dividing by a larger `$f_snx29` slows `Stim1_pm`
assembly, i.e. **more SNX29 = stronger brake**, matching Karnes 2026's
overexpression/silencing results in direction. A 6-condition MaBoSS sweep
({low, normal, high} SNX29 × {AMLO off, on}, steady state, monotherapy)
confirms this quantitatively and reveals the structural point this project
had not yet made explicit:

| SNX29 level | Soce, no AMLO | Soce, chronic AMLO | CREC_1, no AMLO | CREC_1, chronic AMLO |
|---|---|---|---|---|
| low (0.5×) | 0.0615 | 0.378 | 0.158 | 0.293 |
| normal (1.0×) | 0.0419 | 0.338 | 0.112 | 0.286 |
| high (2.0×) | 0.0217 | 0.312 | 0.070 | 0.284 |
| **spread (low/high)** | **2.8×** | **1.21×** | **2.3×** | **1.03×** |

**Without AMLO, SNX29 level is the dominant lever** on `Soce`/`CREC_1`
(2.3-2.8× spread across the three levels) -- the brake works as designed,
and as Karnes 2026's functional data would predict: a patient/cell with
weak SNX29 braking sits closer to a SOCE-active, less vasodilation-capable
state even before any drug. **Under chronic AMLO, that same brake barely
matters** (CREC_1 spread collapses to 1.03×, Soce to 1.21×): the escape
edge (`$u_amlo_escape`, section 1) is written in the `.bnd` as a *separate*
`OR` branch on `Stim1_pm` that is **not divided by `$f_snx29`** -- a
structural choice (store-independence was the literature's explicit claim;
`$f_snx29`-independence was not separately asserted by Johnson et al., but
falls out of giving the escape branch its own, ungated rate) that this
sweep makes visible rather than just assumed from reading the formula.

**Why this matters for the rest of this document.** CCB therapy is, by
design, prescribed to the SNX29-favorable responder subgroup (robust
braking, intact acute vasodilation) -- but this model predicts that
selection does *not* protect them from the chronic AMLO-escape route,
because that route was built to run around the brake entirely. That is a
candidate, model-generated explanation for a known clinical puzzle: why
some well-documented acute CCB "responders" still lose response over years
of chronic treatment. It also reframes section 7's ecological hypothesis:
`ORAI1` subclone status there was treated as *somatic/subclonal*
within-patient heterogeneity; `SNX29` genotype is *germline* --
uniform across all of one patient's SMCs -- so the natural extension is
not a second cell definition but a **second whole-patient parameterization**
(re-running the section-7 ecology study at `$f_snx29` = 0.5 vs 2.0 to ask
whether low-SNX29 patients show a faster or larger-magnitude ORAI1+ sweep),
not yet run.

**Caveats:** the GWAS establishes rs8057488's association with AVD response
and SNX29's gene-level biochemistry (silencing enhances SOCE; overexpression
attenuates constriction) -- it does not, in the abstract-level evidence
available here, establish which *rs8057488 allele* raises vs lowers SNX29
expression/function; the model's three `$f_snx29` levels are illustrative
placeholders from the original curation, not fitted to the GWAS effect
size. `$u_amlo_escape`'s own calibration caveats (section 5) apply
unchanged. Sweep script/cfgs: `results/snx29_sweep/`; data:
`snx29_sweep_results.csv`; figure: `snx29_brake_vs_amlo_bypass.png`.

## 9. A better-grounded alternative: the TRPC6/NFAT feedback loop

Section 1's escape mechanism (Johnson et al. 2020) is scientifically
contested (its own caveat, above). TRPC6 overexpression in PASMC from PAH
patients, and the TRPC6/calcineurin/NFAT transcriptional positive-feedback
loop that drives it, are independently well-replicated and NOT contested
(Kuwahara K et al., *"TRPC6 fulfills a calcineurin signaling circuit during
pathologic cardiac remodeling,"* J Clin Invest 2006, PMC1635163; Yu Y et
al., *"Inhibition of phosphodiesterase-5 suppresses calcineurin/NFAT-
mediated TRPC6 expression in pulmonary artery smooth muscle cells,"*
PMC5519690 -- the latter specifically in PASMC, driven by ET-1 via
calcineurin/NFATc4). Critically, this mechanism requires **no direct
amlodipine-protein interaction at all**: TRPC6 and CaV1.2 are structurally
and pharmacologically distinct channel families, so amlodipine doesn't
touch this loop. The "escape" story becomes conservative rather than
exotic: the drug correctly blocks the one route it targets; the other
route (TRPC6, driven by disease-associated ET-1/hypoxia tone) keeps getting
amplified regardless.

`models/ModeloPASMC_base_v2_TRPC6_feedback.bnd` (built on the **base**
model, not the contested AMLO_escape variant, to isolate this mechanism
cleanly) adds exactly the two edges the base model was missing to realize
this loop -- it already gave `Trpc6`, `Ca_local` (SOCE), and `Ca_ltcc`
(L-type) each their own microdomain ("cada canal el suyo"), but `Trpc6`
had no path back to `Calcineurina`/`Nfat` at all:

1. A new `Ca_trpc6` microdomain, built identically to `Ca_local`/`Ca_ltcc`,
   added as a third, independent input to `Calcineurina`.
2. `Trpc6` gains a second, NFAT-driven transcriptional activation branch
   (`$u_nfat_trpc6`, new and uncalibrated -- set equal to the model's
   existing generic slow-transcription rate `$k_transc`) alongside its
   original fast PLC/DAG-gated branch.

**SNX29 (section 8) cannot protect against this route.** `$f_snx29` only
brakes the `Stim1_pm`/`Soce`/`Ca_local` input to `Calcineurina` -- it has
no connection to `Ca_trpc6`. A patient with a maximally protective SNX29
genotype is fully exposed to this second, SNX29-blind route.

### 9.1 Gradient sweep: AMLO x ET1 x HIPOXIA x K+ (40mM KCl)

72 MaBoSS runs ({base, trpc6_feedback} model x AMLO{0, 0.5, 1.0} (istate
probability, i.e. dose) x ET1{0, 0.5, 1.0} x HIPOXIA{off, on} x K_ALTO{off,
on}), steady state (t=600). Figure: `gradient_sweep_trpc6_vs_base.png`;
full table: `gradient_sweep_results.csv`.

**K_ALTO (the 40mM KCl depolarization test) contributes almost no
independent variation** to the TRPC6 extension's effect -- its main action
(`Em_depol -> Ltcc`) is orthogonal to the `Trpc6`/`Calcineurina` axis, so
the two K_ALTO columns in the figure are nearly identical. One practical
note from building this sweep: the combination HIPOXIA=on AND K_ALTO=on
pushes the Boolean network into a much larger explored state space and
MaBoSS runs ~30x slower under it (~30s vs ~1s per run) -- worth knowing
before scaling this kind of sweep up.

**The TRPC6 extension's effect on `CREC_1` flips sign with hypoxia --
and this needs a caveat, not just a result.** Under normoxia
(`HIPOXIA`=off), adding the TRPC6/NFAT loop increases `CREC_1` by up to
+0.10, largest at moderate ET1 and either AMLO extreme. Under hypoxia
(`HIPOXIA`=on), the SAME extension *decreases* `CREC_1` by up to -0.12 at
high AMLO/high ET1. The most likely explanation is a ceiling effect, not a
sign-reversed causal mechanism: `CREC_1 { logic = Nfat OR Hif1a; rate_up =
(Hif1a ? $u_crec_hif : 0) + (Nfat ? $u_crec_nfat : 0) }` -- under hypoxia,
`Hif1a` alone already drives `CREC_1` to ~0.7, independent of `Nfat`,
leaving little headroom for the new NFAT-mediated contribution to add on
top; the new feedback loop's effect on `Nfat`'s own stochastic dynamics
near that saturating regime can show up as a small net decrease rather
than the increase seen at lower baseline occupancy. **This explanation is
plausible but not independently confirmed** -- flagged explicitly because
this exact model's own version history already documents a near-identical
failure mode (the `.bnd` header's own note: "el problema principal del v7
(resistencia a SNX29 bajo hipoxia) era ESE artefacto" -- v7's main problem
was precisely a hypoxia-dependent artifact of this general shape). Treat
the normoxic-vs-hypoxic sign flip as a real, interesting modeling result
worth further scrutiny (e.g. checking whether it survives a longer
`max_time` or a different `$u_nfat_trpc6` value), not as a settled
biological prediction.

Caveats: `$u_nfat_trpc6` is a new, uncalibrated placeholder (set equal to
`$k_transc` by convention, not fit); this sweep used `sample_count=150`
(reduced from the usual 1000 for the 72-run grid's tractability) and only
3 levels each for AMLO/ET1, so fine dose-response structure between grid
points is not resolved.


## 10. From well-mixed sweep to spatial model: a PhysiBoSS vessel wall with real substrate gradients and closed-loop feedback

The section 9 gradient sweep treats AMLO/ET1/HIPOXIA/K+ as scalar,
cell-independent inputs applied identically to every agent -- there is no
notion of *where* a cell sits relative to the lumen, no diffusion, and no
way for a cell's own intracellular state to feed back onto the substrate
field other agents see. This section replaces that well-mixed
approximation with a spatially-resolved PhysiCell + PhysiBoSS model of a
vessel cross-section, built on `ModeloPASMC_base_v2_TRPC6_feedback` (the
section 9 model) with real diffusing AMLO/ET1/O2/K_ALTO fields and an
explicit closed intracellular -> secretion/uptake -> microenvironment loop.

### 10.1 Geometry

A PAH-thickened vessel cross-section (`generate_physicell_config.py`'s
`build_vessel_cross_section` logic, same convention as the earlier
Control-vs-PAH figures): lumen radius 100 um, media (smooth muscle)
thickness 170 um -> media outer radius 282 um, adventitia outer radius
372 um. Cells placed via the PAH-fitted type-fraction model:
`results/spatial_gradient/cells.csv` holds 1,186 agents -- Smooth_Muscle_Cell
n=722, Fibroblast n=367, Endothelial_Cell n=52 (lining the lumen),
STEAP4hi_SMC_Like n=45 (a PAH-expanded SMC-like population identified in
the snRNA-seq data, wired identically to SMC). Endothelial cells and
fibroblasts are "bystander" cell types in this model: Endothelial_Cell
carries the substrate secretion sources but no Boolean network;
Fibroblast carries neither.

### 10.2 Input mapping: substrate -> Boolean node

PhysiBoSS's native `<input>` mapping (`addons/PhysiBoSS/src/utils.h`)
drives node values directly from local substrate concentration, no custom
C++ required: `action="activation"` turns a node ON above `threshold`;
`action="inhibition"` turns it ON *below* `threshold` (and back off above
`inact_threshold`) -- exactly the logic O2 -> HIPOXIA needs. Four inputs
are wired on every SMC/STEAP4hi_SMC_Like agent: `AMLO` substrate ->
`AMLO` node (activation), `ET1` substrate -> `ET1` node (activation), `O2`
substrate -> `HIPOXIA` node (inhibition), `K_ALTO` substrate -> `K_ALTO`
node (activation).

### 10.3 Output mapping: Boolean node -> secretion/uptake (the closed loop)

PhysiCell's behavior vocabulary natively exposes `"<substrate> secretion"`
and `"<substrate> uptake"` as settable behaviors
(`core/PhysiCell_signal_behavior.cpp`), which PhysiBoSS's `<output>`
mapping can target directly -- so the intracellular network can modify
the microenvironment it reads from, with no custom C++. Two feedback
edges are wired (only in the `feedback=on` arm):

- **`Nfat` -> `"ET1 secretion"`**: an SMC/STEAP4hi cell with `Nfat`=1
  increases its own ET1 secretion rate, closing an autocrine loop.
  Grounding: PASMC are a well-documented autocrine/paracrine source of
  ET-1 in their own right (Fagan et al. 1999; reviewed in MacLean &
  Dempsie 2010), and calcineurin/NFAT-ET-1 coupling in smooth muscle is
  independently described (Yu Y et al., PMC5519690, NFATc4-driven in
  PASMC). `Nfat` is used here as a general stress/activation proxy for
  this edge, **not** a literal claim that NFAT is ET1's transcription
  factor -- the mechanistic specificity of the NFAT->ET1 link itself is
  weaker than the Kuwahara/Yu TRPC6/calcineurin/NFAT loop it sits next to
  (see section 9's caveats), and this is flagged as a modeling choice, not
  a citation-backed edge.
- **`Ca_global` -> `"O2 uptake"`**: a cell with elevated global Ca2+
  increases its own O2 consumption. Grounding: activated-cell
  hypermetabolism raising local O2 demand is an uncontested, general
  cell-biology principle (not PAH-specific), used here mainly to close
  the loop in the direction PAH pathophysiology actually goes: more Ca2+
  signaling -> more local hypoxia -> more HIF1a -> more CREC_1 -> (via
  the section 9 wiring) more Ca2+ signaling again.

### 10.4 Two crashes found and fixed

Building this model surfaced two PhysiBoSS/PhysiCell bugs independent of
the modeling choices above:

1. **`$time_scale` null pointer** (same bug as the section 7 ecology
   model): the sample project's `pre_update_intracellular()` in
   `custom.cpp` references an undefined `$time_scale` MaBoSS parameter.
   Fixed by removing the reference (the `.cfg` already sets
   `intracellular_dt`/`scaling` explicitly).
2. **New: `my_coloring_function()` / `color_node()` segfault on bystander
   cell types.** These helper functions in `custom.cpp` called
   `pCell->phenotype.intracellular->...` unconditionally on every agent
   for visualization coloring; `Fibroblast` and `Endothelial_Cell` have no
   `<intracellular>` block, so `phenotype.intracellular` is `NULL` for
   them, and the uniform-population sample code did not anticipate a
   mixed Boolean/non-Boolean population. Fixed with a null check before
   dereferencing.

Both `feedback=off` and `feedback=on` project variants build clean and run
to completion (exit 0, 4320 simulated minutes, ~6s wall time each; cell
division minimized -- `division_rate` ~1e-5 to 1e-7 -- to isolate the
spatial-substrate response from population growth as a first test, not to
model vessel-wall proliferation).

### 10.5 Calibration note: thresholds set to the achieved range, not literal physiological units

Two parameter issues were found and fixed by inspecting the final
microenvironment state before trusting any Boolean readout from it:

- **Visibility bug (not a dynamics bug):** `AMLO`, `ET1`, `HIPOXIA`, and
  `K_ALTO` are marked `is_internal=TRUE` in the base `.cfg` (so they don't
  clutter the MaBoSS state-trajectory output in well-mixed runs). This
  also hides them from PhysiBoSS's per-agent `*_boolean_intracellular.csv`
  export, making every one of these nodes read as permanently `False`
  regardless of the true simulated state. Fixed by flipping
  `is_internal=FALSE` for these four nodes in the spatial-model `.cfg`
  copy -- a pure export-visibility change, not a logic change.
- **Threshold-vs-achieved-concentration mismatch:** the O2 field's steady
  state (set by SMC uptake, EC secretion, and fast diffusion / slow decay
  over a domain much smaller than the diffusion length constant) settles
  to ~0.10-0.68 (arbitrary density units once rescaled by the secretion
  target), not the literal ~38 mmHg physiological initial condition a
  naive threshold choice assumed; the same was true for the AMLO
  activation threshold (default 0.3, achieved range only 0.035-0.037) and
  the ET1 threshold. All three thresholds (`--o2-threshold`,
  `--o2-inact-threshold`, `--amlo-threshold`, `--et1-threshold`) were
  recalibrated to sit inside the empirically achieved concentration range
  rather than at an absolute physiological value, so that the Boolean
  nodes can actually flip both ways across the vessel wall. **This is a
  methodological convention, not a biophysical calibration** -- the model
  demonstrates the right qualitative structure (steeper drug exposure near
  the lumen, diffusion-limited hypoxia deeper in the media, autocrine
  amplification under feedback) but the absolute substrate-concentration
  scale should not be read as quantitatively calibrated to real tissue
  pharmacokinetics.

### 10.6 Results: spatial structure the well-mixed sweep cannot show

Figure: `spatial_vessel_gradient_feedback.png`; per-agent data table:
`spatial_gradient_per_agent_profiles.csv` (1,205 SMC/STEAP4hi agents
pooled across both arms, with per-agent radial position, local substrate
concentration, and Boolean node state at t=4320 min).

**Feedback OFF** (substrate gradients only, no intracellular->microenvironment
coupling): AMLO exposure is highest near the lumen and falls off toward
the adventitia (Pearson r with radius = -0.76 for the AMLO node's driving
substrate); `HIPOXIA` shows the opposite, expected gradient -- essentially
absent near the endothelium-lined lumen and reaching ~100% of SMCs by
~150-180 um into the media (diffusion-limited medial hypoxia, a
recognized PAH-relevant pathophysiological concept -- the thickened media
outstrips O2 diffusion from the lumen). The convergent readout `CREC_1`
(`Nfat OR Hif1a`) tracks this: 17% ON near the lumen, rising to ~65-70%
deeper in the media, driven mostly by the hypoxia arm of the OR once
`HIPOXIA` saturates.

**Feedback ON** (`Nfat`->`ET1 secretion`, `Ca_global`->`O2 uptake`): the
autocrine ET1 loop does not stay spatially localized -- because ET1's
diffusion length constant is large relative to the ~180 um media
thickness, the self-reinforcing secretion saturates essentially the whole
medial cross-section (ET1 substrate concentration rises ~5x over the
no-feedback baseline, uniformly). Downstream consequences across the
whole SMC population: `Soce` ON-fraction rises from a radially-graded
16-24% (off) to a uniform ~43-48% (on); `Trpc6` rises even more sharply,
from a near-silent 2-9% (off) to a uniform 64-77% (on); `CREC_1` rises
from 58% (off, population mean) to 71% (on). `HIPOXIA` also increases
modestly (86% -> 97% of SMCs) via the `Ca_global`->O2-uptake edge
amplifying local O2 demand on top of the diffusion-limited baseline.

**The qualitative finding the well-mixed sweep could not produce:** a
closed intracellular-to-microenvironment loop does not merely rescale the
well-mixed steady state -- it can convert a radially *graded*,
position-dependent signal (drug exposure, hypoxia) into a spatially
*uniform*, self-sustaining pathway activation once the autocrine
substrate's diffusion range exceeds the tissue compartment size. This is
precisely the kind of emergent, geometry-dependent behavior a
single-compartment Boolean model is structurally unable to represent, and
is the model-level analogue of the clonal-sweep result from section 7
(local perturbations propagating to take over the whole tissue), but
acting through a continuous diffusing field rather than cell division and
competition.

### 10.7 Caveats

- Division/mechanics are minimized by design (see 10.4); this is a
  steady-state substrate-response snapshot, not a model of PAH medial
  thickening or SMC proliferation dynamics.
- `K_ALTO` is modeled as a uniform bath initial condition with no
  secretion source (the 40mM KCl depolarization test is experimentally
  bath-applied, not a physiological spatial gradient), so it is
  deliberately *not* part of the spatial-gradient story here.
- The feedback rate parameters (`Nfat`->ET1 secretion rate, `Ca_global`->O2
  uptake rate) are uncalibrated placeholders, swept for qualitative
  behavior, not fit to any measured autocrine ET1 secretion rate or O2
  consumption rate.
- As in section 9, this model inherits the base model's
  Calcineurina/Nfat/Trpc6 wiring and its own documented hypoxia-interaction
  caveat (9.1) -- the same caution about ceiling effects under hypoxia
  applies here.
- Two independent 4320-minute runs (one per arm, single seed) -- no
  replicate seeds were run, so the reported fractions are point estimates,
  not confidence intervals.


### 10.8 Correction: the original two-arm design never tested amlodipine addition at all

The two arms in sections 10.1-10.7 (`feedback=off` / `feedback=on`) used
an **identical, constant** endothelial AMLO secretion rate in both --
there was no "no drug" vs "chronic AMLO" comparison anywhere in the
spatial model, even though that comparison is the entire premise of this
project. This was caught by inspection, not by construction, and is
corrected here with two additional arms (`amlo_secretion_rate=0`, i.e. a
vehicle/no-drug baseline) crossed with both feedback states, giving a
proper 2x2: {no_drug, chronic_amlo} x {no_fb, fb_on}.

**Per-cell, the drug blockade works exactly as intended:** within the
`chronic_amlo` arms, cells whose local AMLO concentration crosses
threshold (`AMLO`=1) show `Ltcc` ON in only 1.9-5.4% of boolean-trajectory
time, versus 83.9-85.7% for cells below threshold in the same run -- a
>15-fold per-cell drop, confirming `AMLO -| Ltcc` fires correctly.

**But only ~21% of the medial SMC population ever reaches that
threshold.** AMLO's diffusion length constant (sqrt(4000/0.002) = 1414 um)
is much larger than the ~180 um media, so -- despite continuous EC
secretion -- the field settles nearly flat and crosses the activation
threshold only in the innermost ~1-2 radial bins (r < ~160 um). Radially:
`Ltcc` ON-fraction is ~0% at r=115 um and ~3-9% at r=145 um in the
`chronic_amlo` arms, but **jumps back to the undrugged baseline
(~70-87%) by r=175 um and beyond** -- a sharp therapeutic trough confined
to roughly the inner third of the media, not a vessel-wide blockade. The
`no_drug` arms show no such trough (`Ltcc` stays at baseline, ~75-91%,
across the full radius).

**Vessel-wide population averages are therefore diluted by a large,
effectively undosed outer-medial subpopulation:** `Ltcc` ON-fraction
across the whole SMC population drops only from 75.4% (no drug) to 68.4%
(chronic AMLO, no feedback) -- a modest net effect that would look like
"partial drug resistance" if read at the tissue-average level, even
though the per-cell pharmacology is working perfectly on every cell it
actually reaches.

This is a **third, independent candidate explanation for attenuated
amlodipine efficacy in PAH**, orthogonal to both the contested
STIM1/escape mechanism (section 1) and the SNX29 genetic-brake story
(section 8): simple diffusion-limited drug penetration into a
PAH-thickened medial wall, a standard and uncontested
pharmacokinetic/pharmacodynamic principle, not a cell-autonomous
adaptation mechanism at all. It also reframes the original
section 9-12 feedback story: the `feedback=on` autocrine amplification
(section 10.6) still saturates the *whole* vessel wall with ET1/TRPC6
activation regardless of how little of it the drug itself reaches --
i.e. a feedback loop downstream of the drug's target can propagate tissue
-wide even when the drug's own spatial reach is narrow.

Figure: `spatial_amlo_addition_effect.png`. Full 2x2 per-agent dataset:
`spatial_gradient_per_agent_profiles_full2x2.csv`. Project bundles:
`vessel_gradient_novehicle_off_bundle.tar.gz`,
`vessel_gradient_novehicle_on_bundle.tar.gz` (vehicle arms; chronic-AMLO
arms are the section 10.1-10.7 bundles).

**Caveat on this correction itself:** the specific penetration depth
(~160 um) is set by the same uncalibrated `amlo_secretion_rate`/threshold
choice flagged in 10.5, not a fitted pharmacokinetic parameter -- the
qualitative "trough, not uniform blockade" structure is the robust
finding; the exact radius at which it crosses over is not.


### 10.9 Vessel growth, corrected: neointima formation, not medial thickening

A user question ("we see no growth in the vessel, no?") exposed that
division was deliberately minimized in 10.1-10.8 (`division_rate` ~1e-5 to
1e-7/min) specifically to isolate the substrate/feedback result from
population dynamics -- confirmed, growth was real but negligible (1,186 ->
1,205 agents, ~1.6%, over 4,320 min).

**First attempt (medial growth) failed structurally, not numerically.**
Re-running with a realistic ~2.5-day doubling time for `CREC_1`=1 cells
over 20 simulated days grew the population 1,186 -> ~31,000-35,000, but
the media's radial footprint barely expanded (SMC r-range 117-277 ->
109-291 um) -- cells were overlapping (up to 96% below touching distance
by day 7), not pushing the wall outward. Root cause: the PAH media was
deliberately built at `packing_fraction=0.85` (section 4's medial-
hyperplasia design choice) -- already near its mechanical packing limit,
leaving almost no room for additional growth regardless of pressure-gated
contact inhibition (tested down to `pressure_half_max=0.03`,
`pressure_hill_power=6`, with little change in outcome).

**The literature redirected the question, and the redirect fixed the
mechanics problem too.** PAH vascular remodeling has two distinct lesion
types, and the user correctly recalled this: Stacher et al. (Circulation
2021, PMC6656344 lineage) and the broader PAH-pathology literature
describe pulmonary arteries narrowing because of medial thickening AND
occlusion by neointimal lesions, with neointima formation often the
dominant contributor to luminal narrowing -- a small number of SMCs are
selected to establish the neointima and then undergo clonal expansion
to form occlusive lesions, frequently traced to a Notch3-marked founder
subpopulation. This maps directly onto this model's existing
`STEAP4hi_SMC_Like` minority population (6% of PAH SMCs, already in the
geometry from section 4) -- and critically, redirecting growth into the
EMPTY lumen sidesteps the packing-limit problem entirely instead of
fighting it.

**A second, independent bug was found and fixed along the way: mechanics
and motility were silently disabled project-wide.** The
`physiboss_cell_lines` sample template this project is built on sets
`cell_defaults.functions.update_velocity = NULL` (confirmed by diffing
against a motile PhysiCell sample, `sample_projects/biorobots`, which sets
it to `standard_update_cell_velocity`). With this NULL, neither
adhesion/repulsion relaxation after division NOR the motility vector from
any `update_migration_bias` function is ever applied -- this was the
actual reason pressure-gating in the first growth attempt barely helped
(mechanics was never running to relax crowded cells), and also why an
initial founder-cell migration implementation produced zero net movement
over 7 days despite correct Boolean gating (confirmed via targeted
debug instrumentation: `update_migration_bias` was never invoked). Fixed
by setting `update_velocity = standard_update_cell_velocity` in the
generator.

**Implementation:** `STEAP4hi_SMC_Like` gains, only when `CREC_1`=1: (1) a
realistic division rate (~2.5-day doubling, same value tested in the
failed medial-growth attempt) and (2) directed migration ("migration
speed" native output, gated by `CREC_1`) toward the lumen center via a
new custom `update_migration_bias` function (no native PhysiCell XML
primitive exists for "direction toward a point"). Regular
`Smooth_Muscle_Cell` keeps `division_rate=0` (no medial growth, avoiding
the packing-limit problem entirely).

**Result:** both arms (no-drug and chronic-AMLO, `feedback=off`) show
clean, fast neointimal formation -- founder cells (45-46 of them) leave
the media essentially completely within 5 simulated days (46/46 and
46/47 inside the lumen by day 5, final cluster radius ~50-58 um, well
inside the 100 um lumen), visualized as a compact occlusive mass forming
inside the endothelial ring. **The two arms are essentially
indistinguishable** -- AMLO exposure does not measurably alter the
timing or extent of neointima formation at this calibration, most likely
because founder cells start outside the ~160 um AMLO penetration depth
(section 10.8) and, once packed inside the avascular lumen, `CREC_1`
becomes dominated by hypoxia (`Hif1a`) rather than the AMLO-sensitive
`Nfat` arm -- a specific, falsifiable prediction of this particular
wiring, not a general claim that AMLO cannot affect neointima formation.

Figure/GIF: `neointima_formation_timelapse.gif`. Project bundles:
`vessel_gradient_neointima_test2_bundle.tar.gz` (chronic AMLO),
`vessel_gradient_neointima_novehicle_bundle.tar.gz` (no drug).

**Caveats:**
- Founder-cell overlap in the lumen is moderate, not resolved (~94% of
  pairs closer than touching distance at day 7) -- the sustained inward
  migration bias continuously opposes repulsion-driven relaxation; this
  is a known, generally-accepted pattern for directed-migration models
  (comparable to chemotactic spheroid compaction) but is not a
  quantitatively validated neointimal density.
- `STEAP4hi_SMC_Like` as the literature's Notch3+ founder-cell analogue is
  an inferred structural parallel (both are PAH-expanded minority SMC
  subpopulations), not a validated marker correspondence -- this
  project's snRNA-seq/Xenium data was not re-checked for Notch3
  co-expression with STEAP4 before making this assignment.
- Migration speed (0.3 um/min) and bias (0.8) are new, uncalibrated
  placeholders, chosen to produce a legible timescale (~5 days) rather
  than fit to any measured SMC migration rate.
- The update_velocity fix changes mechanics behavior for ALL cells in
  every project generated by this script going forward, including the
  already-delivered section 10.1-10.8 results -- those results used
  near-zero division rates (minimal position change expected either way)
  and are not invalidated, but were not re-run under the corrected
  mechanics; this is disclosed rather than silently left inconsistent.


### 10.10 Correction: mechanical crowding alone doesn't invade the lumen -- replaced the point-attractor with ET1 chemotaxis

A second user critique ("I'd expect these cells to push against one
another, not move toward the center") correctly identified that 10.9's
`neointimal_migration_bias` function was scripted: it unconditionally set
`migration_bias_direction = -position/|position|`, i.e. "always point at
the mathematical origin" -- not a signal any real cell senses, and a much
weaker mechanistic claim than it looked.

**Tested the alternative directly: pure division + mechanics, no directed
bias at all.** Over 7 simulated days, zero founder cells crossed into the
lumen (minimum founder radius stayed frozen at ~120-122 um, the outer
radius instead crept outward). This is itself an informative, reportable
negative result: passive mechanical crowding cannot breach the
endothelial monolayer sitting at the lumen boundary. Biologically this is
consistent with how neointimal invasion is actually described -- crossing
the endothelium/internal elastic lamina is an active, directed process,
not something that happens for free from being pushed by neighbors.

**Fix: replaced the scripted point-attractor with PhysiCell's native
`<chemotaxis>` mechanism**, directing `STEAP4hi_SMC_Like` founder cells up
the existing ET1 gradient (toward its endothelial source at the lumen) --
a real, already-present substrate field, not a new fictional point
attractor, and endothelin-driven VSMC migration is itself independently
documented pharmacology (distinct from the contested AMLO/STIM1 story).
This required no custom C++ at all: PhysiCell's own XML-driven
`chemotaxis_function` replaces the custom `update_migration_bias` pointer
automatically when `<chemotaxis><enabled>true</enabled>` is set -- the
custom `neointimal_migration_bias` function and its `find_cell_definition`
wiring were removed from the generator entirely.

**Result: a more realistic, partial (not 100%-deterministic) invasion.**
35-44 of ~55 founder cells reach the lumen by day 7 in both arms (vs. a
complete, uniform sweep under the old point-attractor), forming an
eccentric, crescent-shaped occlusive mass on one side of the lumen with
visible local endothelial displacement -- closer to the asymmetric
presentation of real neointimal lesions than a circumferentially uniform
ring. The two arms (no-drug, chronic-AMLO) remain essentially
indistinguishable, for the same reason given in 10.9 (founder cells
start outside the AMLO penetration depth; `CREC_1` becomes
hypoxia-dominated once cells pack into the avascular lumen).

Figure/GIF: `neointima_chemotaxis_timelapse.gif` (replaces
`neointima_formation_timelapse.gif`, which used the retracted
point-attractor mechanism and is kept only for the debugging record).
Project bundles: `vessel_gradient_chemotaxis_amlo_bundle.tar.gz`,
`vessel_gradient_chemotaxis_novehicle_bundle.tar.gz`.

**Caveats:** migration speed (0.3 um/min) and `migration_bias` (0.8,
directedness) remain new, uncalibrated placeholders; chemotaxis strength
is not fit to any measured VSMC migration rate. The partial-invasion
fraction (35-44/55) is a function of these same placeholders and the
7-day window chosen for visualization, not an independently predicted
clinical neointimal-occlusion fraction.


### 10.11 Correction: PhysiCell cell volume was inconsistent with the established geometry's realistic diameters

A third correction, caught by inspection: every `smc_phenotype_block`/
`bystander_phenotype_block` cell in sections 10.1-10.10 used PhysiCell's
sample-template default `<volume><total>1000</total></volume>`
(radius ~6.2 um, diameter ~12.4 um) -- but `build_vessel_cross_section.py`
(sections 2/4) was already built around realistic, explicitly-chosen
diameters: `--smc-diameter 15`, `--fibroblast-diameter 15`,
`--ec-diameter 12` (um), used to space the initial, non-overlapping cell
lattice. The spatial PhysiBoSS model's simulated cell size was therefore
~17% smaller in diameter than the geometry it was placed into assumed --
not just a cosmetic mismatch: **PhysiCell's own adhesion/repulsion
mechanics uses the simulated cell radius to decide when cells push
apart**, so the undersized volume meant real mechanical contact was
detected later (cells could get closer) than the geometry's own spacing
implied.

**Fixed** by setting `<volume><total>` to match the established diameters:
1767.15 um^3 for SMC/STEAP4hi/Fibroblast (15 um diameter,
`(4/3)*pi*7.5^3`), 904.78 um^3 for Endothelial_Cell (12 um diameter) --
`bystander_phenotype_block` was parameterized with a `total_volume`
argument (previously hardcoded) so EC and Fibroblast can each get their
own correct value. Re-running the chemotaxis-based neointima arms
(section 10.10) with the corrected volumes showed, if anything, MORE
founder cells reaching the lumen (45/55 and 44/53, vs 35-44/55 before) --
consistent with properly-calibrated repulsion more effectively relieving
local crowding and pushing cells toward the open lumen space.

The figure/GIF was also re-rendered using true-to-scale cell circles
(matplotlib `EllipseCollection` with diameters in data units, not
arbitrary scatter-marker point sizes) so the vessel wall's packing
density and the neointimal mass's true footprint are visually accurate,
not just schematically placed dots.

Figure/GIF: `neointima_chemotaxis_realistic_radius.gif` (supersedes
`neointima_chemotaxis_timelapse.gif` for the cell-size claim; the
substrate-gradient/feedback figures in 10.1-10.8 and the drug-penetration
figure in 10.8 were not re-rendered, since their SMC positions are static
in those runs -- the volume correction only matters where mechanics
actively moves cells, i.e. sections 10.9-10.11's growth/migration
experiments). Project bundles:
`vessel_gradient_chemotaxis_amlo_r_bundle.tar.gz`,
`vessel_gradient_chemotaxis_novehicle_r_bundle.tar.gz`.


### 10.12 Tightening the hypothesis: an impenetrable endothelium, and isolating LCCB/SNX29/NFAT from the hypoxia confound

Two corrections, prompted by user feedback on the biology, replace the
section 10.10 result with a cleaner and more mechanistically grounded one.

**1. Soft mechanics are not enough for "cannot traverse."** Testing
showed that even an 8x boost to EC-EC adhesion/repulsion (section 10.10's
first attempt) was insufficient: the EC ring itself got dragged apart
under sustained founder-cell pressure (EC radial spread widening from
~104-108 um to 25-1278 um), and founder cells still ended up deep in the
lumen. A direct test confirmed the underlying claim this session's
earlier discussion made explicit: pure mechanics with NO directed bias
produces ZERO founder-cell lumen entry over 7 days (min radius frozen at
~120 um) -- crossing an intact endothelium is not something passive
crowding does on its own, matching how neointimal invasion is described
in the literature (an active process).

**2. A hard positional clamp caused a severe, confirmed PhysiCell bug.**
The first fix -- calling `Cell::assign_position()` every mechanics step
to pin founder-cell radius >= a boundary value -- caused a catastrophic
numerical blowup (agent radii reaching 10^4-10^5 um within 1-2 simulated
days) that persisted across every variant tested: moving the clamp
radius outside the EC ring, reducing `dt_mechanics` 10x, damping
velocity after the teleport, and even removing cell division entirely
(fixed 45-cell population). Tracing `assign_position()`'s source
(`core/PhysiCell_cell.cpp`) found the root cause: it calls
`Cell_Container::register_agent()`, which only appends the cell to its
new voxel's neighbor list (`agent_grid[voxel].push_back(agent)`) --
**it never removes the cell from its previous voxel.** Called once per
mechanics step for a cell sitting at a boundary (as the clamp did,
continuously), this duplicates the cell's own neighbor-list registration
without bound, so repulsion forces from that single (growing, duplicated)
entry get counted an ever-increasing number of times -- a textbook
unbounded-feedback numerical blowup, confirmed by the explosion's onset
at day 1-2 and its reduced severity (but not absence) when the clamp fired
less often. The correct, safe way to move an existing cell is PhysiCell's
own `update_position()` -> `update_voxel_in_container()` pathway, which
does remove-then-add correctly -- but that path depends on internal state
(`updated_current_mechanics_voxel_index`) set during the normal velocity
integration, and several of the functions involved (`current_mechanics_
voxel_index`, the `Cell_Container*` member) are private to `Cell`, so a
free function cannot safely replicate it directly either.

**Fixed with a velocity-based barrier instead of a position teleport**:
`founder_barrier_velocity()` calls the standard mechanics update first
(`standard_update_cell_velocity`), then cancels any inward-pointing
component of the resulting velocity when a founder cell is at/inside the
boundary radius, plus a small outward restoring term proportional to
penetration depth (a soft contact spring, not a hard stop). Because this
only ever modifies `velocity` -- never `position` or voxel registration
directly -- all position changes still go through PhysiCell's own correct
`update_position()` pathway. Stable in every test since (fast runtimes
matching the pre-bug baseline, ~35-40s per 7-day run, vs 90-140s and
eventual numerical garbage with the position-teleport version). Result:
founder cells settle into a tight, stable band just outside the
endothelium (r = 115-132 um vs the EC ring's own 104-108 um -> 78-109 um,
visibly compressed inward by the accumulating subendothelial mass) --
genuine neointimal narrowing of the lumen without any cell crossing an
intact endothelial monolayer. See
`subendothelial_neointima_snx29_comparison.png`.

**3. SNX29 was wired into founder cells, and its effect was initially
swamped by the same hypoxia ceiling effect documented in section 9/10.9.**
A 6-run sweep ({SNX29 low/normal/high} x {no drug, chronic AMLO}, via the
existing `$f_snx29` parameter, now exposed as `--f-snx29` in the
generator) first used the same `CREC_1` (`Nfat OR Hif1a`) gate as
section 10.9-10.11's founder division/migration. Result: `Hif1a` sat at
68-81% in ALL SIX conditions regardless of SNX29 or drug, confirming
`CREC_1` was ceiling-dominated by the self-generated hypoxia inside the
packed subendothelial mass -- exactly the confound flagged when the user
asked "why do no-drug and AMLO look the same," now recurring in a
different experiment and masking any SNX29/NFAT signal.

**Isolated the LCCB/SNX29/NFAT axis by re-gating founder cells'
"cycle entry"/"migration speed" outputs on `Nfat` alone** (new
`--founder-pathway-node` argument, default `CREC_1` for backward
compatibility with regular SMC, set to `Nfat` for this experiment).
Re-running the same 6-condition sweep with this change produced a clean,
monotonic, mechanistically consistent result:

| SNX29 | Soce (no drug) | Soce (chronic AMLO) | Nfat (no drug) | Nfat (chronic AMLO) |
|---|---|---|---|---|
| low    | 0.250 | 0.345 | 0.882 | 0.800 |
| normal | 0.094 | 0.156 | 0.734 | 0.750 |
| high   | 0.088 | 0.046 | 0.719 | 0.692 |

**SNX29 main effect (averaged over drug): Soce 0.298 (low) -> 0.125
(normal) -> 0.067 (high)** -- a clean, monotonic brake effect, exactly
matching the Karnes 2026 GWAS direction (low SNX29 = weak brake = more
SOCE; high SNX29 = strong brake = less SOCE), now visible because the
hypoxia ceiling no longer masks it.

**AMLO's effect on Soce depends on SNX29 genotype**: chronic AMLO
*increases* Soce at low SNX29 (+0.095) and normal SNX29 (+0.063) --
consistent with the section 1-disclosed, contested Johnson et al. 2020
STIM1/escape mechanism -- but *decreases* it at high SNX29 (-0.042),
the opposite direction. A plausible reading: a strong SNX29 brake
suppresses the escape route enough that AMLO's intended pharmacology
(blocking Ltcc, reducing calcium influx into the calcineurin/NFAT/SOCE
loop) dominates instead -- i.e., SNX29 genotype could determine whether a
patient is escape-prone (low/normal SNX29, drug-exposure paradoxically
sustains the pathway) or escape-resistant (high SNX29, drug works as
intended). This is the kind of mechanistic link between a real GWAS hit
and chronic-therapy response the project set out to explore -- but it
rests on a SINGLE seed per condition (N=55-68 founder cells each), so the
SNX29 main effect (clean, monotonic, large) should be trusted far more
than the AMLO x SNX29 interaction (smaller delta, direction-sensitive,
not replicated across seeds).

Figures: `snx29_amlo_nfat_pathway.png` (Soce/Nfat vs drug, by SNX29
level), `subendothelial_neointima_snx29_comparison.png` (spatial,
true-to-scale). Data: `snx29_amlo_nfat_sweep_results.csv`. Project
bundles: `nfat_snx29_sweep_bundles.tar.gz` (all 6 runs).

**Caveats:** `--f-snx29` levels (0.5/1.0/2.0) and `--founder-pathway-node`
are both pre-existing/newly-exposed model mechanisms, not new rate
constants, but the founder-specific migration speed/bias and the
velocity-barrier's restoring-force coefficient (0.5, hardcoded) remain
uncalibrated placeholders. Single-seed sweep -- the SNX29 main effect's
monotonicity across 3 independent levels is reassuring but not a
substitute for replicate seeds, which were not run given the compute
budget already spent debugging the barrier mechanics.


### 10.13 Differentiated apoptosis: founder cells are apoptosis-resistant, regular SMC are not

A user observation pointed out a missing mechanistic distinction: regular
SMC and the STEAP4hi_SMC_Like founder subpopulation both had apoptosis
rate = 0 (death model 100, hardcoded) -- i.e. no cell death of either type
was represented at all. This is biologically wrong in a specific,
well-documented way: the literature on PAH neointima formation describes
acquired apoptosis resistance in the proliferative, "cancer-like"
PASMC phenotype as a *necessary* condition for neointima formation, not an
incidental one, and distinguishes it from the more modest apoptosis
reduction seen in bulk PAH PASMCs relative to healthy controls.

**Fix**: added `--smc-apoptosis-rate` (default 5.31667e-05 /min, PhysiCell's
standard template default, ~9-day half-life -- representing partial, not
complete, apoptosis resistance in bulk PAH SMC) and `--founder-apoptosis-rate`
(default 0.0 -- fully apoptosis-resistant) as independent parameters on the
`smc_phenotype_block()` death model, wired separately into the
`Smooth_Muscle_Cell` and `STEAP4hi_SMC_Like` `cell_definition` blocks.
Verified in the generated XML (`Smooth_Muscle_Cell` apoptosis death_rate =
5.31667e-05, `STEAP4hi_SMC_Like` = 0) and confirmed the build+run pipeline
still completes cleanly (26-28s per 7-day run, matching the pre-change
baseline).

**Re-ran the full 6-condition SNX29 x AMLO sweep with this fix.** Two
things changed, one expected and one that forces a retraction of a specific
prior claim:

- **Expected**: the regular SMC population in the media now visibly
  declines over the 7-day run (722 -> 433-472 across arms) due to the
  new nonzero turnover, while the founder population stays flat (45-47) as
  before -- both consistent with the intended design (bulk media cells turn
  over; founder cells accumulate).
- **SNX29 main effect on SOCE is preserved**: low SNX29 (0.289-0.319) remains
  clearly, robustly elevated above normal (0.106-0.130) and high
  (0.109-0.152) SNX29 in both drug conditions -- the core, most
  important finding is NOT an artifact of the missing-apoptosis bug.
- **The AMLO x SNX29 interaction reported in section 10.12 does NOT
  replicate** under this correction. Previously: chronic AMLO increased
  Soce at low/normal SNX29 but decreased it at high SNX29 (a clean sign
  flip). Now: chronic AMLO increases Soce at low (+0.030) AND high
  (+0.043) SNX29, and decreases it slightly at normal (-0.024) -- a
  different, non-monotonic pattern with no obvious mechanistic story. On
  the Nfat readout, the new pattern is in fact simpler and more internally
  consistent than before: chronic AMLO decreases Nfat fraction at ALL
  THREE SNX29 levels (low 0.711->0.681, normal 0.848->0.766, high
  0.870->0.804) -- no sign flip anywhere.

**Conclusion**: the section 10.12 AMLO x SNX29 interaction was correctly
flagged there as preliminary/single-seed, and this is exactly the kind of
non-replication that caveat existed for. It is retracted as a presented
finding. What survives, now on firmer ground (validated against a
biologically-motivated model correction, not just a second identical-model
seed), is the SNX29 main effect on SOCE: low SNX29 genotype produces
meaningfully higher SOCE than normal/high SNX29 in founder cells, in both
drug conditions, matching the Karnes 2026 GWAS direction. The honest
reading of the AMLO interaction is now "not yet established" rather than
"present but sign-dependent on SNX29" -- confirming a true drug-genotype
interaction would need replicate seeds at minimum, and likely a better
mechanistic hypothesis for *why* AMLO's effect on SOCE should depend on
SNX29, before it is worth presenting again.

Figure: `snx29_amlo_nfat_apoptosis_corrected.png`. Data:
`snx29_amlo_nfat_apoptosis_corrected_results.csv`. Bundle:
`apop_snx29_sweep_bundles.tar.gz` (all 6 re-run arms).


### 10.14 ModeloPASMC_v8: Boolean-driven phenotype transformation replaces fixed-fraction "clones"

The user provided a substantially extended Boolean network, `ModeloPASMC_v8`,
adding a **Myocd/Klf4 double-negative toggle** -- the textbook contractile
<-> synthetic VSMC phenotype switch (Myocd maintains the contractile program
and represses Klf4; Klf4 represses the contractile program; Nfat pushes
toward the Klf4/synthetic side) -- plus a **`Resist_apop`** node (`Nfat OR
Klf4`) explicitly documented in the network file's own comments as intended
to be mapped onto the PhysiCell agent's apoptosis rate "once the tissue layer
is built" (i.e., exactly the mechanism section 10.13 had just hand-coded as
two CLI constants). CACNA1C and Kv1_5 are no longer frozen inputs: both now
depend on Myocd with slow (days-scale) kinetics plus a nonzero residual rate,
so a synthetic-state cell loses most, not all, of its CaV1.2/Kv1.5
expression -- distinguishing *channel availability* (slow, transcriptional)
from *channel gating* (fast, voltage/drug-dependent, still handled entirely
by the existing Ltcc/Em_depol nodes).

**This directly enables what the user asked for**: instead of seeding a
fixed ~6% "founder" subpopulation as permanently separate agents (the
section 10.9-10.13 "STEAP4hi_SMC_Like clone" approach), cells now
TRANSFORM between `Smooth_Muscle_Cell` and the renamed `Synthetic_SMC`
(formerly STEAP4hi_SMC_Like -- simplified per the user's request, since the
STEAP4hi marker-gene framing no longer matches the mechanism) at runtime,
driven by their own current Klf4/Myocd Boolean state.

**Implementation, verified to need zero new custom C++** (directly
addressing the "is the istate-override hook too much?" question from
earlier in this session): PhysiCell's `<cell_transformations>` block plus
the SAME `<intracellular><mapping><output>` XML already used for cycle
entry/migration speed. Two native PhysiCell behaviors were confirmed by
direct source inspection (`core/PhysiCell_signal_behavior.cpp`) to be
pre-registered exactly like "cycle entry": `"transform to " + <cell type
name>` and `"apoptosis"`. Added:
- `<output physicell_name="transform to Synthetic_SMC" intracellular_name="Klf4">`
  on `Smooth_Muscle_Cell` (new `--transformation-rate`, default 1.6e-4/min,
  ~3-day half-life once Klf4=1 -- UNCALIBRATED).
- the reverse `<output physicell_name="transform to Smooth_Muscle_Cell" intracellular_name="Myocd">`
  on `Synthetic_SMC` (new `--reverse-transformation-rate`, default 1.6e-5/min,
  10x slower, reflecting the toggle's designed hysteresis -- UNCALIBRATED).
- `<output physicell_name="apoptosis" intracellular_name="Resist_apop">`
  (action=activation, value=0 when ON, base_value=`--smc-apoptosis-rate`
  when OFF) on BOTH cell types -- replacing the two separate
  `--smc-apoptosis-rate`/`--founder-apoptosis-rate` CLI constants from
  section 10.13 with a single continuous, per-agent, Boolean-network-driven
  mechanism.
- All SMC-lineage agents now seed as `Smooth_Muscle_Cell` only
  (`--steap4-fraction 0.0` in `build_vessel_cross_section.py`) -- the literal
  implementation of "no more clones": `Synthetic_SMC` cells exist ONLY via
  transformation during the run.

**A fourth real PhysiCell bug found and fixed this session**: the first
validation run showed ZERO cells ever transforming into `Synthetic_SMC`
despite Klf4 reaching 87% ON in surviving `Smooth_Muscle_Cell` agents.
Tracing the issue to ground truth (not guessing) by testing an aggressive
0.05/min transformation rate revealed cells WERE transforming -- into
**`Endothelial_Cell`**, silently. Root cause, confirmed by reading
`core/PhysiCell_cell.cpp` and `core/PhysiCell_signal_behavior.cpp` directly:
PhysiCell maintains TWO parallel cell-type indexing schemes that are
supposed to agree but don't when a cell_definition's declared XML `ID`
attribute doesn't match its file declaration position --
`cell_definitions_by_index` (a plain vector, populated by constructor-time
`push_back`, so indexed by FILE ORDER) vs `cell_definitions_by_type` (an
`unordered_map<int,Cell_Definition*>` keyed by the literal `ID` attribute,
used to build the GLOBAL "transform to X" behavior dictionary via
`cell_definitions_by_type[i]` for `i=0..N-1`). This project's generator
declared `Smooth_Muscle_Cell` (ID=0), `Synthetic_SMC` (ID=3, a historical
artifact of STEAP4hi being bolted on after Fibroblast=1/EC=2 were already
assigned), `Fibroblast` (ID=1), `Endothelial_Cell` (ID=2) -- ID 0,3,1,2 in
FILE POSITIONS 0,1,2,3. The behavior-dictionary WRITE path resolves
"transform to Synthetic_SMC" to slot 3 (its ID), but
`standard_cell_transformations`'s READ path (`core/PhysiCell_cell.cpp`)
interprets slot 3 via `cell_definitions_by_index[3]` -- FILE POSITION 3,
which is `Endothelial_Cell`. Any PhysiCell project using `<cell_transformations>`
with non-sequential `ID` attributes (IDs that don't match declaration order)
will silently convert cells to the WRONG type with no error or warning.
**Fixed** by reordering the cell_definition declaration sequence so each
type's `ID` matches its file position exactly (`Smooth_Muscle_Cell`=0,
`Fibroblast`=1, `Endothelial_Cell`=2, `Synthetic_SMC`=3, in that file order)
-- verified correct by confirming newly-transformed agents carry the
`Synthetic_SMC` cell type AND an intact MaBoSS Boolean state (prior to the
fix, "transformed" agents landed in the EC bucket with NO intracellular
state at all, since EC has no `<intracellular>` block).

**Validation (7-day run, default rates, f_snx29=1.0, chronic AMLO)**:
`Smooth_Muscle_Cell` declines 767 -> 238 while `Synthetic_SMC` emerges
0 -> 414 -- entirely through Boolean dynamics, no population ever seeded
directly. Total SMC-lineage population only declines mildly (1186 -> 1071,
~10%) despite the ongoing apoptosis baseline, consistent with Resist_apop
rising to 70-79% in both surviving populations by day 7 (apoptosis-prone
cells are selectively lost; the rest increasingly resist). An unplanned but
biologically sensible spatial pattern emerged: `Synthetic_SMC` forms a dense
ring immediately adjacent to the endothelium, with regular
`Smooth_Muscle_Cell` pushed to the outer media -- consistent with the
AMLO/ET1/hypoxia gradients driving the Nfat/Klf4 switch earliest closest to
the lumen.

Figure: `v8_transformation_validation.png` (population dynamics + day-7
spatial snapshot). Data: `v8_transformation_population_timeseries.csv`.
Model: `ModeloPASMC_v8.bnd`/`.cfg`. Bundle: `v8_validate_bundle.tar.gz`.

**Caveats**: `--transformation-rate`/`--reverse-transformation-rate` are
uncalibrated placeholders, matching the network file's own "POR CALIBRAR"
flags on the Klf4/Myocd toggle parameters ($u_klf4, $d_klf4, $leak_klf4).
Single-seed run. The EC-barrier/chemotaxis/SNX29-sweep machinery from
sections 10.12-10.13 is unaffected by this change (still operates on
whatever agents currently hold the `Synthetic_SMC` type, regardless of
whether they arrived there by seeding or transformation) but has not yet
been re-run end-to-end against v8 -- the SNX29xAMLO sweep in section 10.13
still used the pre-v8 fixed-fraction model.


### 10.15 Per-agent continuous heterogeneity: real STIM1/ORAI1/ATP2A2 distributions, no new C++

Resumed the "SMC agent representation" discussion from earlier: each
Smooth_Muscle_Cell/Synthetic_SMC agent now gets its OWN STIM1/ORAI1/ATP2A2
level, sampled independently from the REAL snRNA-seq-fitted Beta_minmax
distributions in `smc_heterogeneity_summary.csv` (PAH condition, n=3510
real SMCs) -- replacing the prior all-agents-identical istate (`=1`,
deterministic, for all three). CACNA1C is deliberately excluded: v8 made it
Myocd-derived, not a frozen input, so personalizing it directly would fight
the new toggle's own dynamics.

**Implementation, confirmed to need zero new C++** (settling the question
from earlier this session about whether a C++ istate-override hook was
necessary): verified via direct source inspection that `MaBoSSIntracellular
::update_inputs()` calls `PhysiCell::get_selected_signals()`, the SAME
generic signal dictionary that resolves `custom:<name>` custom_data values
-- meaning a per-agent custom_data value can drive a Boolean node's input
exactly like the existing AMLO/ET1/HIPOXIA/K_ALTO mappings. Real per-gene
levels are sampled in Python (`sample_gene_levels`, matching the fitted
Beta(alpha,beta) on [min,max]) and written directly into cells.csv as three
new v2-format columns (`custom:stim1_level`, `custom:orai1_level`,
`custom:atp2a2_level`, populated only for SMC-lineage rows), then mapped
onto their corresponding Boolean nodes via
`<input physicell_name="custom:...">` (action=activation, threshold=0.5)
-- gated behind a new `--per-agent-heterogeneity` flag, off by default
(confirmed byte-identical XML output to all prior runs when unset).

**A fifth real PhysiBoSS bug found and fixed**: the third (and only the
third) custom-data-sourced `<input>` in a cell_definition's mapping
consistently read as if its signal were always below threshold, REGARDLESS
of which gene occupied that slot, its variable name, or its position
relative to the pre-existing AMLO/ET1/HIPOXIA/K_ALTO inputs -- isolated
through a systematic elimination: reordering which gene was 3rd (always
that one failed, not a specific gene), renaming the custom_data variable
in that slot (still failed), and moving the whole group to a different
place in the mapping (still failed, as long as it stayed 3rd within its
own group). The deciding experiment: adding a 4th, genuinely harmless
custom-data input after the real 3 made ALL THREE real ones work
correctly. Further isolated to: a custom_data variable that is declared in
the cell_definition's `<custom_data>` block but never actually written via
a `cells.csv` column (relying on its XML default) does not get correctly
wired into the per-agent signal-reading path when it is the last such
input -- confirmed by adding an actual (dummy, constant 0.0) `cells.csv`
column for the padding variable too, which fixed it. **Fixed** with a
permanent workaround: `heterogeneity_inputs_xml`/`heterogeneity_custom_data_xml`
always append a disposable 4th entry (`custom:padding_unused`, mapped to
the otherwise-unused `CPA` node with a threshold of 999999 so it can never
actually fire), absorbing the bug and leaving the 3 real genes correctly
read.

**Validated**: at the first intracellular update after seeding, observed
per-agent ON fractions (STIM1 15.4%, ORAI1 2.1%, ATP2A2 24.0%) closely
track the real data's expressing-cell fractions (STIM1 20.5%, ORAI1 1.3%,
ATP2A2 32.2%) -- a reasonable approximation given the simple >0.5 raw-count
threshold used, not an exact reproduction of the discrete zero-inflation a
continuous Beta fit cannot capture exactly. The underlying sampler itself
was separately confirmed to reproduce each gene's real fitted Beta density
exactly (log-density overlay, n=5000 agents).

Figure: `smc_per_agent_heterogeneity.png` (sampler validation + per-agent
ON-fraction vs real sparsity). Bundle: `het_final_bundle.tar.gz`.

**Caveats**: single seed; the 0.5 raw-count activation threshold is a
simple choice, not independently calibrated per gene; this makes
STIM1/ORAI1/ATP2A2 continuously forced toward their per-agent level for the
agent's entire life (a sustained bias), not a one-time t=0 draw that the
network's own logic is then free to evolve away from -- a real, disclosed
design tradeoff versus the (more complex, C++-requiring) per-instance
istate-override alternative considered and set aside earlier in this
session.

### 10.16 Validating against real clinical PAH markers

Computed real Control-vs-PAH log2 fold-change for all 14 genes in
`smc_heterogeneity_summary.csv` (n=1205 Control / 3510 PAH real SMCs,
Beta_minmax fits) and compared the *direction* of change against the
published PAH literature for each gene's corresponding model node.

**Matches** (model's input-personalization and/or node logic is consistent
with real data *and* literature):
- `HIF1A` up (log2FC +0.458, 9.9%->12.4% of cells expressing) -- the
  pseudohypoxic hallmark.
- `NFATC3` up (log2FC +0.193, 7.0%->8.4%) -- Bonnet 2007 calcineurin/NFAT
  activation.
- `TRPC6` up (log2FC +0.179, 9.5%->9.9%) -- Kuwahara 2006 and others.
- `MYLK` down (log2FC -0.491, 53.9%->48.2%) -- contractile->synthetic
  de-differentiation, directly consistent with the v8 Myocd/Klf4 toggle this
  project just added.

**Mismatches** (real-data direction contradicts the literature):
- `KCNA5` (Kv1.5) **up** in this cohort's real snRNA-seq (log2FC +0.869,
  12.0%->21.3% expressing), directly contradicting the robustly-established
  "Kv channelopathy" (Kv1.5 *down*regulation) literature confirmed by an
  independent web search across multiple studies. Treated as a
  data-interpretation caveat (cohort-specific, or a snRNA-seq-transcript-vs-
  ion-channel-protein/function disconnect -- Kv1.5 availability is heavily
  post-translationally regulated), not a modeling bug: `Kv1_5` in
  ModeloPASMC_v8 is governed purely by network logic (`NOT Nfat AND NOT
  HIPOXIA`), never personalized from transcript data, so this mismatch
  cannot corrupt any simulation result in this project.
- `PPP3CA` (calcineurin) and `RHOA` mRNA down (log2FC -0.120, -0.893) vs. the
  literature's pathway-*activity* up -- plausibly explained by
  post-translational/post-transcriptional activation without a transcript
  increase (calcineurin and ROCK are both regulated well beyond mRNA level);
  softer mismatch than KCNA5 since the literature claim is about activity,
  not expression.
- `ORAI1` down (log2FC -1.057, 3.2%->1.3% expressing) vs. mixed
  literature (some SOCE-proliferation studies report ORAI1 up).

**Quantitative implication acted on in this project**: `ORAI1` is expressed
in only ~1.25% of real PAH SMCs (vs `CACNA1C` ~46.6%) -- since
`Soce = Stim1_pm AND ORAI1`, pre-heterogeneity simulations (deterministic
`ORAI1=1` for 100% of agents) likely overestimated SOCE prevalence. This
was the direct motivation for the per-agent heterogeneity work in section
10.15 and the sweep rerun in section 10.17 below, both of which use the
real ORAI1 sparsity instead of a uniform istate.

Figure: `pah_marker_validation.png`. Table: `pah_marker_validation_table.csv`.

### 10.17 CONTR_1 -> geometry: reversible contraction, clinically-translatable metrics, and the re-run sweep

`CONTR_1` (myosin phosphorylation / MLCK-MLCP balance, requiring `Myocd` --
only contractile-phenotype SMCs can do this) was previously a readout-only
node. Real vasoconstriction is a distinct, reversible, functional narrowing
mechanism separate from the slow/irreversible *structural* narrowing already
modeled via Synthetic_SMC accumulation -- so it needed to map onto actual
vessel geometry, not just be read off the Boolean state.

**Mechanism** (`apply_contraction_pull` + `apply_lumen_barrier`, both new,
small, `update_velocity`-pattern C++ functions, same precedented pattern as
the already-validated `founder_barrier_velocity`):
- While `CONTR_1=1`, a small inward radial velocity (`--contraction-speed`,
  default 0.5 um/min, uncalibrated placeholder) is added to the agent's
  velocity, pulling it toward the vessel centerline.
- While `CONTR_1=0`, no inward term is added -- existing repulsion/pressure
  mechanics relax the agent back outward on their own. This is a VELOCITY
  term with **no accumulator or hysteresis**: the pull is recomputed from
  the live Boolean state every mechanics step.
- Wired to **both** `Smooth_Muscle_Cell` (new `contraction_velocity`
  function) and `Synthetic_SMC` (added into the existing
  `founder_barrier_velocity`), sharing one `apply_lumen_barrier` helper that
  keeps both from ever crossing into the lumen (basement
  membrane/IEL, impenetrable to mechanical force of any kind, same as
  before). The `Endothelial_Cell` ring itself is NOT clamped -- it is free
  to be mechanically compressed inward by pressure from the SMC-lineage mass
  outside it, which is exactly how narrowing is read out.

**A real bug found during validation testing**: without the lumen-barrier
extension, contracting `Smooth_Muscle_Cell` agents (which previously had no
barrier at all, only `Synthetic_SMC`/founder cells did) tunnelled straight
through the endothelial ring under sustained contraction (observed minimum
radius 3.4 um, i.e. essentially at the vessel centerline) -- unphysical, since
media SMCs should never enter the flowing-blood lumen. Fixed by extracting
the barrier clamp into a shared `apply_lumen_barrier` applied to both cell
types.

**Validation** (`contraction_mechanism_validation.png`): (a) population
level -- `--contraction-speed 0.5` compresses the EC lumen (5th-pct radial
position) from 104.5 to ~65-67 um over 1-2 days, while `--contraction-speed
0.0` leaves it completely unchanged (104.5 -> 104.5); (b) per-agent level --
individual `Smooth_Muscle_Cell` agents' radial position tracks their own
live `CONTR_1` state with no hysteresis (e.g. one agent: r=100.0 um while
`CONTR_1=off` -> r=117.9 while `CONTR_1=ON` is immediately followed by
r=138.3 once `CONTR_1` turns `off` again -- genuine relaxation, not a
one-way ratchet).

**A second real bug found while re-running the sweep**: the
`--amlo-secretion-rate 0.02` value carried over from the earlier
(base_v2-era) SNX29xAMLO sweep never crossed the AMLO Boolean activation
threshold (0.3) anywhere in the current v8 `cells_pah_v8.csv` geometry --
the microenvironment field plateaued at ~0.032-0.035 domain-wide (diffusion
length `sqrt(D/decay) = sqrt(4000/0.002) ~ 1414` um vastly exceeds the
~840 um domain, so the field equilibrates to a low uniform level set by
total secreting area vs. decay, not by distance from the source). Both
"novehicle" and "chronic" arms were silently running drug-free. Rate-scanned
{0.1, 0.5, 2.0, 10.0} and selected **0.5**, which produces a near-uniform,
above-threshold field (~0.40-0.43) across the whole domain within ~1 day.
This is a deliberate deviation from section 10.8's earlier "shallow
penetration" framing (where AMLO stayed sub-threshold beyond ~160 um of the
lumen) -- at this geometry and these diffusion parameters, the drug field
cannot be made spatially graded and clinically relevant (i.e. crossing
threshold at all) simultaneously without further retuning of
diffusion/decay, which was out of scope here. The choice is best read as "a
systemic, whole-wall-penetrating chronic dose" rather than a locally-limited
one.

**Clinically-translatable metrics** (`clinical_metrics.py`):

| Metric | Definition | Clinical correlate |
|---|---|---|
| Effective lumen radius | 5th-pct `Endothelial_Cell` radial position | Pathology luminal-occlusion grading |
| % diameter stenosis | `100*(1 - r_t/r_0)` | Angiographic/pathology stenosis grading |
| Relative PVR proxy | `(r_0/r_t)^4` (Poiseuille, fixed flow/viscosity/length) | RHC-measured PVR (one-segment proxy, NOT whole-lung) |
| `Soce_frac` / `Nfat_frac` | fraction of SMC-lineage agents with node ON | Molecular pathway activation |
| `CONTR_1` fraction | fraction of SMC-lineage agents with node ON | Acute vasoreactivity testing (Sitbon/Humbert CCB-responder criteria) |
| AMLO efficacy index | `1 - (metric_chronic / metric_novehicle)` at matched time, computed for both CONTR_1 and stenosis | Clinical drug response magnitude |
| AMLO escape index | `efficacy_acute(day1) - efficacy_chronic(day7)` | Erosion of benefit, cf. "time to clinical worsening" trial endpoints |

**Re-run sweep results** (all 6 arms: {low=0.5, normal=1.0, high=2.0}
`f_snx29` x {novehicle, chronic=0.5} AMLO, now with
`--per-agent-heterogeneity` and `--contraction-speed 0.5` on top of the
section 10.13/10.14 Nfat-pathway-gated, apoptosis-corrected, v8 config; 7-day
runs, single seed):

- Day-7 novehicle stenosis is essentially flat across SNX29 levels
  (37.4-38.1%), as is day-7 chronic stenosis (28.0-29.1%) -- chronic AMLO
  reduces stenosis by ~9-10 percentage points at every SNX29 level in this
  single-seed run.
- `Nfat_frac` falls from ~50-53% (novehicle) to ~18-22% (chronic) at day 7;
  `CONTR_1_frac` falls from ~17-18% to ~12-15%.
- **AMLO efficacy (CONTR_1-based) erodes with chronic dosing at every SNX29
  level**: low 0.541->0.246, normal 0.603->0.107, high 0.443->0.315 (acute
  day-1 -> chronic day-7). The resulting **escape index** is lowest at high
  SNX29 (0.128) and highest at normal SNX29 (0.496, low SNX29 0.294) --
  directionally consistent with the SNX29-brake hypothesis at the high end
  (more brake = less SOCE-driven escape) but non-monotonic between low and
  normal, which single seed per arm cannot distinguish from noise.
- The stenosis-based efficacy/escape indices are much smaller and in two of
  three arms slightly *negative* (sustained or marginally improved benefit:
  normal -0.022, high -0.059) -- because structural narrowing (accumulated
  `Synthetic_SMC` mass, slow-changing over 7 days) dominates the stenosis
  metric and dilutes the functional-contraction-specific escape signal that
  `CONTR_1` captures directly. This divergence is itself informative: a
  purely structural/hemodynamic readout (closer to what RHC-measured PVR
  would show) would understate the erosion of the drug's acute functional
  benefit that a vasoreactivity-style readout (`CONTR_1`) reveals.

Figure: `snx29_amlo_clinical_metrics.png`. Tables:
`snx29_amlo_clinical_metrics.csv` (endpoint + efficacy/escape summary),
`snx29_amlo_clinical_metrics_timeseries.csv` (full per-frame series).
Bundle: `clin_sweep_bundles.tar.gz`.

**Caveats**: single idealized 2D vessel cross-section, not whole-lung
hemodynamics -- the relative-PVR proxy reflects one vessel segment's
resistance contribution, not a calibrated absolute mPAP/PVR prediction;
`contraction_speed` and the AMLO efficacy/escape magnitudes are uncalibrated
placeholders (no published PhysiCell-scale SMC contraction velocity or
Boolean-to-geometry mapping exists to fit against); single random seed per
arm throughout -- the non-monotonic low-vs-normal SNX29 escape ordering
should not be over-interpreted without replicate seeds; the AMLO dose
(0.5) was chosen to produce a detectable above-threshold effect in this
geometry, not derived from a pharmacokinetic model.

### 10.18 Correction to section 10.16/10.17: ORAI1 dropout, and why the spatial sweep's escape-kinetics fix was abandoned

**ORAI1 is NOT rare -- the ~1-3% detection was snRNA-seq dropout, not biology.**
Section 10.16 flagged ORAI1's real-data detection rate (~1.3-3.2% of SMCs) as a
possible mismatch; follow-up established this more firmly from two
independent angles: (a) the PAH literature treats STIM1/Orai1-mediated SOCE
as a broadly-active, often-*upregulated* pathway in PASMC -- a central
mechanism, not a rare-subpopulation one; (b) ORAI1's gene length (16.6 kb) is
14-44x shorter than STIM1 (238.7 kb), ATP2A2 (70.5 kb), and CACNA1C (727.3
kb) -- single-**nucleus** RNA-seq relies heavily on intronic reads for
capture, so a gene with a small intron budget is structurally disadvantaged
for detection regardless of true expression. Note: the raw `PAH_snRNAseq.h5ad`
was never actually available to this session (only a pre-computed, 14-gene
summary CSV was shared by a collaborator) -- the exact requested QC
(length/GC-matched housekeeping-gene comparison in the same dataset) could not
be run; the gene-length argument above is the best available substitute.

**Fix**: ORAI1 keeps per-agent heterogeneity (not deterministic ON=1 for
everyone either) but with a corrected baseline. A smooth, unimodal
Beta(alpha=3, beta=2, vmin=0, vmax=1) -- no boundary pile-up, unlike the
original fit's very small alpha/beta -- solved to give ~65-69% of agents ON
at the shared 0.5 raw-count threshold (validated: 69.6% sampled in a test
run), sitting in the literature's "broadly active" regime instead of the
dropout-suppressed ~1-3%. This is a deliberately uncalibrated placeholder (no
paper reports "fraction of PASMC with functional SOCE") -- chosen only to
correct the *order of magnitude*, not to claim a specific true value.
`REAL_GENE_BETA_PARAMS_PAH` in `generate_physiboss_vessel_gradient_config.py`
now carries STIM1, ATP2A2, and this corrected ORAI1.

**Why this mattered**: with the original ~1-3% ORAI1, `Soce = Stim1_pm AND
ORAI1` was capped near zero regardless of `Stim1_pm`'s own state -- a direct
test (forcing `Stim1_pm`'s rate constants down ~100,000x via a new
`--escape-rate-scale` CLI arg) showed `Nfat`/`CONTR_1`/stenosis were
*unchanged* whether SOCE was fully active or fully suppressed, proving
`Ca_ltcc`/`Ca_trpc6` (not SOCE) dominated those readouts under the original
heterogeneity. With ORAI1 corrected, `Soce` tracks `Stim1_pm` almost 1:1
again, restoring real flux for the SNX29 brake to act on.

**Why the spatial sweep's "escape over time" framing was abandoned (this
session) rather than fixed**: chasing a genuine, gradually-escalating escape
signal in the spatial PhysiCell model surfaced three more layers of
complication, each real and each documented here for future reference:

1. **MaBoSS/PhysiCell time-compression.** PhysiBoSS runs a *full* MaBoSS
   simulation out to `max_time=200` internal units on *every*
   `intracellular_dt=60`-minute PhysiCell step. A rate constant needs
   `rate * (max_time/intracellular_dt) * T_target << 1` to represent a
   genuinely slow, multi-day process -- for a ~20-30 day timescale this means
   rates on the order of `1e-5`, not the ~10-100x reduction that would be the
   naive guess from the nominal rate values alone.
2. **The calcium-signaling layer is deliberately self-limiting.** Even with
   ORAI1 corrected and `Stim1_pm`'s rates scaled to the correct order of
   magnitude, `Soce` rose then *decayed* back toward zero over ~20 days under
   chronic AMLO, rather than escalating. Root cause: `Serca`'s refill
   capacity is gated only by each agent's own (fixed, AMLO-independent)
   `ATP2A2` level -- `Sr_lleno` (the SR calcium store) refills regardless of
   sustained Ltcc blockade, closing `Stim1_pm`'s gate again. The network's
   own comment on the sibling node `Calcineurina` states this is intentional
   ("R2: NO ABSORBENTE" -- phosphatase activity decays when calcium drops, it
   doesn't stay fixed): the calcium layer is built as a non-ratcheting,
   self-limiting system throughout, not a one-way escalation.
3. **A CaMKII/phospholamban-motivated fix was tried and reverted.** Gating
   `Serca`'s `rate_up` on `AMLO` (chronic Ltcc blockade -> less CaMKII
   activity -> phospholamban stays inhibitory -> degraded SERCA refill) had
   *zero* effect, because `Serca.istate=1` and its `logic` (`ATP2A2 AND NOT
   CPA`) is permanently true for personalized agents -- making `rate_down=0`
   always, an absorbing ON-state from t=0 that a `rate_up` edit can never
   reach. A real fix needs `rate_down` made AMLO-responsive too, a materially
   larger change than approved; reverted rather than compounding further
   unreviewed edits to a shared model file.

Given four non-trivial, interacting findings in one investigation thread, the
spatial PhysiCell sweep's "escape index" framing (section 10.17) is
**withdrawn as a demonstrated result** -- not because the underlying
SNX29/AMLO mechanism is wrong, but because reproducing it cleanly in the
spatial model's current parameterization needs more rework (the `rate_down`
fix above, re-validated geometry/contraction interactions, a new sweep) than
fits in one session on top of everything else already changed. Section 10.19
below demonstrates the same mechanism in a simpler, well-mixed setting where
it resolves cleanly.

### 10.19 A well-mixed toy demonstration of the SNX29 x amlodipine mechanism

Stepping back from the spatial model's accumulated complexity, this directly
tests the core hypothesis -- does SNX29's brake strength change how much
amlodipine drives compensatory `Stim1_pm` activation -- using
`ModeloPASMC_v8.bnd`/`.cfg` with the **standalone MaBoSS binary** (no
PhysiCell, no spatial geometry, no per-agent sampling): 3 SNX29 levels
(`$f_snx29`: low=0.5, normal=1.0, high=2.0) x 2 AMLO conditions (off/on) = 6
runs, each an ensemble of 1000 independent stochastic trajectories (the
cfg's own `sample_count` default) giving a clean marginal P(node=1)
trajectory over MaBoSS's own internal time (max_time=500, uncalibrated to
real clinical time). `ORAI1.istate` uses the corrected population split
(31% off / 69% on, matching section 10.18's Beta(3,2) target) via MaBoSS's
probabilistic-istate syntax (`[ORAI1].istate = 0.31 [0], 0.69 [1];` --
**note the required square brackets around the node name on the LHS**, a
MaBoSS syntax quirk documented earlier in this project).

**Result** (`toy_snx29_amlo_demo.png`, `toy_snx29_amlo_summary.csv`):

- `Stim1_pm`'s baseline level (independent of AMLO) decreases monotonically
  with SNX29 level, exactly as the brake mechanism is designed: low
  SNX29 ~0.066-0.076, normal ~0.043-0.054, high ~0.024-0.037 -- a direct,
  clean confirmation that `$f_snx29` divides `Stim1_pm`'s rate_up as
  intended.
- The single-timepoint AMLO-minus-novehicle delta is noisy at the raw
  sampling resolution (SEM ~0.01 on probabilities ~0.05-0.1 from a
  1000-trajectory ensemble) -- windowed averaging resolves it cleanly.
  Averaged over t=0-50 ("acute") and t=400-500 ("chronic"):

  | SNX29 level | acute \u0394P(Stim1_pm) | chronic \u0394P(Stim1_pm) | cumulative \u222b\u0394 dt (full window) |
  |---|---|---|---|
  | low (weak brake) | +0.0015 \u00b1 0.0003 | +0.0018 \u00b1 0.0005 | 0.940 |
  | normal | +0.0008 \u00b1 0.0002 | +0.0016 \u00b1 0.0003 | 0.521 |
  | high (strong brake) | +0.0006 \u00b1 0.0001 | +0.0010 \u00b1 0.0002 | 0.233 |

  Every entry is several SEM from zero. The effect is monotonic in SNX29
  level at both timepoints, and grows from acute to chronic at every SNX29
  level -- a clean, well-resolved demonstration that **AMLO drives more
  compensatory Stim1_pm activation when SNX29's brake is weaker, and this
  effect accumulates with sustained exposure**, exactly the mechanism this
  whole project is built around, now shown without the spatial model's
  confounds.

**Caveats**: well-mixed (no spatial/mechanical structure at all -- not a
replacement for the spatial stenosis/contraction work, a complementary,
cleaner mechanistic check); MaBoSS's own internal time units, uncalibrated to
real clinical time (same caveat as throughout this project); single ensemble
per arm (no replicate seeds -- the SEMs above are the ensemble's own sampling
SEM, not a between-seed variance); `Stim1_pm` itself (not a downstream
functional readout like `CONTR_1` or stenosis) is the chosen readout
specifically because it is the most direct target of `$f_snx29` and was
shown to carry the cleanest signal -- downstream readouts dilute it through
additional pathways (`Ca_ltcc`, `Ca_trpc6`) not gated by SNX29 at all.

Script: `toy_snx29_amlo_experiment.py`. Data:
`toy_snx29_amlo_timeseries.csv` (24,000 rows: 6 arms x 4 tracked nodes x
1000 timepoints), `toy_snx29_amlo_summary.csv`.

### 10.20 A well-mixed demonstration of the SMC phenotype transition (Myocd<->Klf4 toggle)

A second, complementary well-mixed (standalone MaBoSS) demonstration, this
time of the *other* major mechanism in ModeloPASMC_v8: the `Myocd`/`Klf4`
double-negative toggle that drives the contractile-to-synthetic SMC
phenotype switch (the Boolean logic behind the `"transform to Synthetic_SMC"`
native PhysiCell behavior validated spatially in section 10.14).

**Design**: does a *transient* depolarization pulse (`K_ALTO`, representing
an acute vasospasm/40 mM-KCl-type event) commit a cell irreversibly to the
synthetic phenotype, or does it revert once the pulse ends? Standalone
MaBoSS has no native time-varying input within one run, so this chains two
sequential simulations: **Stage A** runs with `K_ALTO.istate=1` (pulse ON,
all else at resting istates) to `max_time = pulse_duration`, and extracts the
marginal P(node=1) for every node in the network at that instant. **Stage
B** runs with `K_ALTO.istate=0` (pulse OFF) and every other node's istate
set to Stage A's end marginal (probabilistic istate, treating the joint
state as independent per-node margins -- a standard approximation for
chaining MaBoSS runs), for a long (400-unit) observation window. Swept
`pulse_duration` over {0, 5, 10, 15, 20, 30, 40, 60, 80, 120, 160, 200}.

**Result** (`toy_phenotype_toggle_demo.png`,
`toy_phenotype_toggle_results.csv`): two findings, one more interesting than
the bistable-switch hypothesis this was designed to test:

1. **The resting state itself is not a stable attractor over long
   timescales.** Starting fully contractile (`Myocd=1`, no pulse at all,
   `pulse_duration=0`) and simply observing for 400 units, the population
   spontaneously drifts from `Myocd=1.0`/`Klf4=0.0` at t=0 to
   `Myocd=0.46`/`Klf4=0.54` by t=400 -- already past the 50/50 point, with no
   external perturbation whatsoever. This is itself a notable property of the
   current network parameterization: the "contractile" resting state is
   metastable (slowly draining), not a hard fixed point, at these (uncalibrated)
   basal rate constants.
2. **The K_ALTO pulse accelerates this drift rather than creating a sharp,
   separate attractor.** Longer pulses push the Stage-B starting point
   further toward `Klf4` (e.g. `pulse_duration=200` starts Stage B at
   `Myocd=0.064`), but all trajectories partially reconverge over the
   400-unit window rather than snapping to two cleanly separated endpoints
   -- the final `Myocd`/`Klf4` vs. `pulse_duration` curve is **graded and
   monotonic, not step-like**. End-of-window values range smoothly from
   `Myocd=0.456`/`Klf4=0.544` (no pulse) to `Myocd=0.242`/`Klf4=0.757`
   (200-unit pulse), saturating at long pulse durations rather than
   exhibiting a sharp threshold.

This is an honest, more nuanced finding than "clean bistable switch": the
pulse's role here is to give the synthetic fraction a head start on an
ongoing, graded transition that the network already undergoes at baseline,
not to flip a hard toggle with two well-separated stable states. Caveats:
uncalibrated MaBoSS time units (same caveat as section 10.19); single
ensemble per condition (no replicate seeds); the independent-marginals
approximation used to chain Stage A into Stage B's istates discards
joint-state correlations, a standard but real approximation.

Script: `toy_phenotype_toggle_experiment.py`. Data/figure:
`toy_phenotype_toggle_results.csv`, `toy_phenotype_toggle_demo.png`.

**Companion spatial visualization** (`smc_phenotype_transition_spatial.gif`):
the well-mixed toggle demo above shows *that* the transition happens and
under what conditions; the already-validated v8 spatial run (section 10.14,
`v8_validate_fixed`, re-visualized here with no new simulation) shows
*where*. A 29-frame time-lapse (day 0 to day 7, phenotype-colored:
`Smooth_Muscle_Cell` blue, `Synthetic_SMC` red, `Endothelial_Cell` orange,
`Fibroblast` grey) makes the spatial pattern immediately legible: a dense
ring of `Synthetic_SMC` forms directly against the endothelium and thickens
over the week, pushing the remaining contractile `Smooth_Muscle_Cell`
population outward toward the adventitia -- the same emergent pattern
reported qualitatively in section 10.14, now the explicit subject of the
figure. Counts verified unchanged from the original validation (767->238
`Smooth_Muscle_Cell`, 0->414 `Synthetic_SMC`, `Fibroblast`/`Endothelial_Cell`
constant at 367/52). **v2**: cells rendered as true-radius circles (7.5 um
for `Smooth_Muscle_Cell`/`Fibroblast`/`Synthetic_SMC`, 6.0 um for
`Endothelial_Cell`, read directly from the PhysiCell volume parameters --
1767.15 / 904.78 um^3 respectively -- rather than fixed-size scatter dots)
so the packed, touching tissue layer and the degree of lumen crowding by the
`Synthetic_SMC` ring are both visible at their correct relative scale.

**Companion substrate-field GIF** (`substrate_fields_spatial.gif`): the same
29 frames (day 0 to day 7, same `v8_validate_fixed` run), now showing the
diffusible-substrate side of the model instead of the cells -- a 2x2 panel
grid (`AMLO`, `ET1`, `O2`, `K_ALTO`), each a heatmap on the native 42x42,
20 um-spaced voxel grid, with a FIXED color scale per substrate (set from
that substrate's own global min/max across all 29 frames, so color encodes
genuine temporal change rather than being re-normalized away each frame).
`AMLO` and `ET1` both build up visibly from zero: `ET1` (secreted by
`Endothelial_Cell`, consumed by the founder-cell chemotaxis mechanism)
forms a persistent, strongly-localized hot spot centered on the lumen,
while `AMLO` (secreted at this run's original `--amlo-secretion-rate 0.02`,
later shown in section 10.17 to be too low to cross the Boolean activation
threshold anywhere) stays well below 0.035 and spreads nearly uniformly
across the whole domain -- a direct visual confirmation of the diffusion-
length argument in section 10.17 (`sqrt(D/decay) ~ 1414 um >> ~820 um`
domain width). `O2` shows a visible local depletion patch forming where the
cell mass is densest. `K_ALTO` stays flat at exactly zero throughout --
correctly reflecting that this substrate is not used in this scenario
(secretion rate 0), not a rendering artifact.

Figures/GIFs this section: `smc_phenotype_transition_spatial.gif` (cells,
true radius), `substrate_fields_spatial.gif` (substrates, fixed per-field
color scale) -- same underlying run, same 29 frames, complementary views.

### 10.21 Steeper AMLO/O2 gradients, and fixing the unrealistic media-layer "holes"

Two requested corrections to the spatial example, both now the generator
script's defaults:

**Steeper substrate gradients.** `AMLO` (diffusion_coefficient 4000->600,
decay_rate 0.002->0.03) and `O2` (diffusion_coefficient 100000->20000,
decay_rate 0.1->0.4) were re-tuned so each substrate's decay length
(`sqrt(D/decay)`) sits within the vessel-wall thickness (~141 um for AMLO,
~224 um for O2) instead of vastly exceeding the ~820 um domain -- giving a
real radial gradient across the SMC layer instead of near-uniform
saturation. Verified at day 7: `AMLO` now declines ~9x from center to edge
(vs ~6% in the original parameterization); `O2` shows real radial structure
peaking near the EC/lumen source and declining outward. New CLI args:
`--amlo-diffusion-coefficient`, `--amlo-decay-rate`,
`--o2-diffusion-coefficient`, `--o2-decay-rate`. Both figures/GIFs in this
section were regenerated with these settings.

**Caveat**: the steeper O2 parameters, combined with the (much larger, see
below) cell population now consuming it, deplete `O2` severely almost
everywhere (~0.02-0.6 mmHg vs the ~38 mmHg EC-side target) rather than
producing a graded high-to-low profile across a physiological range -- real
structure is present (the gradient direction and shape are correct) but the
absolute magnitude may be over-corrected for a "mostly normoxic media,
hypoxic adventitia" picture. Worth revisiting if a specific physiological
O2 range (e.g. 5-40 mmHg) is wanted.

**No more unrealistic "holes" in the SMC layer.** The original
`Smooth_Muscle_Cell` division rate (`--division-rate`, 1e-5/min, ~48-day
doubling) was far slower than the ~7-day transformation/apoptosis turnover
driving the population down -- cells were lost to `Synthetic_SMC`
transformation or death much faster than remaining wild-type SMC could
proliferate to backfill the gaps, leaving a visibly sparser, holey media
layer by day 7 (`r_mean` thinning from 209 to 171 um over the run). Raised
the default to 2e-4/min (~2.4-day doubling). Isolated test (both runs with
`--contraction-speed 0`, since the CONTR_1 contraction mechanism from
section 10.17 defaults to an active 0.5 um/min that would otherwise confound
a pure division-rate comparison):

| | day 0 | day 7 (`Smooth_Muscle_Cell`/`Synthetic_SMC`) | `r_mean` day 0 -> day 7 |
|---|---|---|---|
| original (1e-5/min) | 767 / 0 | 238 / 414 (652 total) | 209 -> 171 um (thinning) |
| fixed (2e-4/min) | 767 / 0 | 478 / 668 (1146 total) | 209 -> 211 um (stable) |

The fixed run's total SMC-lineage population nearly doubles (1146 vs 652 at
day 7) and the tissue's mean radius stays essentially constant throughout
-- the media layer stays fully packed instead of thinning, visually
confirmed in the regenerated `smc_phenotype_transition_spatial.gif`
(compare the sparse, gap-containing original to the dense, continuous
fixed version).

Bundle: `v9_steep_validate_bundle.tar.gz`.

### 10.22 The 6-month PhysiBoSS campaign: periodic AMLO dosing, real O2 gradient, 4-gene heterogeneity, and a chain of mechanics bugs

This section documents the largest single buildout this project has done:
folding spatial gradients, per-agent cell heterogeneity, and a literature-
grounded drug-dosing regimen into one 3-6 month PhysiBoSS campaign.

**1. Literature-grounded periodic AMLO dosing.** Replaced AMLO's constant EC
secretion_rate with a one-compartment oral-absorption PK model (new
`update_amlodipine_pk` function in `custom.cpp`, called once per diffusion
step from `main.cpp`): `dA_gut/dt = -ka*A_gut`, `dC/dt = ka*A_gut - ke*C`,
redosed every `--amlo-dose-interval-min` (1440 = 24h default). `ka`/`ke` are
solved in Python at config-generation time from real amlodipine PK
(terminal half-life 30-50h -> `ke`; Tmax 6-12h -> `ka`, PMID 16296697 and
the FDA label), with `--amlo-dose-mg` (5/10/20mg) linearly scaling the dose
amount relative to a 10mg reference. Validated on a 14-day run: near-lumen
tissue AMLO concentration rises to steady state by day ~9-10 with a 24.3%
peak-to-trough swing, matching the literature's reported ~20-25% -- see
`amlo_pk_profile_validation.png`.

**2. Fixed the O2 radial gradient.** The Dirichlet-0 wall boundary condition
added in section 10.21 alone was not enough -- two further corrections were
needed: (a) the hardcoded per-cell O2 uptake rates (10/min SMC, 5/min
Fibroblast) gave a local depletion length (`sqrt(D/uptake)`) of only
~45-100 um, far shorter than the packed tissue's radial extent, crashing O2
to near-zero almost everywhere regardless of the boundary condition --
exposed as `--o2-baseline-uptake-rate` (1.5) / `--o2-fibro-uptake-rate`
(0.75); (b) only 52 `Endothelial_Cell` agents secrete O2 against a much
larger consuming population, so `--o2-secretion-rate` needed raising from
1.0 to 300 to approximate blood's effectively-unlimited O2 reservoir.
Combined with reverting `--o2-diffusion-coefficient`/`--o2-decay-rate` to
the standard PhysiCell tissue-oxygen values (100000, 0.1) now that the
Dirichlet wall + corrected uptake/secretion do the gradient-shaping work,
this gives a real, physiologically plausible profile: ~30-35 mmHg near the
lumen declining smoothly to ~4-5 mmHg at the outer wall, crossing the
HIPOXIA threshold (15 mmHg) partway through the media -- see
`o2_gradient_fix_validation.png`.

**3. A 4th personalized real gene: TRPC6.** Extended the per-agent
heterogeneity pipeline (previously STIM1/ATP2A2/ORAI1) with TRPC6, using
the real snRNA-seq PAH Beta_minmax fit (alpha=0.093, beta=4.656, vmax=15)
directly, WITHOUT the dropout correction ORAI1 needed -- TRPC6's gene span
(425 kb) is longer than STIM1's (238 kb), arguing against snRNA-seq dropout,
and the raw fit's implied ON-fraction (~13.6%) is already close to the
measured pct_expressing (9.9%). Required a genuine `.bnd` edit: added a new
frozen `TRPC6` gene-dosage node gating the dynamic `Trpc6` channel node
(`Trpc6.logic = TRPC6 AND (...)`, mirroring how `STIM1`/`ORAI1`/`ATP2A2`
already gate `Stim1_pm`/`Soce`/`Serca`) -- validated to be behaviorally
equivalent to the standalone-MaBoSS config-only knockout used earlier in
section 10.21's escape-route experiments.

**4. A chain of spatial-mechanics bugs, found by the user reviewing early
long-run output, not by any automated check:**

- **Tissue emptying under "no AMLO".** `--founder-division-rate` (Synthetic_
  SMC) defaulted to 1e-5 (~48-day doubling) while `--division-rate`
  (Smooth_Muscle_Cell) had been raised to 2e-4 (~2.4-day) in section 10.21 --
  once most agents transformed to Synthetic_SMC, the tissue's effective
  division capacity collapsed ~20x. Fixed: `--founder-division-rate` now
  defaults to `--division-rate`'s value unless explicitly overridden --
  there is no literature basis for synthetic SMC dividing slower than
  contractile SMC in situ.
- **Chemotaxis structurally enabled but functionally inert.**
  `--founder-migration-speed` defaulted to 0.0 -- direction without speed
  is a no-op. Raised to a modest 0.05 um/min; confirmed Synthetic_SMC mean
  radial position drifts inward over time (148.7 -> 138.7 um over 15 days
  in one validation run).
- **Transformation rate far too fast.** `--transformation-rate` (1.6e-4,
  ~3-day half-life once Klf4=1) transformed ~91% of the SMC lineage within
  3 weeks even with NO drug at all, given Klf4 reaches ~90%+ quickly near
  the ET1-rich lumen (chronic `Endothelial_Cell` ET1 secretion). Slowed
  ~15x to 1.0667e-5 (~45-day half-life), keeping the same 10x forward:
  reverse ratio -- now spreads gradually (0% -> 24% by day 30 in
  validation) instead of saturating in 3 weeks.
- **Local pressure-gating alone could not bound growth.** After fixing the
  above three, population exploded to ~24,000 agents by day 36 (from a
  previous steady state of ~1,700) -- the (accidental) brake from
  Synthetic_SMC dividing far slower than Smooth_Muscle_Cell had been doing
  double duty as a population limiter. Diagnosed with a custom
  `pressure_debug` custom_data field (which itself had a bug -- it
  captured something tracking simulation time, not real pressure) but
  cross-checked against PhysiCell's own standard `pressure` export field,
  confirming REAL geometric pressure stays genuinely low (mean 0.2-0.3,
  max 0.8->3.6 over 50 days) even as the population grows several-fold --
  a continuously-expanding low-pressure tissue EDGE keeps outrunning local
  suppression. Tried `pressure_half_max` 0.1 (delayed the explosion to
  day ~45-65) and 0.3 with a sharper Hill power (made it WORSE -- 10,304
  agents by day 25). Fixed with a second, GLOBAL Hill-type gate on total
  agent count (`--population-half-max`=2000, `--population-hill-power`=8),
  multiplicative with the local pressure gate -- local contact inhibition
  still shapes WHERE growth concentrates (the low-density edge), while the
  global gate caps HOW MANY cells can exist in total. Validated: genuine
  logistic-like saturation over 90 days (`population_cap_validation.png`),
  settling toward ~2,500 agents.

**5. The campaign.** Three 180-day (6-month) arms, each with the real
4-gene per-agent heterogeneity, daily `full_data` saves (~2.1-3MB/frame,
~475MB/run): `no_amlo` (natural history), `amlo_snx29_normal` (periodic
10mg dosing, `$f_snx29`=1.0), `amlo_snx29_low` (periodic 10mg dosing,
`$f_snx29`=0.5, less brake on SOCE). All three completed in ~22 minutes
each. Key results (`longrun_campaign_summary.png`):

- **Population** saturates almost identically across all three arms
  (~2,450-2,510 by day 180) -- the global population cap is robust to the
  AMLO/SNX29 perturbations tested.
- **Transformation fraction** (`Synthetic_SMC`/(`Synthetic_SMC`+
  `Smooth_Muscle_Cell`)) reaches 65.8% (no_amlo), 68.1% (amlo_snx29_normal),
  66.6% (amlo_snx29_low) by day 180 -- all three arms converge to a similar
  final extent, with only modest timing differences between them. Neither
  AMLO nor SNX29 level meaningfully changes the 6-month ENDPOINT, though
  they do shift the trajectory's timing slightly.
- **Radial profile at day 180 confirms the two-zone prediction** from
  section 10.20's mechanistic audit: a sharp inner zone (~100-180 um from
  the lumen) reaches 80-100% transformation, largely independent of AMLO
  or SNX29 (the ET1/TRPC6-driven route dominates there), dropping steeply
  to near-zero beyond ~200 um. The outer zone -- where SNX29-modulated
  SOCE-driven escape would be expected to matter, per the well-mixed
  ET1=0 experiment in section 10.21 -- simply has not had enough time to
  engage meaningfully even at 6 months, at the current (uncalibrated)
  rate constants. This is a genuine model prediction, not a null result:
  it says SNX29's modulatory effect on CCB resistance, if real, would
  only become visible over a LONGER follow-up than 6 months, or require
  faster outer-zone kinetics than currently calibrated -- both plausible
  refinements, not currently implemented.

**Caveats, explicitly**: the transformation-rate and population-cap
constants are tuned to produce qualitatively sensible, bounded dynamics --
they are not calibrated against any published rate. The PK model is a
standard one-compartment oral-absorption approximation, not a validated
amlodipine PK/PD model (no renal/hepatic clearance variation, no protein
binding, no explicit Cmax/Css in real concentration units). The population
cap is a deliberately crude GLOBAL backstop, not a fix to whatever
underlying local-mechanics subtlety let a continuously-expanding low-
pressure edge outrun local contact inhibition -- that root cause was not
fully resolved, only worked around robustly.

Deliverables: `amlo_pk_profile_validation.png`, `o2_gradient_fix_validation.png`,
`population_cap_validation.png`, `population_stability_check.png`,
`longrun_campaign_summary.png`, `longrun_6month_phenotype.gif`,
`longrun_no_amlo_bundle.tar.gz`, `longrun_amlo_snx29_normal_bundle.tar.gz`,
`longrun_amlo_snx29_low_bundle.tar.gz`.
