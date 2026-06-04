# Table & Figure Lab (Applied-Economics Exhibit Builder)

Generate the **full, methodology-aware sequence of tables and figures** for an applied-economics paper — not just one regression table, but the canonical exhibit set in the right order, in both **paper** and **presentation** variants, with an explicit **standard-error** choice. Emits Stata (`reghdfe`/`esttab`/`coefplot`) and R (`fixest`/`modelsummary`) plus AEJ booktabs LaTeX.

Companion interactive prototype (simulated data, live SE lab, paper-vs-presentation toggle):
**https://jkitmang.github.io/replication-katrina/**

This command **composes with** the existing skills — call them for the individual pieces:
- `reg-table` → the main regression / robustness / heterogeneity / IV tables (AEJ esttab templates)
- `balance-table` → the RCT balance table
- `attrition-check` → the attrition analysis
- `beamer-slides` → the presentation-variant figures/tables (colorblind-safe palette)

This command's job is the **orchestration**: which exhibits, in what order, for the chosen methodology; the figures `reg-table` does not cover (event-study, RD plot, first-stage scatter, forest plots); and the standard-error decision.

## Instructions

When the user runs this command, gather:

1. **Methodology**: DiD / event study · RCT · IV · RDD (drives the exhibit sequence).
2. **Key variables**: outcome(s), treatment/endogenous regressor, instrument, running variable, unit/time/cluster ids, strata, covariates.
3. **Clustering level** and whether **few clusters** (G ≲ 40 → wild cluster bootstrap).
4. **Target**: paper (`.tex` fragments, monochrome booktabs) and/or presentation (Beamer, colored, one-message-per-slide).
5. **Output dir** (default `output/tables` and `output/figures`).
6. **Language**: Stata (default) or R.

Then produce the exhibits **in canonical order**, each as its own file, with a short note on the storytelling role and the SE choice.

---

## Canonical exhibit order (house style of a top applied-micro paper)

1. **Table 1 — Summary / descriptive statistics**
2. **Balance / pre-treatment** (RCT: balance table → `balance-table`; quasi-exp: determinants of treatment timing)
3. **Figure 1 — Motivating descriptive** (raw trends, map, distribution)
4. **Figure 2 — Identification figure** (the centerpiece; differs by method, see below)
5. **Main results table** (progressive columns: naive → +controls → +FE → preferred; report DV mean, FE Yes/No → `reg-table`)
6. **Robustness table** (alt. estimators / samples / **clustering & SE definitions** / placebos)
7. **Heterogeneity** (subgroups/interactions; prefer a **forest plot** over a wide interaction table)
8. **Sensitivity figure** (HonestDiD / RD bandwidth / weak-IV AR confidence sets)
9. **Mechanisms / secondary outcomes**
10. **Appendix** (variable definitions, extra robustness, permutation/randomization inference)

### Identification figure (#4) by methodology
| Methodology | Identification figure | Key robustness |
|---|---|---|
| **DiD / event study** | Event-study plot (leads≈0, lags = dynamic effect; mark ref period −1, shade post) | Goodman-Bacon, CS/SA/dCDH estimators, HonestDiD |
| **RCT** | Covariate-balance plot / outcome distribution by arm | Attrition (Lee bounds), multiple-testing (Romano-Wolf/q-values), RI |
| **IV** | First-stage scatter with fitted line + first-stage F | Effective F (MOP), reduced form, over-ID (Hansen J), AR CIs |
| **RDD** | RD plot (binned means + low-order polynomial) | Bandwidth/polynomial grid, density test (McCrary/CJM), placebo cutoffs, CCT robust SE |

---

## Standard-error decision guide (the inference choice is an exhibit choice)

| Situation | Use | Stata | R |
|---|---|---|---|
| No grouped assignment, large N | Robust HC1 | `, robust` | `vcov = "HC1"` |
| Small sample / influential obs | HC2 / HC3 | `, vce(hc3)` | `vcov = "HC3"` |
| Grouped/serially-correlated assignment | Cluster CR1 at assignment level | `, cluster(g)` | `cluster = ~g` |
| **Few clusters (G ≲ 40)** | **Wild cluster bootstrap** | `boottest x, reps(9999) cluster(g)` | `fwildclusterboot::boottest()` |
| Known/randomized design | Randomization inference | `ritest`/`randcmd` | `ri2` |
| Two-way / spatial / panel-time | Multiway / Conley / Driscoll-Kraay | `vce(cluster g1 g2)`, `acreg`, `xtscc` | `cluster = ~g1+g2`, `vcov = "DK"` |

Report the **same point estimate** under the chosen definition; if few clusters, the wild-bootstrap CI is the honest one (naive cluster SE over-rejects). The live lab demonstrates this interactively.

---

## Paper vs. presentation variants (always produce both when asked)

| Dimension | Paper | Presentation |
|---|---|---|
| Density | All columns/rows, footnotes | Drop R²/extra columns; one message |
| Color | Monochrome, print-safe | Colorblind-safe (Okabe-Ito / `beamer-slides` palette: blue `0072B2`, green `009E73`, orange `D55E00`) |
| Type size | Small | Large, legible from the back |
| Emphasis | Stars + SE + CI | Highlight the one coefficient; build incrementally |
| FE rows | Yes/No kept | Collapse to a note |

---

## Template — DiD event-study figure (the centerpiece)

**Stata**
```stata
* Two-way FE event study, reference period -1, clustered SE
reghdfe y ib(-1).evt, absorb(unit period) cluster(group)
* paper figure (monochrome)
coefplot, keep(*.evt) vertical yline(0) xline(`=`reftick'') ///
    ciopts(recast(rcap)) msymbol(O) mcolor(black) ///
    xtitle("Event time (periods rel. to adoption; -1 = reference)") ///
    ytitle("Effect on y") graphregion(color(white))
graph export "output/figures/event_study.pdf", replace

* Heterogeneity-robust alternatives (staggered adoption):
csdid y, ivar(unit) time(period) gvar(first_treat)      // Callaway-Sant'Anna
eventstudyinteract y leads* lags*, ...                    // Sun-Abraham
* Few clusters? wild bootstrap the post-period ATT:
reghdfe y D, absorb(unit period) cluster(group)
boottest D, reps(9999) cluster(group)
```

**R**
```r
library(fixest); library(ggiplot)
m <- feols(y ~ i(evt, ref = -1) | unit + period, data = df, cluster = ~group)
ggiplot(m, ref.line = -1) +                         # paper: theme_minimal()
  ggplot2::labs(x = "Event time (-1 = reference)", y = "Effect on y")
# presentation variant: larger base_size, color = "#0072B2", highlight post window
# staggered-robust: did::att_gt() (Callaway-Sant'Anna); HonestDiD::createSensitivityResults()
```

---

## Template — Main results & SE robustness

Use `reg-table` Template 2 (progressive specifications) for the main table. For the **SE-robustness** exhibit, hold the spec fixed and vary inference:

```stata
reghdfe y D, absorb(unit period) cluster(group)
estadd local se "Cluster (CR1)"
eststo cr1
* HC, wild bootstrap, RI rows:
reg y D i.unit i.period, robust ; eststo hc1
reghdfe y D, absorb(unit period) cluster(group) ; boottest D, reps(9999) cluster(group)
ritest D _b[D], reps(2000): reghdfe y D, absorb(unit period)
* Tabulate estimate / SE / 95% CI / p across rows -> table_se_robustness.tex
```

---

## Template — RD plot & robustness grid

```stata
rdplot y run, c(0) p(1) graph_options(graphregion(color(white)) ///
    xtitle("Running variable") ytitle("Outcome"))
graph export "output/figures/rd_plot.pdf", replace
* robustness across bandwidths / polynomial + manipulation + placebo:
foreach h in 0.5 0.3 0.2 { rdrobust y run, c(0) h(`h') ; eststo h_`h' }
rdrobust y run, c(0) p(2)                 // quadratic
rddensity run, c(0)                        // manipulation (Cattaneo-Jansson-Ma)
rdrobust cov run, c(0)                     // covariate placebo (jump ≈ 0)
```

```r
library(rdrobust)
rdplot(df$y, df$run, c = 0, p = 1)
rdrobust(df$y, df$run, c = 0)              # CCT robust bias-corrected
rddensity::rddensity(df$run, c = 0)        # manipulation test
```

---

## Template — IV table (OLS vs 2SLS + first stage)

Use `reg-table` Template 5 (IV/LATE). Add the first-stage scatter figure and report the **effective F**:

```stata
ivreg2 y (d = z), robust first             // first-stage F, 2SLS
weakivtest                                  // Montiel-Olea-Pflueger effective F
twoway (scatter d z, mcolor(%30)) (lfit d z), ///
    legend(off) xtitle("Instrument z") ytitle("Endogenous d") ///
    graphregion(color(white))
graph export "output/figures/first_stage.pdf", replace
```

---

## Workflow

1. Confirm methodology + variables + clustering + target(s) + language.
2. Lay out the exhibit checklist (the canonical order above), marking which apply.
3. For each exhibit: call the relevant existing skill (`reg-table`, `balance-table`, `attrition-check`) or use the templates here for the figures those skills don't cover.
4. Stamp each output file with a one-line note: storytelling role + SE choice.
5. If presentation requested, regenerate figures/tables with the `beamer-slides` palette and one-message framing.
6. Point the user to the live lab to explore SE/paper-vs-presentation trade-offs before finalizing.

### File naming
`table01_summary.tex`, `fig02_event_study.pdf`, `table03_main.tex`, `table04_se_robustness.tex`, `fig03_heterogeneity.pdf`, `table05_iv.tex`, `fig02_rd_plot.pdf` — zero-padded, ordered, descriptive.

### Required packages
```stata
ssc install estout reghdfe coefplot boottest rdrobust rddensity ftools ivreg2 ranktest weakivtest
* staggered DiD: ssc install csdid drdid eventstudyinteract
```
```r
install.packages(c("fixest","modelsummary","ggiplot","rdrobust","rddensity","fwildclusterboot","did","HonestDiD"))
```
