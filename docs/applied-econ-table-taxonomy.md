# Diagnosis: Tables & Figures in Top Applied-Economics Journals

*A catalog of the canonical exhibit types, their ordering, and how they vary by
identification strategy — the empirical grounding for the Table & Figure Lab.*

**Scope.** AER, AER: Insights, AEA P&P, the four AEJ titles (Applied, Policy,
Micro, Macro), REStud, REStat, QJE, JPE, plus the econometrics anchors.
**Ground-truth source.** [EJD — Find Economic Articles with Data](https://ejd.econ.mathematik.uni-ulm.de/)
(S. Kranz, Uni Ulm), a searchable index of articles with replication packages,
queryable by URL and by methodological keyword (`DID`, `IV`, `RCT`, `RDD`, …).
Links below are EJD searches you can re-run to pull current exemplars.

---

## 1. The canonical "house order" of exhibits

Across applied-micro papers there is a remarkably stable narrative sequence. Not
every paper has every item, and appendices carry the overflow, but the spine is:

| # | Exhibit | Purpose / what the reader is checking |
|---|---------|----------------------------------------|
| 1 | **Table 1 — Summary / descriptive statistics** | Sample size, key variables, group composition. Orientation. |
| 2 | **Balance / pre-treatment comparison** | RCT: covariate balance. Quasi-exp: correlates of treatment/timing. |
| 3 | **Figure 1 — Motivating descriptive figure** | Raw trends, a map, or a distribution. Earns trust *before* any modeling. |
| 4 | **Figure 2 — Identification figure** | Event study (DiD), first-stage scatter (IV), RD plot (RDD), balance CDF (RCT). The visual heart of the design. |
| 5 | **Main results table** | Progressive columns: bivariate → +controls → +FE → preferred. DV mean, FE rows, N, R². |
| 6 | **Robustness table** | Alternative estimators, samples, clustering/SE definitions, placebos. |
| 7 | **Heterogeneity** | Subgroups / interactions. Increasingly a **forest/coefficient plot** rather than a wide table. |
| 8 | **Sensitivity figure** | HonestDiD (Rambachan–Roth), RD bandwidth curves, weak-IV Anderson–Rubin sets. |
| 9 | **Mechanisms / secondary outcomes** | Why the effect happens; intermediate outcomes. |
| 10 | **Appendix** | Variable definitions, extra robustness, permutation/randomization inference, full coefficient tables. |

**Trend (2018→):** the center of gravity has moved from tables to **figures** —
especially the identification figure and coefficient/forest plots — and from a
single point estimate to **a curve of estimates** (dynamic effects, sensitivity
ranges). The AEA Data Editor's mandatory replication package has also made
*"every exhibit ⇄ one script with a fixed seed"* a de facto formatting rule.

---

## 2. Variation by identification strategy

### 2.1 Difference-in-Differences / Event Study  ·  [EJD: difference-in-differences](https://ejd.econ.mathematik.uni-ulm.de/?search=difference-in-differences&sortby=date)

The most active methodological frontier of the last decade. The naive two-way
fixed-effects (TWFE) estimator is biased under staggered adoption with
heterogeneous effects (negative weighting), which reshaped the exhibit set.

- **Identification figure = event-study plot.** Leads (≈ 0) defend parallel
  trends; lags trace the dynamic ATT. Reference period (usually −1) omitted and
  marked; CIs always shown; post window shaded.
- **Main table** reports a single post ATT across progressive FE columns.
- **Estimators now expected** alongside TWFE: Callaway–Sant'Anna, Sun–Abraham,
  de Chaisemartin–D'Haultfœuille, Borusyak–Jaravel–Spiess.
- **Diagnostics:** Goodman-Bacon decomposition (what TWFE is averaging),
  explicit pre-trend tests, and **HonestDiD** (Rambachan–Roth) sensitivity.
- **Inference:** cluster at the treatment level; with **few clusters** (≲ 40)
  the **wild cluster bootstrap** (Cameron–Gelbach–Miller; `boottest`) is
  expected because CR1 over-rejects.
- EJD: [event study](https://ejd.econ.mathematik.uni-ulm.de/?search=event%20study&sortby=date) ·
  [parallel trends](https://ejd.econ.mathematik.uni-ulm.de/?search=parallel%20trends&sortby=date) ·
  [Sun Abraham](https://ejd.econ.mathematik.uni-ulm.de/?search=Sun%20Abraham&sortby=date)

### 2.2 Randomized Controlled Trial  ·  [EJD: RCT](https://ejd.econ.mathematik.uni-ulm.de/?search=randomized%20controlled%20trial&sortby=date)

- **Balance table is mandatory** and usually Table 1: control vs treatment
  means, normalized differences (flag |ND| > 0.25, Imbens–Wooldridge), a joint
  orthogonality F-test.
- **Attrition table** when follow-up is incomplete (differential attrition;
  Lee bounds).
- **Main estimates:** ITT (OLS/ANCOVA) and **LATE/2SLS** under non-compliance.
- **Pre-specified vs exploratory** clearly separated; **multiple-testing**
  corrections (Romano–Wolf, sharpened FDR q-values) for outcome families.
- **Inference:** robust or clustered at the unit of randomization;
  **randomization inference** (Fisher) for design-based exact p-values.

### 2.3 Instrumental Variables  ·  [EJD: IV](https://ejd.econ.mathematik.uni-ulm.de/?search=instrumental%20variables&sortby=date)

A triptych of exhibits:
- **First stage** (often a figure + a table) with an **F / effective F**
  statistic (Montiel-Olea–Pflueger). A weak first stage is shown, not hidden.
- **Reduced form** (the instrument's effect on the outcome directly).
- **2SLS** as the headline.
- **Over-identification:** Hansen J. **Weak-IV-robust inference:**
  Anderson–Rubin confidence sets when relevance is marginal.

### 2.4 Regression Discontinuity  ·  [EJD: RDD](https://ejd.econ.mathematik.uni-ulm.de/?search=regression%20discontinuity&sortby=date)

- **The RD plot is the paper:** binned means + a low-order polynomial on each
  side. Bin count is a real choice (too few hides noise; too many hides signal).
- **Robustness:** across **bandwidths** and **polynomial orders** (`rdrobust`,
  Calonico–Cattaneo–Titiunik robust bias-corrected SEs).
- **Validity:** density/manipulation test (McCrary; Cattaneo–Jansson–Ma),
  covariate continuity at the cutoff, **placebo cutoffs**.

---

## 3. Standard-error taxonomy (the simulation lab)

Inference is where applied papers most often go wrong, so the Lab makes the SE
choice an interactive knob. The point estimate is held fixed; only the inference
moves.

| Definition | Use it when | Notes |
|------------|-------------|-------|
| **Classical (iid)** | Almost never | Homoskedastic baseline; shown for contrast. |
| **Robust HC0–HC3** | No grouped assignment, large N | HC1 = Stata `, robust`; HC3 best for small samples / high leverage. |
| **Cluster-robust (CR1)** | Grouped assignment (DiD/RCT) | Cluster **at the level of treatment assignment**. The applied default. |
| **Wild cluster bootstrap** | **Few clusters** (G ≲ 40) | Cameron–Gelbach–Miller; `boottest`. CR1 over-rejects here. |
| **CR2 / Bell–McCaffrey** | Few clusters, analytic alt. | `clubSandwich`; degrees-of-freedom correction. |
| **Randomization inference** | Known/randomized assignment | Fisher exact p-values; design-based. |
| **Two-way / multiway** | Two clustering dimensions | e.g. firm × year. |
| **Conley (spatial HAC)** | Spatial correlation | Distance-decayed. |
| **Driscoll–Kraay** | Panel with cross-sectional dependence | Time-series-robust panel SEs. |

The Lab implements iid, HC0–HC3, CR1, the wild cluster bootstrap
(percentile-*t*), and randomization inference from first principles; the rest
are documented here as the natural extension set.

---

## 4. Communication & data-storytelling principles

Distilled from the AEA Data Editor guidance, Schwabish's *Better Data
Visualizations*, and the conventions visible across the journals above.

1. **One message per exhibit.** If you can't state the takeaway in a sentence,
   split it.
2. **Figure > table** when the point is a relationship or a trend; **table >
   figure** when exact magnitudes matter.
3. **The identification figure is the emotional core** — make it the cleanest
   object on the page.
4. **Tables:** booktabs (no vertical rules), align on the decimal, 2–3
   significant figures, units in the header, notes in `threeparttable`, report
   the **DV mean** and **FE rows** as Yes/No, stars *and* SEs/CIs.
5. **Color:** colorblind-safe (Okabe–Ito); one consistent color for the
   treatment throughout the paper.
6. **Replicability is a communication property:** every exhibit maps to one
   script with a fixed seed.

### Paper vs. presentation

| Dimension | Paper | Presentation |
|-----------|-------|--------------|
| Density | Complete: all columns/rows, footnotes | Sparse: drop R² and extra columns; one message |
| Color | Monochrome, print-safe | Colorblind-safe accents (Beamer palette) |
| Type size | Small; reader controls zoom | Large; legible from the back row |
| Emphasis | Stars + SE + CI; reader judges | Highlight the one coefficient; build incrementally |
| FE / controls | Yes/No rows kept | Often collapsed to a note |

---

## 5. How this maps to the prototype

Each row of §1 is a **card** in `index.html`; each strategy in §2 is the
**methodology switcher**; §3 is the **SE selector**; §4 is the **paper /
presentation toggle** and the per-card storytelling callouts. Exhibits are drawn
from simulated DGPs (no proprietary data) so the whole thing is self-contained
and reproducible — EJD links point to the real papers as visual ground-truth.
