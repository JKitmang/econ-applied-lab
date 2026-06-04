# Table & Figure Lab

An interactive prototype for **diagnosing and replicating the canonical tables &
figures of applied-economics papers** — across methodologies, with a live
standard-error lab, a paper-vs-presentation toggle, and one-click LaTeX export.

Built to be **simple and reproducible**: a single static page, **no build step,
no dependencies, no real/proprietary data**. Every exhibit is drawn from a
simulated data-generating process with a fixed seed, so the same inputs always
produce the same numbers. Replication packages from journals (indexed by
[EJD](https://ejd.econ.mathematik.uni-ulm.de/)) are used as *visual
ground-truth* — the page links out to them per exhibit.

---

## Open it

**Option A — just open the file.** Double-click `index.html`. Because there are
no module imports or `fetch` calls, it runs directly from `file://` (e.g. from
Google Drive).

**Option B — local server** (nicer URLs, console available):

```bash
cd "Table-Figure Lab"
python3 -m http.server 8732
# → http://localhost:8732
```

In Claude Code, the launch config `tflab` in `.claude/launch.json` starts this
for the preview tools.

---

## What's inside

- **Methodology switcher** — DiD / event study (full depth), plus RCT, IV, RDD
  (identification figure + the strategy's exhibit grammar).
- **Exhibit cards** following the canonical paper order (see
  [`docs/applied-econ-table-taxonomy.md`](docs/applied-econ-table-taxonomy.md)):
  summary stats → raw trends → event study → main results → SE robustness →
  heterogeneity.
- **Live standard-error lab** — switch between iid, HC1/2/3, cluster-robust,
  **wild cluster bootstrap**, and **randomization inference**; SEs, stars, and
  CIs recompute everywhere. With few clusters you can watch the wild bootstrap CI
  widen and the stars drop — the whole point.
- **Simulation controls** — clusters G, units/cluster, periods T, treatment
  effect, ICC, pre-trend violation, adoption timing, bootstrap reps, and seed.
- **Paper ↔ Presentation** themes — booktabs/monochrome vs. the colorblind-safe
  Beamer palette with the key coefficient highlighted.
- **Replication recipe** per card — the matching Stata (`reghdfe`/`esttab`/
  `boottest`) and R (`feols`/`modelsummary`/`boottest`) snippets, plus
  **Copy-LaTeX** export in booktabs/`threeparttable` form.

---

## Architecture

```
Table-Figure Lab/
  index.html          # shell, controls, storytelling panel
  css/styles.css      # paper vs presentation themes; booktabs table styling
  js/econ.js          # linear algebra, OLS, within-FE, SE estimators, wild bootstrap, seeded RNG
  js/dgp.js           # data-generating processes (DiD staggered, RCT, IV, RDD)
  js/models.js        # estimation pipelines: event study, progressive FE columns, SE dispatch
  js/render.js        # table renderers (HTML + LaTeX export) and hand-rolled SVG figures
  js/taxonomy.js      # the diagnosed exhibit catalog + storytelling notes + EJD links
  js/app.js           # state, wiring, per-card rendering
  docs/applied-econ-table-taxonomy.md   # the diagnosis (Part 1 deliverable)
```

Everything is plain ES5-ish JavaScript attached to globals (`Econ`, `DGP`,
`Models`, `Render`, `TAXONOMY`) and loaded with `<script>` tags — deliberately,
so it works from `file://`. The econometrics (OLS, sandwich/cluster variance,
wild cluster bootstrap, randomization inference, two-way within transformation)
is implemented from scratch in `js/econ.js` so every number is auditable.

> Note: the data catalog lives in `js/taxonomy.js` (a global) rather than a
> `data/*.json` file so the page can open offline without `fetch`/CORS issues.

---

## How it maps to your existing skills

This prototype is the **groundwork for a Phase-2 Claude Code skill** that emits
tables/figures matching your installed toolkit. The conventions already line up:

| Lab feature | Your skill | Convention matched |
|-------------|------------|--------------------|
| Main results table | `reg-table` | booktabs 3-line table, SE in parens, stars `* .10 / ** .05 / *** .01`, DV mean, FE as Yes/No, `threeparttable` |
| Balance table (RCT) | `balance-table` | control/treatment means, normalized difference, |ND| > 0.25 flag |
| Presentation theme | `beamer-slides` | Okabe–Ito palette (blue #0072B2, green #009E73, orange #D55E00, yellow #F0E442), no vertical rules |
| SE lab | `power-calc` / PAP | cluster level = assignment level; few-cluster wild bootstrap |
| Copy-LaTeX export | `reg-table` | drop-in `\begin{table} … threeparttable … booktabs` fragment |

---

## Extending it

- **Add an exhibit:** append an entry to the methodology's `exhibits` array in
  `js/taxonomy.js` (give it a `render` key), implement the model in
  `js/models.js` and the renderer in `js/render.js`, then dispatch it in
  `renderExhibit()` in `js/app.js`.
- **Deepen a stub (RCT/IV/RDD):** the DGPs already exist in `js/dgp.js`; add the
  remaining exhibit cards the same way DiD's six are built.
- **Add an SE definition:** add a branch in `seForTarget()` (`js/models.js`) and
  an `<option>` in `index.html`.

---

## Roadmap (Phase 2)

1. Package the generator as a **Claude Code skill** (`table-figure-lab`) that
   takes a methodology + a dataset/spec and emits the Stata/R script **and** the
   LaTeX/figure, reusing `reg-table` and `beamer-slides`.
2. Flesh out RCT/IV/RDD to the same depth as DiD.
3. Optional **hybrid data mode:** load a real public replication dataset
   (AEA/openICPSR/Dataverse, via EJD) to reproduce a specific paper's exhibit.
