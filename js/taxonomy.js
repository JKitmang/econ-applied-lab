/* ============================================================================
   taxonomy.js — The diagnosed catalog of canonical exhibits per methodology,
   their position in a paper, storytelling notes, and EJD reference searches.
   Loaded as a <script> (global TAXONOMY) so the page also opens from file://.
   ========================================================================== */
(function (global) {
  "use strict";
  const EJD = "https://ejd.econ.mathematik.uni-ulm.de/?search=";
  const ejd = (q) => EJD + encodeURIComponent(q) + "&sortby=date";

  const TAXONOMY = {
    // ---- canonical "house order" of exhibits in an applied-micro paper -----
    canonicalOrder: [
      "Table 1 — Summary / descriptive statistics",
      "Balance / pre-treatment comparison",
      "Figure 1 — Motivating descriptive figure (raw trends, map, distribution)",
      "Figure 2 — Identification figure (event study / first stage / RD plot)",
      "Main results table (progressive columns: +controls, +FE, preferred)",
      "Robustness table (alt. estimators, samples, clustering, placebos)",
      "Heterogeneity (subgroups / interactions — forest plot)",
      "Sensitivity figure (HonestDiD / bandwidth / AR confidence sets)",
      "Mechanisms / secondary outcomes",
      "Appendix (variable definitions, extra robustness, permutation inference)",
    ],

    seGlossary: [
      { key: "iid", label: "Classical (iid)", when: "Homoskedasticity assumed. Almost never defensible in applied micro; shown as a baseline." },
      { key: "HC1", label: "Robust HC1", when: "Stata's `, robust`. Heteroskedasticity-robust; fine with no grouped assignment and large N." },
      { key: "HC2", label: "Robust HC2", when: "Leverage-adjusted; better small-sample behavior than HC1." },
      { key: "HC3", label: "Robust HC3", when: "Most conservative HC; preferred for small samples / influential points." },
      { key: "cluster", label: "Cluster-robust (CR1)", when: "Cluster at the level of treatment assignment. The default in DiD/RCT with grouped designs." },
      { key: "wildboot", label: "Wild cluster bootstrap", when: "Few clusters (G ≲ 40). Cameron–Gelbach–Miller; `boottest`. CR1 over-rejects here." },
      { key: "ri", label: "Randomization inference", when: "Design-based exact p-values (Fisher). Natural when assignment is known/randomized." },
    ],

    storytelling: {
      general: [
        "One message per exhibit — if you can't state the takeaway in a sentence, split it.",
        "Prefer a figure over a table when the point is a relationship or a trend.",
        "The identification figure (event study / first stage / RD plot) is the emotional core — make it the cleanest object.",
        "Round to 2–3 significant figures; ragged precision signals noise as if it were signal.",
        "Every exhibit maps to one script with a fixed seed — replicability is a communication property, not just an ethics one.",
      ],
      paperVsPresent: [
        ["Density", "Complete: all columns, all rows, footnotes.", "Sparse: drop R², extra columns; one message."],
        ["Color", "Monochrome, print-safe.", "Colorblind-safe accents (Okabe–Ito / your Beamer palette)."],
        ["Type size", "Small; reader controls zoom.", "Large; legible from the back row."],
        ["Emphasis", "Stars + SEs + CI; let the reader judge.", "Highlight the one coefficient; build incrementally."],
        ["FE / controls", "Yes/No rows kept.", "Often collapsed to a note."],
      ],
    },

    methodologies: {
      // ================= DiD / event study (lead, full depth) =============
      did: {
        name: "Difference-in-Differences / Event Study",
        blurb: "The workhorse of modern applied micro. Staggered adoption broke naive TWFE; the field moved to event-study plots plus heterogeneity-robust estimators (Callaway–Sant'Anna, Sun–Abraham, de Chaisemartin–D'Haultfœuille, Borusyak–Jaravel–Spiess).",
        estimators: ["TWFE", "Callaway–Sant'Anna", "Sun–Abraham", "de Chaisemartin–D'Haultfœuille", "Borusyak–Jaravel–Spiess"],
        ejdSearch: ejd("difference-in-differences"),
        diagnostics: ["Goodman-Bacon decomposition", "Pre-trend / parallel-trends test", "HonestDiD (Rambachan–Roth) sensitivity"],
        exhibits: [
          { id: "summary", n: "Table 1", title: "Summary statistics", type: "table", render: "summary",
            story: "Orient the reader: sample size, key variables, and how ever-treated vs never-treated units compare before any modeling.",
            ejd: [{ label: "EJD: DiD papers", url: ejd("difference-in-differences") }] },
          { id: "trends", n: "Figure 1", title: "Motivating raw trends", type: "figure", render: "trends",
            story: "Show the unadjusted series for treated vs control. This is the 'before the magic' figure — it earns the reader's trust that parallel trends is plausible.",
            ejd: [{ label: "EJD: parallel trends", url: ejd("parallel trends") }] },
          { id: "eventstudy", n: "Figure 2", title: "Event-study plot", type: "figure", render: "eventstudy",
            story: "The centerpiece. Leads ≈ 0 support parallel trends; lags trace the dynamic effect. Always show CIs, mark the reference period (−1), and shade the post window.",
            ejd: [{ label: "EJD: event study", url: ejd("event study") }, { label: "EJD: Sun Abraham", url: ejd("Sun Abraham") }] },
          { id: "main", n: "Table 2", title: "Main results", type: "table", render: "main",
            story: "Build the spec in columns: naive → +unit FE → +time FE (preferred) → +covariate. Report the DV mean and FE rows so the reader can locate identifying variation.",
            ejd: [{ label: "EJD: two-way fixed effects", url: ejd("two-way fixed effects") }] },
          { id: "robustness", n: "Table 3", title: "Inference under different SEs", type: "table", render: "robustness",
            story: "The point estimate is fixed; only inference moves. With few clusters, the wild cluster bootstrap CI is the honest one — naive cluster SEs over-reject.",
            ejd: [{ label: "EJD: wild bootstrap", url: ejd("wild bootstrap") }] },
          { id: "heterogeneity", n: "Figure 3", title: "Heterogeneity (forest plot)", type: "figure", render: "heterogeneity",
            story: "A forest plot beats a wide interaction table: one row per subgroup, point + CI, zero line. Reader sees the pattern in one glance.",
            ejd: [{ label: "EJD: heterogeneous effects", url: ejd("heterogeneous treatment effects") }] },
        ],
      },

      // ================= RCT (stub: taxonomy + identification) ============
      rct: {
        name: "Randomized Controlled Trial",
        blurb: "Cleanest identification. The exhibit grammar centers on a balance table, ITT (and LATE under non-compliance), pre-specified heterogeneity, attrition, and multiple-testing corrections (Romano–Wolf, sharpened q-values).",
        estimators: ["ITT (OLS)", "LATE / 2SLS", "ANCOVA", "Lee bounds (attrition)"],
        ejdSearch: ejd("randomized controlled trial"),
        diagnostics: ["Balance / joint orthogonality F-test", "Attrition analysis", "Randomization inference", "Multiple-testing adjustment"],
        exhibits: [
          { id: "balance", n: "Table 1", title: "Balance table", type: "table", render: "balance",
            story: "Mandatory first exhibit. Control vs treatment means, normalized differences (flag |ND|>0.25), and a joint test. The reader is checking that randomization 'took'.",
            ejd: [{ label: "EJD: balance table / RCT", url: ejd("balance randomized") }] },
        ],
      },

      // ================= IV (stub: taxonomy + first stage) ================
      iv: {
        name: "Instrumental Variables",
        blurb: "The exhibit set is a triptych: first stage (with an F / effective-F statistic), reduced form, and 2SLS. Weak instruments demand identification-robust inference (Anderson–Rubin); over-ID needs a Hansen J.",
        estimators: ["2SLS", "LIML", "Anderson–Rubin CIs", "Montiel-Olea–Pflueger effective F"],
        ejdSearch: ejd("instrumental variables"),
        diagnostics: ["First-stage F / effective F", "Over-identification (Hansen J)", "Weak-IV-robust confidence sets"],
        exhibits: [
          { id: "firststage", n: "Figure 2", title: "First-stage scatter", type: "figure", render: "firststage",
            story: "Make the instrument's relevance visible: scatter of endogenous regressor on instrument with the fitted line and the first-stage F. A weak first stage is a fatal flaw shown honestly.",
            ejd: [{ label: "EJD: first stage / IV", url: ejd("first stage instrument") }] },
        ],
      },

      // ================= RDD (stub: taxonomy + RD plot) ===================
      rdd: {
        name: "Regression Discontinuity",
        blurb: "The RD plot is the paper. Binned means + a low-order polynomial each side, then robustness across bandwidths and polynomial orders, a manipulation/density test, and placebo cutoffs. Inference uses robust bias-corrected SEs (CCT).",
        estimators: ["Local linear (CCT)", "rdrobust", "Density test (Cattaneo–Jansson–Ma)"],
        ejdSearch: ejd("regression discontinuity"),
        diagnostics: ["Bandwidth / polynomial robustness", "McCrary / CJM density test", "Placebo cutoffs & covariate continuity"],
        exhibits: [
          { id: "rdplot", n: "Figure 2", title: "RD plot", type: "figure", render: "rdplot",
            story: "Binned scatter so the eye sees the discontinuity, with the fitted jump annotated. Don't over-smooth: too few bins hides noise, too many hides the signal.",
            ejd: [{ label: "EJD: RD design", url: ejd("regression discontinuity") }] },
        ],
      },
    },
  };

  global.TAXONOMY = TAXONOMY;
})(typeof window !== "undefined" ? window : this);
