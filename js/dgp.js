/* ============================================================================
   dgp.js — Data-generating processes per methodology (pure JS, seedable)
   ----------------------------------------------------------------------------
   Each generator returns a flat, column-oriented dataset (arrays aligned by
   observation) plus metadata the models/renderers need. No real data is used:
   self-contained and 100% reproducible from (params, seed).
   Exposed as global `DGP`.
   ========================================================================== */
(function (global) {
  "use strict";
  const RNG = global.Econ.makeRNG;

  /* ---------------------------------------------------------------------- *
   *  Staggered Difference-in-Differences panel
   * ---------------------------------------------------------------------- */
  function did(p) {
    p = Object.assign({
      nClusters: 6,        // G groups (e.g. states) — clustering level
      unitsPerCluster: 8,  // units within a group
      T: 12,               // periods
      tau: 1.0,            // treatment effect (long-run)
      ramp: 3,             // periods to reach full effect (1 = instant)
      preTrend: 0.0,       // pre-trend slope among treated (violation knob)
      icc: 0.2,            // intra-cluster correlation knob (group random effect)
      sigma: 1.0,          // idiosyncratic noise sd scale
      timing: "staggered", // "staggered" | "single"
      shareTreated: 0.66,  // fraction of clusters ever treated
      seed: 12345,
    }, p || {});

    const rng = RNG(p.seed);
    const G = p.nClusters, T = p.T;
    const nTreatedG = Math.max(1, Math.round(G * p.shareTreated));

    // --- adoption period per cluster (Infinity = never-treated control) ---
    const adopt = [];
    const firstAdopt = Math.max(2, Math.round(T * 0.3));
    const lastAdopt = Math.max(firstAdopt + 1, Math.round(T * 0.75));
    for (let g = 0; g < G; g++) {
      if (g < nTreatedG) {
        if (p.timing === "single") adopt[g] = Math.round((firstAdopt + lastAdopt) / 2);
        else adopt[g] = firstAdopt + Math.round(((lastAdopt - firstAdopt) * g) / Math.max(nTreatedG - 1, 1));
      } else adopt[g] = Infinity;
    }

    const groupSd = Math.sqrt(p.icc) * 1.5;
    const unitSd = 1.0;
    const epsSd = Math.sqrt(1 - Math.min(p.icc, 0.95)) * p.sigma;
    const timeFE = Array.from({ length: T }, (_, t) => 0.15 * t + 0.4 * rng.normal()); // common trend + shock

    const effect = (e) => (e < 0 ? 0 : p.tau * Math.min((e + 1) / Math.max(p.ramp, 1), 1));

    const unit = [], time = [], group = [], y = [], x = [], D = [], evt = [], everTreated = [], adoptCol = [];
    let uid = 0;
    for (let g = 0; g < G; g++) {
      const gEff = groupSd * rng.normal();
      for (let u = 0; u < p.unitsPerCluster; u++) {
        const aEff = unitSd * rng.normal();
        const xUnit = 0.5 * rng.normal(); // a (mostly time-invariant) covariate
        for (let t = 0; t < T; t++) {
          const e = adopt[g] === Infinity ? -Infinity : t - adopt[g];
          const post = e >= 0 ? 1 : 0;
          const pre = adopt[g] !== Infinity && e < 0 ? e : 0; // negative event time among treated
          const yi = aEff + gEff + timeFE[t] + effect(e) + p.preTrend * pre + 0.3 * xUnit + epsSd * rng.normal();
          unit.push(uid); time.push(t); group.push(g);
          y.push(yi); x.push(xUnit + 0.1 * rng.normal());
          D.push(post); everTreated.push(adopt[g] === Infinity ? 0 : 1);
          adoptCol.push(adopt[g]);
          // event time, binned to [-5, +5] for the event-study plot
          let ev = adopt[g] === Infinity ? null : Math.max(-6, Math.min(6, t - adopt[g]));
          evt.push(ev);
        }
        uid++;
      }
    }
    return {
      kind: "did",
      n: unit.length, G, T, adopt,
      cols: { unit, time, group, y, x, D, evt, everTreated, adopt: adoptCol },
      params: p,
      eventTimes: range(-5, 5), // displayed leads/lags (-1 omitted as reference)
    };
  }

  /* ---------------------------------------------------------------------- *
   *  RCT (single cross-section, possibly clustered assignment)
   * ---------------------------------------------------------------------- */
  function rct(p) {
    p = Object.assign({
      n: 1000, tau: 0.3, nClusters: 40, icc: 0.05, sigma: 1.0, seed: 777, controls: 3,
      nStrata: 4, complyRate: 0.7, attritBase: 0.10, attritDiff: 0.04,
    }, p || {});
    const rng = RNG(p.seed);
    const perCl = Math.max(1, Math.round(p.n / p.nClusters));
    const N = perCl * p.nClusters;
    const unit = [], cluster = [], strata = [], W = [], D = [], y = [], X = [], attrit = [], complier = [];
    // stratified assignment: strata = g % nStrata; treatment varies WITHIN stratum
    // (alternate by the cluster's index within its stratum) so W is not collinear with strata.
    const stratumOf = (g) => g % p.nStrata;
    const clTreat = Array.from({ length: p.nClusters }, (_, g) => (Math.floor(g / p.nStrata) % 2 === 0 ? 1 : 0));
    const groupSd = Math.sqrt(p.icc) * 1.5, epsSd = Math.sqrt(1 - Math.min(p.icc, 0.95)) * p.sigma;
    let id = 0;
    for (let g = 0; g < p.nClusters; g++) {
      const gEff = groupSd * rng.normal();
      const st = stratumOf(g);
      for (let u = 0; u < perCl; u++) {
        const xs = Array.from({ length: p.controls }, () => rng.normal());
        const isComplier = rng.uniform() < p.complyRate;
        const assigned = clTreat[g];
        const received = assigned * (isComplier ? 1 : 0);      // one-sided noncompliance
        const yi = 1 + p.tau * received + 0.4 * xs[0] + 0.2 * st + gEff + epsSd * rng.normal();
        const pAtt = p.attritBase + p.attritDiff * assigned;   // differential attrition
        const att = rng.uniform() < pAtt ? 1 : 0;
        unit.push(id++); cluster.push(g); strata.push(st);
        W.push(assigned); D.push(received); y.push(yi); X.push(xs);
        attrit.push(att); complier.push(isComplier ? 1 : 0);
      }
    }
    return { kind: "rct", n: N, cols: { unit, cluster, strata, W, D, y, X, attrit, complier }, clTreat, params: p };
  }

  /* ---------------------------------------------------------------------- *
   *  IV (single endogenous regressor, one instrument)
   * ---------------------------------------------------------------------- */
  function iv(p) {
    p = Object.assign({ n: 800, beta: 0.5, strength: 0.6, sigma: 1.0, seed: 99 }, p || {});
    const rng = RNG(p.seed);
    const z = [], d = [], y = [];
    for (let i = 0; i < p.n; i++) {
      const zi = rng.normal();
      const u = rng.normal();                         // confounder
      const di = p.strength * zi + 0.6 * u + 0.5 * rng.normal(); // first stage + endogeneity
      const yi = 1 + p.beta * di + 0.8 * u + p.sigma * rng.normal();
      z.push(zi); d.push(di); y.push(yi);
    }
    return { kind: "iv", n: p.n, cols: { z, d, y }, params: p };
  }

  /* ---------------------------------------------------------------------- *
   *  Sharp RDD
   * ---------------------------------------------------------------------- */
  function rdd(p) {
    p = Object.assign({ n: 1500, cutoff: 0, jump: 0.6, slope: 0.8, curve: 0.3, sigma: 0.4, seed: 2024 }, p || {});
    const rng = RNG(p.seed);
    const run = [], y = [], T = [], cov = [];
    for (let i = 0; i < p.n; i++) {
      const r = 2 * (rng.uniform() - 0.5) * 1.0; // running var in [-1,1]
      const above = r >= p.cutoff ? 1 : 0;
      const yi = 0.5 + p.slope * r + p.curve * r * r + p.jump * above + p.sigma * rng.normal();
      const ci = 0.3 + 0.5 * r + 0.4 * rng.normal();  // a covariate, CONTINUOUS at cutoff (placebo)
      run.push(r); y.push(yi); T.push(above); cov.push(ci);
    }
    return { kind: "rdd", n: p.n, cols: { run, y, T, cov }, params: p };
  }

  function range(a, b) { const o = []; for (let i = a; i <= b; i++) o.push(i); return o; }

  global.DGP = { did, rct, iv, rdd };
})(typeof window !== "undefined" ? window : this);
