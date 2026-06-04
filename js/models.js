/* ============================================================================
   models.js — Estimation pipelines that turn a DGP dataset + an SE method into
   the exact objects the renderers draw. Bridges DGP ⇄ Econ ⇄ render.
   Exposed as global `Models`.
   ========================================================================== */
(function (global) {
  "use strict";
  const E = global.Econ;
  const Z = 1.959964; // normal 97.5% quantile for 95% CIs

  /* ---- one-way demean helper ------------------------------------------- */
  function oneWayDemean(v, id) {
    const sum = {}, cnt = {};
    for (let i = 0; i < v.length; i++) { sum[id[i]] = (sum[id[i]] || 0) + v[i]; cnt[id[i]] = (cnt[id[i]] || 0) + 1; }
    return v.map((x, i) => x - sum[id[i]] / cnt[id[i]]);
  }
  function nUnique(arr) { const s = new Set(arr); return s.size; }
  function withinR2(yd, resid) {
    const tss = yd.reduce((s, v) => s + v * v, 0);
    const rss = resid.reduce((s, v) => s + v * v, 0);
    return tss > 0 ? 1 - rss / tss : 0;
  }

  /* ---- SE dispatch for a single target coefficient --------------------- *
     Returns { se, lo, hi, p, stars, label, note } for fit.beta[target]. */
  function seForTarget(opts) {
    const { fit, method, clusterIds, dfAbsorb, target, buildFitY, rng, B } = opts;
    const b = fit.beta[target];
    let V;
    if (method === "iid") V = E.vcovClassical(fit);
    else if (method === "HC1") V = E.vcovHC(fit, "HC1");
    else if (method === "HC2") V = E.vcovHC(fit, "HC2");
    else if (method === "HC3") V = E.vcovHC(fit, "HC3");
    else if (method === "cluster") V = E.vcovCluster(fit, clusterIds, dfAbsorb);
    else if (method === "wildboot") {
      const r = E.wildClusterBootCI(buildFitY, fit, clusterIds, target, B || 199, rng, dfAbsorb);
      return { se: r.impliedSE, lo: r.ci[0], hi: r.ci[1], p: r.pvalue, stars: E.starsFromP(r.pvalue),
               label: METHOD_LABEL.wildboot, note: "percentile-t, " + (B || 199) + " reps" };
    } else V = E.vcovCluster(fit, clusterIds, dfAbsorb);
    const se = Math.sqrt(Math.max(V[target][target], 0));
    const p = E.pValueFromZ(b / se);
    const Gnote = V._G != null ? ("G = " + V._G + " clusters") : "";
    return { se, lo: b - Z * se, hi: b + Z * se, p, stars: E.starsFromP(p), label: METHOD_LABEL[method] || method, note: Gnote };
  }

  const METHOD_LABEL = {
    iid: "Classical (iid)", HC1: "Robust HC1", HC2: "Robust HC2", HC3: "Robust HC3",
    cluster: "Cluster-robust", wildboot: "Wild cluster bootstrap", ri: "Randomization inference",
  };

  /* ====================================================================== *
   *  DiD: event study (two-way FE via demeaning; event-time dummies)
   * ====================================================================== */
  function didEventStudy(data, method, rng, B) {
    const c = data.cols;
    const evTimes = data.eventTimes.filter(e => e !== -1); // -1 omitted (reference)
    const clampE = (e) => (e == null ? null : Math.max(-5, Math.min(5, e)));
    // build raw event-time dummies (treated units only carry 1s; controls all 0)
    const cols = evTimes.map(e => c.evt.map(ev => (clampE(ev) === e ? 1 : 0)));
    // two-way demean y and each dummy
    const yd = E.twoWayDemean(c.y, c.unit, c.time);
    const Xd = cols.map(col => E.twoWayDemean(col, c.unit, c.time));
    const X = []; for (let i = 0; i < data.n; i++) X.push(Xd.map(col => col[i]));
    const fit = E.ols(X, yd);
    const dfAbsorb = nUnique(c.unit) + nUnique(c.time) - 1;
    const buildFitY = (ystar) => E.ols(X, ystar);
    const seMethod = (method === "ri" ? "cluster" : method); // RI not defined per-lead → show cluster band
    const points = evTimes.map((e, j) => {
      const r = seForTarget({ fit, method: seMethod, clusterIds: c.group, dfAbsorb, target: j, buildFitY, rng, B });
      return { e, coef: fit.beta[j], se: r.se, lo: r.lo, hi: r.hi };
    });
    // inject reference period (-1) at coef 0
    points.push({ e: -1, coef: 0, se: 0, lo: 0, hi: 0, ref: true });
    points.sort((a, b2) => a.e - b2.e);
    return { points, method: seMethod, refPeriod: -1, trueEffect: data.params.tau };
  }

  /* ====================================================================== *
   *  DiD: main results — progressive FE columns
   * ====================================================================== */
  function didMain(data, method, rng, B) {
    const c = data.cols;
    const specs = [];

    // (1) No FE: y ~ a + D
    {
      const X = c.D.map(d => [1, d]);
      const fit = E.ols(X, c.y);
      const r = seForTarget({ fit, method: method === "ri" ? "cluster" : method, clusterIds: c.group, dfAbsorb: 0, target: 1, buildFitY: (ys) => E.ols(X, ys), rng, B });
      specs.push(col("(1) OLS", fit.beta[1], r, false, false, false, fit.n, fit.r2));
    }
    // (2) Unit FE
    {
      const yd = oneWayDemean(c.y, c.unit), Dd = oneWayDemean(c.D, c.unit);
      const X = Dd.map(d => [d]); const fit = E.ols(X, yd);
      const r = seForTarget({ fit, method: method === "ri" ? "cluster" : method, clusterIds: c.group, dfAbsorb: nUnique(c.unit), target: 0, buildFitY: (ys) => E.ols(X, ys), rng, B });
      specs.push(col("(2) Unit FE", fit.beta[0], r, true, false, false, fit.n, withinR2(yd, fit.resid)));
    }
    // (3) Unit + Time FE (preferred)
    let preferred;
    {
      const yd = E.twoWayDemean(c.y, c.unit, c.time), Dd = E.twoWayDemean(c.D, c.unit, c.time);
      const X = Dd.map(d => [d]); const fit = E.ols(X, yd);
      const dfAbsorb = nUnique(c.unit) + nUnique(c.time) - 1;
      const buildFitY = (ys) => E.ols(X, ys);
      let r;
      if (method === "ri") r = riForDiD(data, fit.beta[0], rng, B);
      else r = seForTarget({ fit, method, clusterIds: c.group, dfAbsorb, target: 0, buildFitY, rng, B });
      preferred = col("(3) Two-way FE", fit.beta[0], r, true, true, false, fit.n, withinR2(yd, fit.resid));
      preferred.preferred = true;
      specs.push(preferred);
    }
    // (4) Two-way FE + covariate
    {
      const yd = E.twoWayDemean(c.y, c.unit, c.time), Dd = E.twoWayDemean(c.D, c.unit, c.time), xd = E.twoWayDemean(c.x, c.unit, c.time);
      const X = Dd.map((d, i) => [d, xd[i]]); const fit = E.ols(X, yd);
      const dfAbsorb = nUnique(c.unit) + nUnique(c.time) - 1;
      const r = seForTarget({ fit, method: method === "ri" ? "cluster" : method, clusterIds: c.group, dfAbsorb, target: 0, buildFitY: (ys) => E.ols(X, ys), rng, B });
      specs.push(col("(4) + Covariate", fit.beta[0], r, true, true, true, fit.n, withinR2(yd, fit.resid)));
    }
    return { columns: specs, target: "Treated × Post", trueEffect: data.params.tau, dvMean: mean(c.y) };
  }
  function col(label, coef, r, feUnit, feTime, ctrl, n, r2) {
    return { label, coef, se: r.se, lo: r.lo, hi: r.hi, p: r.p, stars: r.stars, feUnit, feTime, ctrl, n, r2 };
  }

  /* ---- Randomization inference for the two-way-FE ATT ------------------- */
  function riForDiD(data, betaHat, rng, B) {
    const c = data.cols;
    const adoptByGroup = data.adopt.slice();
    const buildFromAdopt = (perm) => {
      const Dperm = c.time.map((t, i) => (perm[c.group[i]] !== Infinity && t >= perm[c.group[i]] ? 1 : 0));
      const yd = E.twoWayDemean(c.y, c.unit, c.time), Dd = E.twoWayDemean(Dperm, c.unit, c.time);
      const X = Dd.map(d => [d]); return E.ols(X, yd).beta[0];
    };
    const ri = E.randInference((perm) => buildFromAdopt(perm), betaHat, adoptByGroup, B || 499, rng);
    const se = NaN; // RI yields a p-value, not an SE
    return { se: se, lo: NaN, hi: NaN, p: ri.pvalue, stars: E.starsFromP(ri.pvalue) };
  }

  /* ====================================================================== *
   *  DiD: robustness = SAME estimate under every SE definition (the SE lab)
   * ====================================================================== */
  function didRobustnessSE(data, rng, B) {
    const c = data.cols;
    const yd = E.twoWayDemean(c.y, c.unit, c.time), Dd = E.twoWayDemean(c.D, c.unit, c.time);
    const X = Dd.map(d => [d]); const fit = E.ols(X, yd);
    const dfAbsorb = nUnique(c.unit) + nUnique(c.time) - 1;
    const buildFitY = (ys) => E.ols(X, ys);
    const methods = ["iid", "HC1", "HC3", "cluster", "wildboot"];
    const rows = methods.map(m => {
      const r = seForTarget({ fit, method: m, clusterIds: c.group, dfAbsorb, target: 0, buildFitY, rng, B });
      return { method: m, label: METHOD_LABEL[m], coef: fit.beta[0], se: r.se, lo: r.lo, hi: r.hi, p: r.p, stars: r.stars, note: r.note };
    });
    // add randomization inference row (p-value only)
    const ri = riForDiD(data, fit.beta[0], rng, B);
    rows.push({ method: "ri", label: METHOD_LABEL.ri, coef: fit.beta[0], se: NaN, lo: NaN, hi: NaN, p: ri.p, stars: ri.stars, note: "p-value only" });
    return { rows, coef: fit.beta[0], nClusters: nUnique(c.group) };
  }

  /* ====================================================================== *
   *  DiD: heterogeneity by covariate median → forest plot
   * ====================================================================== */
  function didHeterogeneity(data, method, rng, B) {
    const c = data.cols;
    const med = median(c.x);
    const groupsDef = [{ label: "Below median x", keep: (i) => c.x[i] <= med }, { label: "Above median x", keep: (i) => c.x[i] > med }];
    const out = groupsDef.map(gd => {
      const idx = []; for (let i = 0; i < data.n; i++) if (gd.keep(i)) idx.push(i);
      const yy = idx.map(i => c.y[i]), uu = idx.map(i => c.unit[i]), tt = idx.map(i => c.time[i]), DD = idx.map(i => c.D[i]), gg = idx.map(i => c.group[i]);
      const yd = E.twoWayDemean(yy, uu, tt), Dd = E.twoWayDemean(DD, uu, tt);
      const X = Dd.map(d => [d]); const fit = E.ols(X, yd);
      const dfAbsorb = nUnique(uu) + nUnique(tt) - 1;
      const r = seForTarget({ fit, method: method === "ri" || method === "wildboot" ? "cluster" : method, clusterIds: gg, dfAbsorb, target: 0, buildFitY: (ys) => E.ols(X, ys), rng, B });
      return { label: gd.label, coef: fit.beta[0], se: r.se, lo: r.lo, hi: r.hi, n: idx.length };
    });
    return { groups: out };
  }

  /* ====================================================================== *
   *  DiD: summary statistics (Table 1)
   * ====================================================================== */
  function didSummary(data) {
    const c = data.cols;
    const treatedIdx = [], controlIdx = [];
    for (let i = 0; i < data.n; i++) (c.everTreated[i] ? treatedIdx : controlIdx).push(i);
    const rowFor = (name, vals) => {
      const t = treatedIdx.map(i => vals[i]), k = controlIdx.map(i => vals[i]);
      return { name, allMean: mean(vals), allSd: sd(vals), tMean: mean(t), tSd: sd(t), cMean: mean(k), cSd: sd(k) };
    };
    return {
      rows: [rowFor("Outcome y", c.y), rowFor("Covariate x", c.x), rowFor("Treated × Post (share)", c.D)],
      nTreatedUnits: nUnique(treatedIdx.map(i => c.unit[i])),
      nControlUnits: nUnique(controlIdx.map(i => c.unit[i])),
      G: data.G, T: data.T, n: data.n,
    };
  }

  /* ====================================================================== *
   *  RCT stub: balance + ITT
   * ====================================================================== */
  function rctBalance(data) {
    const c = data.cols;
    const t = [], k = [];
    for (let i = 0; i < data.n; i++) (c.W[i] ? t : k).push(i);
    const nx = c.X[0].length;
    const rows = [];
    for (let j = 0; j < nx; j++) {
      const xt = t.map(i => c.X[i][j]), xk = k.map(i => c.X[i][j]);
      const diff = mean(xt) - mean(xk);
      const sePooled = Math.sqrt(variance(xt) / xt.length + variance(xk) / xk.length);
      const p = E.pValueFromZ(diff / sePooled);
      const nd = diff / Math.sqrt((variance(xt) + variance(xk)) / 2); // normalized difference
      rows.push({ name: "X" + (j + 1), c1: mean(xt), c1sd: sd(xt), c0: mean(xk), c0sd: sd(xk), diff, p, nd });
    }
    return { rows, nT: t.length, nC: k.length };
  }
  function rctITT(data, method, rng, B) {
    const c = data.cols;
    const X = c.W.map((w, i) => [1, w].concat(c.X[i]));
    const fit = E.ols(X, c.y);
    const r = seForTarget({ fit, method: method === "ri" || method === "wildboot" ? "cluster" : method, clusterIds: c.cluster, dfAbsorb: 0, target: 1, buildFitY: (ys) => E.ols(X, ys), rng, B });
    return { coef: fit.beta[1], se: r.se, lo: r.lo, hi: r.hi, p: r.p, stars: r.stars };
  }

  /* ====================================================================== *
   *  IV stub: first stage / reduced form / 2SLS (just-identified)
   * ====================================================================== */
  function ivFit(data) {
    const c = data.cols;
    // first stage: d ~ a + z
    const Xz = c.z.map(z => [1, z]);
    const fs = E.ols(Xz, c.d);
    const Vfs = E.vcovHC(fs, "HC1");
    const Ffirst = (fs.beta[1] * fs.beta[1]) / Vfs[1][1];
    // reduced form: y ~ a + z
    const rf = E.ols(Xz, c.y);
    // 2SLS just-identified: beta = cov(y,z)/cov(d,z)
    const beta2sls = cov(c.y, c.z) / cov(c.d, c.z);
    return {
      firstStage: { slope: fs.beta[1], F: Ffirst },
      reducedForm: { slope: rf.beta[1] },
      tsls: beta2sls,
      scatter: c.z.map((z, i) => ({ z, d: c.d[i] })),
      fsLine: { a: fs.beta[0], b: fs.beta[1] },
    };
  }

  /* ====================================================================== *
   *  RDD stub: binned means + local-linear jump
   * ====================================================================== */
  function rddFit(data, nbins) {
    nbins = nbins || 20;
    const c = data.cols, cut = data.params.cutoff;
    const lo = Math.min(...c.run), hi = Math.max(...c.run), w = (hi - lo) / nbins;
    const bins = Array.from({ length: nbins }, () => ({ n: 0, sy: 0, sx: 0 }));
    for (let i = 0; i < data.n; i++) {
      let b = Math.floor((c.run[i] - lo) / w); if (b >= nbins) b = nbins - 1; if (b < 0) b = 0;
      bins[b].n++; bins[b].sy += c.y[i]; bins[b].sx += c.run[i];
    }
    const binned = bins.filter(b => b.n > 0).map(b => ({ x: b.sx / b.n, y: b.sy / b.n, side: b.sx / b.n >= cut ? 1 : 0 }));
    // local linear each side within bandwidth
    const fitSide = (sidePred) => {
      const X = [], Y = [];
      for (let i = 0; i < data.n; i++) if (sidePred(c.run[i])) { X.push([1, c.run[i] - cut]); Y.push(c.y[i]); }
      return E.ols(X, Y).beta[0]; // intercept = level at cutoff
    };
    const right = fitSide(r => r >= cut && r <= cut + 0.3);
    const left = fitSide(r => r < cut && r >= cut - 0.3);
    return { binned, jump: right - left, cutoff: cut, left, right, trueJump: data.params.jump };
  }

  /* ====================================================================== *
   *  RCT: ITT progressive columns (unadjusted → +covariates → +strata FE)
   * ====================================================================== */
  function rctMainITT(data, method, rng, B) {
    const c = data.cols;
    const obs = []; for (let i = 0; i < data.n; i++) if (!c.attrit[i]) obs.push(i); // observed (post-attrition)
    const y = obs.map(i => c.y[i]), W = obs.map(i => c.W[i]), cl = obs.map(i => c.cluster[i]),
      strata = obs.map(i => c.strata[i]), X = obs.map(i => c.X[i]);
    const m = (method === "ri" ? "cluster" : method);
    const ctrlMean = mean(y.filter((v, i) => W[i] === 0));
    function run(buildX) {
      const Xm = buildX(); const fit = E.ols(Xm, y);
      const r = seForTarget({ fit, method: m, clusterIds: cl, dfAbsorb: 0, target: 1, buildFitY: ys => E.ols(Xm, ys), rng, B });
      return { coef: fit.beta[1], se: r.se, p: r.p, stars: r.stars, n: fit.n };
    }
    const lv = uniqueSorted(strata);
    const c1 = run(() => W.map(w => [1, w]));
    const c2 = run(() => W.map((w, i) => [1, w].concat(X[i])));
    const c3 = run(() => W.map((w, i) => { const row = [1, w].concat(X[i]); for (let s = 1; s < lv.length; s++) row.push(strata[i] === lv[s] ? 1 : 0); return row; }));
    const cols = [
      { label: "(1)", coef: c1.coef, se: c1.se, p: c1.p, stars: c1.stars },
      { label: "(2)", coef: c2.coef, se: c2.se, p: c2.p, stars: c2.stars },
      { label: "(3)", coef: c3.coef, se: c3.se, p: c3.p, stars: c3.stars, pref: true },
    ];
    const rows = [
      { name: "Covariates", vals: ["No", "Yes", "Yes"] },
      { name: "Strata FE", vals: ["No", "No", "Yes"] },
      { name: "Observations", vals: [c1.n, c2.n, c3.n] },
      { name: "Control mean", vals: [f3(ctrlMean), f3(ctrlMean), f3(ctrlMean)] },
    ];
    const trueITT = data.params.tau * data.params.complyRate;
    const note = "Intention-to-treat (assignment). SEs " + (METHOD_LABEL[m] || m).toLowerCase() + ". True ITT (DGP) = " + f3(trueITT) + ".";
    return { spec: { target: "Assigned to treatment", columns: cols, rows: rows, note: note }, caption: "RCT: intention-to-treat estimates" };
  }

  function rctLATE(data) {
    const c = data.cols;
    const obs = []; for (let i = 0; i < data.n; i++) if (!c.attrit[i]) obs.push(i);
    const y = obs.map(i => c.y[i]), W = obs.map(i => c.W[i]), D = obs.map(i => c.D[i]);
    const late = cov(y, W) / cov(D, W);
    const fsT = mean(D.filter((v, i) => W[i] === 1)), fsC = mean(D.filter((v, i) => W[i] === 0));
    return { late, compliance: fsT - fsC, trueLATE: data.params.tau };
  }

  function rctAttrition(data) {
    const c = data.cols;
    const byArm = [{ arm: "Control", w: 0 }, { arm: "Treatment", w: 1 }].map(a => {
      let n = 0, nAtt = 0;
      for (let i = 0; i < data.n; i++) if (c.W[i] === a.w) { n++; nAtt += c.attrit[i]; }
      return { arm: a.arm, n, nAtt, rate: nAtt / n };
    });
    const X = c.W.map(w => [1, w]); const fit = E.ols(X, c.attrit);
    const V = E.vcovCluster(fit, c.cluster, 0); const se = Math.sqrt(Math.max(V[1][1], 0));
    return { byArm, diff: fit.beta[1], p: E.pValueFromZ(fit.beta[1] / se) };
  }

  function rctHeterogeneity(data, method, rng, B) {
    const c = data.cols;
    const obs = []; for (let i = 0; i < data.n; i++) if (!c.attrit[i]) obs.push(i);
    const x1 = obs.map(i => c.X[i][0]); const med = median(x1);
    const defs = [{ label: "Below median X1", keep: (k) => c.X[obs[k]][0] <= med }, { label: "Above median X1", keep: (k) => c.X[obs[k]][0] > med }];
    const m = (method === "ri" || method === "wildboot" ? "cluster" : method);
    const groups = defs.map(gd => {
      const idx = []; for (let k = 0; k < obs.length; k++) if (gd.keep(k)) idx.push(obs[k]);
      const y = idx.map(i => c.y[i]), W = idx.map(i => c.W[i]), cl = idx.map(i => c.cluster[i]);
      const Xm = W.map(w => [1, w]); const fit = E.ols(Xm, y);
      const r = seForTarget({ fit, method: m, clusterIds: cl, dfAbsorb: 0, target: 1, buildFitY: ys => E.ols(Xm, ys), rng, B });
      return { label: gd.label, coef: fit.beta[1], se: r.se, lo: r.lo, hi: r.hi, n: idx.length };
    });
    return { groups };
  }

  /* ====================================================================== *
   *  IV: OLS vs 2SLS (just-identified) + first stage / reduced form
   * ====================================================================== */
  function ivResults(data) {
    const c = data.cols;
    // OLS y ~ d
    const Xd = c.d.map(d => [1, d]); const ols = E.ols(Xd, c.y); const Vo = E.vcovHC(ols, "HC1");
    const seO = Math.sqrt(Vo[1][1]), pO = E.pValueFromZ(ols.beta[1] / seO);
    // First stage d ~ z ; reduced form y ~ z
    const Xz = c.z.map(z => [1, z]); const fs = E.ols(Xz, c.d); const Vf = E.vcovHC(fs, "HC1");
    const F = (fs.beta[1] * fs.beta[1]) / Vf[1][1];
    const rf = E.ols(Xz, c.y); const Vr = E.vcovHC(rf, "HC1"); const seR = Math.sqrt(Vr[1][1]);
    // 2SLS just-identified with robust SE: b = (Z'X)^-1 Z'y
    const Z = c.z.map(z => [1, z]); const Xm = c.d.map(d => [1, d]);
    const ZX = E.matMul(E.transpose(Z), Xm); const ZXinv = E.inverse(ZX);
    const b = E.matVec(ZXinv, E.matVec(E.transpose(Z), c.y));
    const e = c.y.map((yi, i) => yi - (b[0] + b[1] * c.d[i]));
    const meat = [[0, 0], [0, 0]];
    for (let i = 0; i < data.n; i++) { const zi = Z[i], w = e[i] * e[i]; for (let a = 0; a < 2; a++) for (let bb = 0; bb < 2; bb++) meat[a][bb] += w * zi[a] * zi[bb]; }
    const V2 = E.matMul(E.matMul(ZXinv, meat), E.transpose(ZXinv));
    const se2 = Math.sqrt(Math.max(V2[1][1], 0)), p2 = E.pValueFromZ(b[1] / se2);
    const cols = [
      { label: "(1) OLS", coef: ols.beta[1], se: seO, p: pO, stars: E.starsFromP(pO) },
      { label: "(2) 2SLS", coef: b[1], se: se2, p: p2, stars: E.starsFromP(p2), pref: true },
    ];
    const rows = [
      { name: "Estimator", vals: ["OLS", "2SLS (IV)"] },
      { name: "First-stage F", vals: ["—", f1(F)] },
      { name: "Observations", vals: [data.n, data.n] },
    ];
    const note = "OLS is biased by endogeneity; 2SLS uses z. Reduced-form slope = " + f3(rf.beta[1]) + ". True β (DGP) = " + f3(data.params.beta) + ". " + (F < 10 ? "First-stage F < 10 → weak instrument; use Anderson–Rubin CIs." : "First stage strong.");
    return { spec: { target: "Endogenous regressor d", columns: cols, rows: rows, note: note }, caption: "Instrumental-variables estimates", F, weak: F < 10 };
  }

  /* ====================================================================== *
   *  RDD: bandwidth × polynomial robustness + density / covariate placebo
   * ====================================================================== */
  function rdLocal(c, cut, h, poly, outcome) {
    const yv = outcome || c.y;
    const X = [], Y = [];
    for (let i = 0; i < yv.length; i++) {
      const rr = c.run[i] - cut; if (Math.abs(rr) > h) continue;
      const t = c.run[i] >= cut ? 1 : 0; const row = [1, t];
      for (let pp = 1; pp <= poly; pp++) { row.push(Math.pow(rr, pp)); row.push(t * Math.pow(rr, pp)); }
      X.push(row); Y.push(yv[i]);
    }
    const fit = E.ols(X, Y); const V = E.vcovHC(fit, "HC1");
    const jump = fit.beta[1], se = Math.sqrt(Math.max(V[1][1], 0));
    return { jump, se, p: E.pValueFromZ(jump / se), n: X.length };
  }
  function rddEstimates(data) {
    const c = data.cols, cut = data.params.cutoff;
    const specsDef = [
      { label: "(1)", h: 0.5, poly: 1 }, { label: "(2)", h: 0.3, poly: 1 },
      { label: "(3)", h: 0.2, poly: 1 }, { label: "(4)", h: 0.3, poly: 2 },
    ];
    const est = specsDef.map(s => rdLocal(c, cut, s.h, s.poly));
    const cols = est.map((e, i) => ({ label: specsDef[i].label, coef: e.jump, se: e.se, p: e.p, stars: E.starsFromP(e.p), pref: i === 1 }));
    const rows = [
      { name: "Bandwidth (h)", vals: specsDef.map(s => f2(s.h)) },
      { name: "Polynomial", vals: specsDef.map(s => (s.poly === 1 ? "Linear" : "Quadratic")) },
      { name: "Effective N", vals: est.map(e => e.n) },
    ];
    // covariate continuity placebo (jump in a covariate should be ~0)
    const placebo = rdLocal(c, cut, 0.3, 1, c.cov);
    const note = "Local-polynomial RD with robust SEs. True jump (DGP) = " + f3(data.params.jump) + ". Covariate placebo jump = " + f3(placebo.jump) + " (p = " + f3(placebo.p) + ").";
    return { spec: { target: "Discontinuity at cutoff", columns: cols, rows: rows, note: note }, caption: "RD estimates across bandwidth and polynomial order", placebo };
  }
  function rddDensity(data) {
    const c = data.cols, cut = data.params.cutoff, h = 0.15;
    let left = 0, right = 0;
    for (let i = 0; i < data.n; i++) { const r = c.run[i]; if (r >= cut - h && r < cut) left++; else if (r >= cut && r <= cut + h) right++; }
    const logdiff = Math.log((right + 0.5) / (left + 0.5));
    const se = Math.sqrt(1 / (right + 0.5) + 1 / (left + 0.5));
    return { left, right, logdiff, z: logdiff / se, p: E.pValueFromZ(logdiff / se) };
  }

  /* ---- small stats ----------------------------------------------------- */
  function f3(v) { return (v == null || isNaN(v)) ? "—" : Number(v).toFixed(3); }
  function f2(v) { return (v == null || isNaN(v)) ? "—" : Number(v).toFixed(2); }
  function f1(v) { return (v == null || isNaN(v)) ? "—" : Number(v).toFixed(1); }
  function uniqueSorted(a) { return Array.from(new Set(a)).sort((x, y) => x - y); }
  function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
  function variance(a) { const m = mean(a); return a.reduce((s, v) => s + (v - m) * (v - m), 0) / Math.max(a.length - 1, 1); }
  function sd(a) { return Math.sqrt(variance(a)); }
  function cov(a, b) { const ma = mean(a), mb = mean(b); let s = 0; for (let i = 0; i < a.length; i++) s += (a[i] - ma) * (b[i] - mb); return s / (a.length - 1); }
  function median(a) { const s = a.slice().sort((x, y) => x - y); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }

  global.Models = {
    METHOD_LABEL,
    didSummary, didEventStudy, didMain, didRobustnessSE, didHeterogeneity,
    rctBalance, rctITT, rctMainITT, rctLATE, rctAttrition, rctHeterogeneity,
    ivFit, ivResults, rddFit, rddEstimates, rddDensity,
  };
})(typeof window !== "undefined" ? window : this);
