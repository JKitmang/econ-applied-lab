/* ============================================================================
   econ.js — Minimal, transparent econometrics engine (pure JS, no deps)
   ----------------------------------------------------------------------------
   Everything here is implemented from first principles so that every number
   the prototype shows is auditable. Scope kept small on purpose (small k).

   Exposed as a global `Econ` object (no modules, so it works from file://).
   ========================================================================== */
(function (global) {
  "use strict";

  /* ---------- Seeded RNG (mulberry32) → reproducibility ------------------- */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function makeRNG(seed) {
    const u = mulberry32(seed);
    let spare = null;
    function normal() {
      if (spare !== null) { const s = spare; spare = null; return s; }
      // Box–Muller
      let u1 = 0, u2 = 0;
      while (u1 <= 1e-12) u1 = u();
      u2 = u();
      const mag = Math.sqrt(-2 * Math.log(u1));
      spare = mag * Math.sin(2 * Math.PI * u2);
      return mag * Math.cos(2 * Math.PI * u2);
    }
    return {
      uniform: u,
      normal,
      rademacher: () => (u() < 0.5 ? -1 : 1),
      intBelow: (n) => Math.floor(u() * n),
    };
  }

  /* ---------- Tiny linear algebra ---------------------------------------- */
  function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }
  function transpose(A) {
    const r = A.length, c = A[0].length, T = [];
    for (let j = 0; j < c; j++) { T[j] = []; for (let i = 0; i < r; i++) T[j][i] = A[i][j]; }
    return T;
  }
  function matMul(A, B) {
    const r = A.length, n = B.length, c = B[0].length, M = [];
    for (let i = 0; i < r; i++) {
      M[i] = new Array(c).fill(0);
      for (let k = 0; k < n; k++) { const a = A[i][k]; if (a === 0) continue; for (let j = 0; j < c; j++) M[i][j] += a * B[k][j]; }
    }
    return M;
  }
  function matVec(A, x) { return A.map(row => dot(row, x)); }
  // Gauss–Jordan inverse (small symmetric PD matrices)
  function inverse(A) {
    const n = A.length;
    const M = A.map((row, i) => row.concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
    for (let col = 0; col < n; col++) {
      let piv = col;
      for (let r = col + 1; r < n; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
      if (Math.abs(M[piv][col]) < 1e-12) throw new Error("Matrix is singular (collinear regressors?)");
      [M[col], M[piv]] = [M[piv], M[col]];
      const d = M[col][col];
      for (let j = 0; j < 2 * n; j++) M[col][j] /= d;
      for (let r = 0; r < n; r++) {
        if (r === col) continue;
        const f = M[r][col];
        if (f === 0) continue;
        for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[col][j];
      }
    }
    return M.map(row => row.slice(n));
  }
  function diag(A) { return A.map((row, i) => row[i]); }

  /* ---------- OLS -------------------------------------------------------- */
  // X: n×k design matrix (caller includes intercept/dummies as needed). y: n.
  function ols(X, y) {
    const n = X.length, k = X[0].length;
    const Xt = transpose(X);
    const XtX = matMul(Xt, X);
    const XtXinv = inverse(XtX);
    const beta = matVec(XtXinv, matVec(Xt, y));
    const fitted = matVec(X, beta);
    const resid = y.map((yi, i) => yi - fitted[i]);
    const rss = dot(resid, resid);
    const ybar = y.reduce((s, v) => s + v, 0) / n;
    const tss = y.reduce((s, v) => s + (v - ybar) * (v - ybar), 0);
    return { X, y, n, k, beta, fitted, resid, rss, XtXinv, r2: tss > 0 ? 1 - rss / tss : 0 };
  }

  /* ---------- Variance estimators (return full k×k vcov) ----------------- */
  function vcovClassical(fit) {
    const { rss, n, k, XtXinv } = fit;
    const s2 = rss / Math.max(n - k, 1);
    return XtXinv.map(row => row.map(v => v * s2));
  }
  // meat = sum_i w_i * x_i x_i' ; sandwich = XtXinv * meat * XtXinv
  function sandwich(fit, weights) {
    const { X, k, XtXinv } = fit;
    const meat = Array.from({ length: k }, () => new Array(k).fill(0));
    for (let i = 0; i < X.length; i++) {
      const xi = X[i], w = weights[i];
      for (let a = 0; a < k; a++) { const xa = xi[a]; if (xa === 0) continue; for (let b = 0; b < k; b++) meat[a][b] += w * xa * xi[b]; }
    }
    return matMul(matMul(XtXinv, meat), XtXinv);
  }
  function hatValues(fit) {
    const { X, XtXinv, n } = fit;
    const h = new Array(n);
    for (let i = 0; i < n; i++) { const xi = X[i]; h[i] = dot(xi, matVec(XtXinv, xi)); }
    return h;
  }
  // type: "HC0" | "HC1" | "HC2" | "HC3"
  function vcovHC(fit, type) {
    const { resid, n, k } = fit;
    let w;
    if (type === "HC0" || type === "HC1") {
      w = resid.map(e => e * e);
    } else {
      const h = hatValues(fit);
      if (type === "HC2") w = resid.map((e, i) => (e * e) / (1 - h[i]));
      else /* HC3 */ w = resid.map((e, i) => (e * e) / ((1 - h[i]) * (1 - h[i])));
    }
    const V = sandwich(fit, w);
    if (type === "HC1") { const c = n / (n - k); return V.map(r => r.map(v => v * c)); }
    return V;
  }
  // Cluster-robust CR1. clusterIds: length n. dfAbsorb = #FE params absorbed (for k).
  function vcovCluster(fit, clusterIds, dfAbsorb) {
    const { X, resid, n, k, XtXinv } = fit;
    const kEff = k + (dfAbsorb || 0);
    const groups = {};
    for (let i = 0; i < n; i++) (groups[clusterIds[i]] || (groups[clusterIds[i]] = [])).push(i);
    const meat = Array.from({ length: k }, () => new Array(k).fill(0));
    let G = 0;
    for (const g in groups) {
      G++;
      const idx = groups[g];
      const s = new Array(k).fill(0);              // X_g' e_g
      for (const i of idx) { const e = resid[i], xi = X[i]; for (let a = 0; a < k; a++) s[a] += xi[a] * e; }
      for (let a = 0; a < k; a++) { if (s[a] === 0) continue; for (let b = 0; b < k; b++) meat[a][b] += s[a] * s[b]; }
    }
    const c = (G / (G - 1)) * ((n - 1) / (n - kEff)); // CR1 small-sample correction
    const V = matMul(matMul(XtXinv, meat), XtXinv).map(r => r.map(v => v * c));
    V._G = G;
    return V;
  }

  /* ---------- Wild cluster bootstrap (percentile-t, unrestricted) --------
     Returns a CI for a single target coefficient and an "implied SE" so it can
     be displayed in a table cell. Demonstrates CI widening with few clusters. */
  function wildClusterBootCI(buildFit, fit, clusterIds, targetIdx, B, rng, dfAbsorb, alpha) {
    alpha = alpha || 0.05;
    const betaHat = fit.beta[targetIdx];
    const Vcl = vcovCluster(fit, clusterIds, dfAbsorb);
    const seHat = Math.sqrt(Math.max(Vcl[targetIdx][targetIdx], 0));
    const { fitted, resid } = fit;
    const n = fit.n;
    const groups = {};
    for (let i = 0; i < n; i++) (groups[clusterIds[i]] || (groups[clusterIds[i]] = [])).push(i);
    const groupKeys = Object.keys(groups);
    const tstars = [];
    for (let b = 0; b < B; b++) {
      const w = {}; for (const g of groupKeys) w[g] = rng.rademacher();
      const ystar = new Array(n);
      for (let i = 0; i < n; i++) ystar[i] = fitted[i] + w[clusterIds[i]] * resid[i];
      let f2;
      try { f2 = buildFit(ystar); } catch (e) { continue; }
      const Vb = vcovCluster(f2, clusterIds, dfAbsorb);
      const seb = Math.sqrt(Math.max(Vb[targetIdx][targetIdx], 1e-18));
      tstars.push((f2.beta[targetIdx] - betaHat) / seb);
    }
    tstars.sort((a, b) => a - b);
    const q = (p) => {
      if (tstars.length === 0) return 1.96;
      const pos = p * (tstars.length - 1), lo = Math.floor(pos), hi = Math.ceil(pos);
      return tstars[lo] + (tstars[hi] - tstars[lo]) * (pos - lo);
    };
    const ciLo = betaHat - seHat * q(1 - alpha / 2);
    const ciHi = betaHat - seHat * q(alpha / 2);
    // p-value: share of |t*| exceeding |t_hat| with t_hat = betaHat/seHat
    const tHat = betaHat / seHat;
    let extreme = 0; for (const t of tstars) if (Math.abs(t) >= Math.abs(tHat)) extreme++;
    const pvalue = (extreme + 1) / (tstars.length + 1);
    return { beta: betaHat, ci: [ciLo, ciHi], impliedSE: (ciHi - ciLo) / (2 * 1.959964), pvalue, seCluster: seHat };
  }

  /* ---------- Randomization inference ------------------------------------
     Permute a cluster-level binary treatment, recompute the target coef under
     the sharp null, build a two-sided p-value. */
  function randInference(buildFitFromTreat, betaHat, clusterTreatVector, B, rng) {
    const T = clusterTreatVector.slice();
    let extreme = 0, used = 0;
    for (let b = 0; b < B; b++) {
      // Fisher–Yates shuffle of the treatment labels across clusters
      for (let i = T.length - 1; i > 0; i--) { const j = rng.intBelow(i + 1); const t = T[i]; T[i] = T[j]; T[j] = t; }
      let bperm;
      try { bperm = buildFitFromTreat(T); } catch (e) { continue; }
      used++;
      if (Math.abs(bperm) >= Math.abs(betaHat)) extreme++;
    }
    return { pvalue: (extreme + 1) / (used + 1), draws: used };
  }

  /* ---------- Two-way within transformation (balanced panel) ------------- */
  // Demean a vector by unit and by period, add back grand mean. Exact two-way FE.
  function twoWayDemean(v, unitId, timeId) {
    const n = v.length;
    const uSum = {}, uCnt = {}, tSum = {}, tCnt = {};
    let grand = 0;
    for (let i = 0; i < n; i++) {
      const ui = unitId[i], ti = timeId[i];
      uSum[ui] = (uSum[ui] || 0) + v[i]; uCnt[ui] = (uCnt[ui] || 0) + 1;
      tSum[ti] = (tSum[ti] || 0) + v[i]; tCnt[ti] = (tCnt[ti] || 0) + 1;
      grand += v[i];
    }
    grand /= n;
    const out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = v[i] - uSum[unitId[i]] / uCnt[unitId[i]] - tSum[timeId[i]] / tCnt[timeId[i]] + grand;
    return out;
  }

  /* ---------- Helpers ---------------------------------------------------- */
  // Normal CDF / two-sided p-value from z
  function normCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
  function erf(x) {
    const t = 1 / (1 + 0.3275911 * Math.abs(x));
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return x >= 0 ? y : -y;
  }
  function pValueFromZ(z) { return 2 * (1 - normCdf(Math.abs(z))); }
  function starsFromP(p) { return p < 0.01 ? "***" : p < 0.05 ? "**" : p < 0.1 ? "*" : ""; }
  function ones(n) { return Array.from({ length: n }, () => [1]); }
  function colBind() { const cols = Array.from(arguments); const n = cols[0].length; const out = []; for (let i = 0; i < n; i++) { out[i] = []; for (const c of cols) { if (Array.isArray(c[i])) out[i].push(...c[i]); else out[i].push(c[i]); } } return out; }

  global.Econ = {
    makeRNG, mulberry32,
    transpose, matMul, matVec, inverse, diag, dot, colBind, ones,
    ols, vcovClassical, vcovHC, vcovCluster, hatValues, sandwich,
    wildClusterBootCI, randInference, twoWayDemean,
    normCdf, pValueFromZ, starsFromP,
  };
})(typeof window !== "undefined" ? window : this);
