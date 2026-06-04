/* ============================================================================
   render.js — Pure renderers: tables (HTML + LaTeX export) and hand-rolled SVG
   figures. No dependencies. `view` is "paper" | "present" and switches theme.
   Exposed as global `Render`.
   ========================================================================== */
(function (global) {
  "use strict";

  /* ---- number formatting ----------------------------------------------- */
  function f(x, d) { if (x == null || isNaN(x)) return "—"; return Number(x).toFixed(d == null ? 3 : d); }
  function coefCell(coef, stars, d) { return f(coef, d) + (stars ? '<sup class="st">' + stars + "</sup>" : ""); }
  function paren(x, d) { return x == null || isNaN(x) ? "" : "(" + f(x, d) + ")"; }

  const PALETTE = {
    paper:   { pt: "#222", ci: "#555", zero: "#222", onset: "#999", grid: "#e3e3e3", pre: "#f4f4f4", post: "#ededed", treated: "#222", control: "#9a9a9a", truth: "#bbb", axis: "#333" },
    present: { pt: "#0072B2", ci: "#0072B2", zero: "#333", onset: "#D55E00", grid: "#ececec", pre: "#fff7e6", post: "#eaf3fb", treated: "#0072B2", control: "#9a9a9a", truth: "#009E73", axis: "#333" },
  };
  const pal = (v) => PALETTE[v] || PALETTE.paper;

  /* ====================================================================== *
   *  TABLES (HTML)
   * ====================================================================== */
  function table1(sum, view) {
    const d = 2;
    let h = '<table class="econ-table ' + view + '"><thead><tr>' +
      "<th>Variable</th><th>All</th><th>Ever-treated</th><th>Never-treated</th></tr></thead><tbody>";
    for (const r of sum.rows) {
      h += "<tr><td class='vn'>" + r.name + "</td>" +
        "<td>" + f(r.allMean, d) + "<br><span class='sd'>" + paren(r.allSd, d) + "</span></td>" +
        "<td>" + f(r.tMean, d) + "<br><span class='sd'>" + paren(r.tSd, d) + "</span></td>" +
        "<td>" + f(r.cMean, d) + "<br><span class='sd'>" + paren(r.cSd, d) + "</span></td></tr>";
    }
    h += "</tbody><tfoot><tr><td colspan='4' class='note'>" +
      "Means with SD in parentheses. Units: " + sum.nTreatedUnits + " ever-treated, " +
      sum.nControlUnits + " never-treated; " + sum.G + " clusters × " + sum.T + " periods (N = " + sum.n + ").</td></tr></tfoot></table>";
    return h;
  }

  function mainTable(main, view) {
    const present = view === "present";
    const d = 3;
    const cols = main.columns;
    let h = '<table class="econ-table ' + view + '"><thead><tr><th></th>';
    cols.forEach(c => { h += "<th" + (c.preferred ? " class='pref'" : "") + ">" + c.label + "</th>"; });
    h += "</tr></thead><tbody>";
    // coefficient row
    h += "<tr class='coefrow'><td class='vn'>" + main.target + "</td>";
    cols.forEach(c => { h += "<td" + (c.preferred ? " class='pref'" : "") + ">" + coefCell(c.coef, c.stars, d) + "</td>"; });
    h += "</tr><tr class='serow'><td></td>";
    cols.forEach(c => { h += "<td" + (c.preferred ? " class='pref'" : "") + ">" + (isNaN(c.se) ? "[p=" + f(c.p, 3) + "]" : paren(c.se, d)) + "</td>"; });
    h += "</tr>";
    // FE / controls rows
    const yn = (b) => (b ? "Yes" : "No");
    h += feRow("Unit FE", cols.map(c => yn(c.feUnit)));
    h += feRow("Time FE", cols.map(c => yn(c.feTime)));
    h += feRow("Covariate", cols.map(c => yn(c.ctrl)));
    h += feRow("Observations", cols.map(c => String(c.n)));
    h += feRow("R² (within)", cols.map(c => f(c.r2, 3)));
    h += "<tr class='dvmean'><td class='vn'>Mean of DV</td><td colspan='" + cols.length + "'>" + f(main.dvMean, 3) +
      "  ·  True effect (DGP): " + f(main.trueEffect, 3) + "</td></tr>";
    h += "</tbody><tfoot><tr><td colspan='" + (cols.length + 1) +
      "' class='note'>Standard errors in parentheses (or randomization-inference p-value in brackets). " +
      "<sup>*</sup> p&lt;0.10, <sup>**</sup> p&lt;0.05, <sup>***</sup> p&lt;0.01.</td></tr></tfoot></table>";
    return h;
    function feRow(name, vals) {
      let r = "<tr class='ferow'><td class='vn'>" + name + "</td>";
      vals.forEach((v, i) => { r += "<td" + (cols[i].preferred ? " class='pref'" : "") + ">" + v + "</td>"; });
      return r + "</tr>";
    }
  }

  function robustnessTable(rob, view) {
    const d = 3;
    let h = '<table class="econ-table ' + view + '"><thead><tr>' +
      "<th>SE definition</th><th>Estimate</th><th>Std. error</th><th>95% CI</th><th>p-value</th><th>Note</th></tr></thead><tbody>";
    for (const r of rob.rows) {
      const ci = isNaN(r.lo) ? "—" : "[" + f(r.lo, 2) + ", " + f(r.hi, 2) + "]";
      h += "<tr><td class='vn'>" + r.label + "</td>" +
        "<td>" + coefCell(r.coef, r.stars, d) + "</td>" +
        "<td>" + (isNaN(r.se) ? "—" : f(r.se, d)) + "</td>" +
        "<td>" + ci + "</td>" +
        "<td>" + f(r.p, 3) + "</td>" +
        "<td class='note'>" + (r.note || "") + "</td></tr>";
    }
    h += "</tbody><tfoot><tr><td colspan='6' class='note'>Same point estimate (" + f(rob.coef, 3) +
      "), " + rob.nClusters + " clusters; only the <em>inference</em> changes. With few clusters the wild bootstrap CI is wider than the naive cluster CI.</td></tr></tfoot></table>";
    return h;
  }

  function balanceTable(bal, view) {
    const d = 3;
    let h = '<table class="econ-table ' + view + '"><thead><tr>' +
      "<th>Variable</th><th>Control</th><th>Treatment</th><th>Diff.</th><th>Norm. diff.</th><th>p</th></tr></thead><tbody>";
    for (const r of bal.rows) {
      const flag = Math.abs(r.nd) > 0.25 ? " class='warn'" : "";
      h += "<tr><td class='vn'>" + r.name + "</td>" +
        "<td>" + f(r.c0, d) + " <span class='sd'>" + paren(r.c0sd, 2) + "</span></td>" +
        "<td>" + f(r.c1, d) + " <span class='sd'>" + paren(r.c1sd, 2) + "</span></td>" +
        "<td>" + f(r.diff, d) + "</td><td" + flag + ">" + f(r.nd, 3) + "</td><td>" + f(r.p, 3) + "</td></tr>";
    }
    h += "</tbody><tfoot><tr><td colspan='6' class='note'>N: " + bal.nC + " control, " + bal.nT +
      " treatment. |Norm. diff.| > 0.25 flagged (Imbens–Wooldridge).</td></tr></tfoot></table>";
    return h;
  }

  // Generic coefficient table: columns of estimates + indicator rows.
  // spec = { target, columns:[{label,coef,se,p,stars,pref,brackets}], rows:[{name,vals:[]}], note }
  function resultsTable(spec, view) {
    const d = 3, cols = spec.columns;
    let h = '<table class="econ-table ' + view + '"><thead><tr><th></th>';
    cols.forEach(c => { h += "<th" + (c.pref ? " class='pref'" : "") + ">" + c.label + "</th>"; });
    h += "</tr></thead><tbody>";
    h += "<tr class='coefrow'><td class='vn'>" + spec.target + "</td>";
    cols.forEach(c => { h += "<td" + (c.pref ? " class='pref'" : "") + ">" + coefCell(c.coef, c.stars, d) + "</td>"; });
    h += "</tr><tr class='serow'><td></td>";
    cols.forEach(c => { h += "<td" + (c.pref ? " class='pref'" : "") + ">" + (isNaN(c.se) ? (c.brackets || "") : paren(c.se, d)) + "</td>"; });
    h += "</tr>";
    (spec.rows || []).forEach(r => {
      h += "<tr class='ferow'><td class='vn'>" + r.name + "</td>";
      r.vals.forEach((v, i) => { h += "<td" + (cols[i].pref ? " class='pref'" : "") + ">" + v + "</td>"; });
      h += "</tr>";
    });
    h += "</tbody><tfoot><tr><td colspan='" + (cols.length + 1) + "' class='note'>" +
      (spec.note || "") + " <sup>*</sup> p&lt;0.10, <sup>**</sup> p&lt;0.05, <sup>***</sup> p&lt;0.01.</td></tr></tfoot></table>";
    return h;
  }

  function attritionTable(att, view) {
    let h = '<table class="econ-table ' + view + '"><thead><tr>' +
      "<th>Group</th><th>N assigned</th><th>Attrited</th><th>Attrition rate</th></tr></thead><tbody>";
    att.byArm.forEach(a => {
      h += "<tr><td class='vn'>" + a.arm + "</td><td>" + a.n + "</td><td>" + a.nAtt + "</td><td>" + f(a.rate, 3) + "</td></tr>";
    });
    const flag = att.p < 0.10 ? " class='warn'" : "";
    h += "</tbody><tfoot><tr><td colspan='4' class='note'>Differential attrition (T − C): <strong" + flag + ">" +
      f(att.diff, 3) + "</strong> (p = " + f(att.p, 3) + "). " +
      (att.p < 0.10 ? "Differential attrition detected — consider Lee bounds." : "No significant differential attrition.") +
      "</td></tr></tfoot></table>";
    return h;
  }

  /* ====================================================================== *
   *  LaTeX EXPORT (booktabs / threeparttable — matches reg-table skill)
   * ====================================================================== */
  function latexResults(spec, caption) {
    const d = 3, cols = spec.columns, k = cols.length;
    let L = "\\begin{table}[htbp]\\centering\n\\caption{" + (caption || "Results") + "}\n";
    L += "\\begin{threeparttable}\n\\begin{tabular}{l" + "c".repeat(k) + "}\n\\toprule\n";
    L += " & " + cols.map((c, i) => "(" + (i + 1) + ")").join(" & ") + " \\\\\n";
    L += " & " + cols.map(c => escapeAmp(c.label)).join(" & ") + " \\\\\n\\midrule\n";
    L += escapeAmp(spec.target) + " & " + cols.map(c => "$" + f(c.coef, d) + (c.stars ? "^{" + c.stars + "}" : "") + "$").join(" & ") + " \\\\\n";
    L += " & " + cols.map(c => (isNaN(c.se) ? "" : "(" + f(c.se, d) + ")")).join(" & ") + " \\\\\n\\addlinespace\n";
    (spec.rows || []).forEach(r => { L += escapeAmp(r.name) + " & " + r.vals.map(escapeAmp).join(" & ") + " \\\\\n"; });
    L += "\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n";
    L += "\\item " + escapeAmp(spec.note || "") + "\n";
    L += "\\item \\textsuperscript{*} $p<0.10$, \\textsuperscript{**} $p<0.05$, \\textsuperscript{***} $p<0.01$.\n";
    L += "\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n";
    return L;
  }
  function latexMain(main) {
    const d = 3, cols = main.columns, k = cols.length;
    const star = (c) => c.stars;
    const cc = (c) => f(c.coef, d) + (star(c) ? "^{" + star(c) + "}" : "");
    let L = "";
    L += "\\begin{table}[htbp]\\centering\n\\caption{Main results: difference-in-differences}\n";
    L += "\\begin{threeparttable}\n\\begin{tabular}{l" + "c".repeat(k) + "}\n\\toprule\n";
    L += " & " + cols.map(c => c.label.replace(/[()]/g, "")).map(s => "(" + s.trim().split(" ")[0] + ")").join(" & ").replace(/\(\(/g, "(").replace(/\)\)/g, ")");
    L += " \\\\\n\\midrule\n";
    L += escapeAmp(main.target) + " & " + cols.map(c => "$" + cc(c) + "$").join(" & ") + " \\\\\n";
    L += " & " + cols.map(c => (isNaN(c.se) ? "[" + f(c.p, 3) + "]" : "(" + f(c.se, d) + ")")).join(" & ") + " \\\\\n";
    L += "\\addlinespace\n";
    L += "Unit FE & " + cols.map(c => (c.feUnit ? "Yes" : "No")).join(" & ") + " \\\\\n";
    L += "Time FE & " + cols.map(c => (c.feTime ? "Yes" : "No")).join(" & ") + " \\\\\n";
    L += "Covariate & " + cols.map(c => (c.ctrl ? "Yes" : "No")).join(" & ") + " \\\\\n";
    L += "Observations & " + cols.map(c => c.n).join(" & ") + " \\\\\n";
    L += "$R^2$ (within) & " + cols.map(c => f(c.r2, 3)).join(" & ") + " \\\\\n";
    L += "Mean of DV & \\multicolumn{" + k + "}{c}{" + f(main.dvMean, 3) + "} \\\\\n";
    L += "\\bottomrule\n\\end{tabular}\n";
    L += "\\begin{tablenotes}[flushleft]\\footnotesize\n";
    L += "\\item Standard errors in parentheses (randomization-inference $p$-value in brackets).\n";
    L += "\\item \\textsuperscript{*} $p<0.10$, \\textsuperscript{**} $p<0.05$, \\textsuperscript{***} $p<0.01$.\n";
    L += "\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n";
    return L;
  }
  function latexRobustness(rob) {
    const d = 3;
    let L = "\\begin{table}[htbp]\\centering\n\\caption{Inference under alternative standard-error definitions}\n";
    L += "\\begin{threeparttable}\n\\begin{tabular}{lcccc}\n\\toprule\n";
    L += "SE definition & Estimate & Std. error & 95\\% CI & $p$-value \\\\\n\\midrule\n";
    for (const r of rob.rows) {
      const ci = isNaN(r.lo) ? "--" : "[" + f(r.lo, 2) + ",\\," + f(r.hi, 2) + "]";
      L += escapeAmp(r.label) + " & $" + f(r.coef, d) + (r.stars ? "^{" + r.stars + "}" : "") + "$ & " +
        (isNaN(r.se) ? "--" : f(r.se, d)) + " & " + ci + " & " + f(r.p, 3) + " \\\\\n";
    }
    L += "\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n";
    L += "\\item Same point estimate; only inference varies. " + rob.nClusters + " clusters.\n";
    L += "\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n";
    return L;
  }
  function escapeAmp(s) { return String(s).replace(/&/g, "\\&").replace(/%/g, "\\%").replace(/×/g, "$\\times$"); }

  /* ====================================================================== *
   *  SVG FIGURES
   * ====================================================================== */
  function svgFrame(W, H, m) {
    return { W, H, m, iw: W - m.l - m.r, ih: H - m.t - m.b };
  }
  function axes(fr, xd, yd, view, opts) {
    opts = opts || {};
    const p = pal(view), m = fr.m;
    const sx = (x) => m.l + ((x - xd[0]) / (xd[1] - xd[0])) * fr.iw;
    const sy = (y) => m.t + (1 - (y - yd[0]) / (yd[1] - yd[0])) * fr.ih;
    let g = "";
    // y gridlines + ticks
    const yticks = ticks(yd[0], yd[1], 5);
    for (const t of yticks) {
      g += line(m.l, sy(t), m.l + fr.iw, sy(t), p.grid, 1);
      g += text(m.l - 8, sy(t) + 4, f(t, 1), p.axis, "end", opts.fs || 12);
    }
    const xticks = opts.xticks || ticks(xd[0], xd[1], 6);
    for (const t of xticks) {
      g += text(sx(t), m.t + fr.ih + 18, opts.xfmt ? opts.xfmt(t) : f(t, 0), p.axis, "middle", opts.fs || 12);
    }
    g += line(m.l, m.t + fr.ih, m.l + fr.iw, m.t + fr.ih, p.axis, 1.5); // x axis
    g += line(m.l, m.t, m.l, m.t + fr.ih, p.axis, 1.5);                 // y axis
    return { g, sx, sy };
  }

  function eventStudy(es, view, opts) {
    opts = opts || {};
    const present = view === "present";
    const W = present ? 720 : 560, H = present ? 420 : 340;
    const fr = svgFrame(W, H, { l: 56, r: 20, t: 22, b: 46 });
    const p = pal(view);
    const xs = es.points.map(pt => pt.e);
    const allY = es.points.flatMap(pt => [pt.lo, pt.hi, pt.coef]).concat([0]);
    const yd = padDomain(Math.min(...allY), Math.max(...allY), 0.15);
    const xd = [Math.min(...xs) - 0.5, Math.max(...xs) + 0.5];
    let body = "";
    // shaded post region
    const ax0 = 56, axw = W - 56 - 20;
    const sxRaw = (x) => fr.m.l + ((x - xd[0]) / (xd[1] - xd[0])) * fr.iw;
    body += rect(sxRaw(-0.5), fr.m.t, sxRaw(xd[1]) - sxRaw(-0.5), fr.ih, p.post, 1);
    const A = axes(fr, xd, yd, view, { fs: present ? 14 : 12 });
    body += A.g;
    // zero line + onset line
    body += line(fr.m.l, A.sy(0), fr.m.l + fr.iw, A.sy(0), p.zero, 1.2, "4 3");
    body += line(A.sx(-0.5), fr.m.t, A.sx(-0.5), fr.m.t + fr.ih, p.onset, 1.4, "5 4");
    // points + CI bars
    for (const pt of es.points) {
      const x = A.sx(pt.e);
      if (!pt.ref) {
        body += line(x, A.sy(pt.lo), x, A.sy(pt.hi), p.ci, present ? 2.2 : 1.6);
        body += line(x - 4, A.sy(pt.lo), x + 4, A.sy(pt.lo), p.ci, present ? 2.2 : 1.6);
        body += line(x - 4, A.sy(pt.hi), x + 4, A.sy(pt.hi), p.ci, present ? 2.2 : 1.6);
      }
      body += circle(x, A.sy(pt.coef), pt.ref ? 3 : (present ? 5 : 3.5), pt.ref ? "#fff" : p.pt, pt.ref ? p.onset : p.pt);
    }
    body += text(fr.m.l + fr.iw, fr.m.t + fr.ih + 36, "Event time (periods relative to adoption; −1 = reference)", p.axis, "end", present ? 13 : 11);
    return wrap(W, H, body, view);
  }

  function forest(het, view) {
    const present = view === "present";
    const W = present ? 680 : 540, H = 60 + het.groups.length * 54;
    const fr = svgFrame(W, H, { l: 150, r: 24, t: 16, b: 40 });
    const p = pal(view);
    const allX = het.groups.flatMap(g => [g.lo, g.hi]).concat([0]);
    const xd = padDomain(Math.min(...allX), Math.max(...allX), 0.15);
    const sx = (x) => fr.m.l + ((x - xd[0]) / (xd[1] - xd[0])) * fr.iw;
    let body = "";
    body += line(sx(0), fr.m.t, sx(0), fr.m.t + fr.ih, p.zero, 1.2, "4 3");
    for (const t of ticks(xd[0], xd[1], 5)) {
      body += line(sx(t), fr.m.t, sx(t), fr.m.t + fr.ih, p.grid, 1);
      body += text(sx(t), fr.m.t + fr.ih + 20, f(t, 1), p.axis, "middle", present ? 13 : 11);
    }
    het.groups.forEach((g, i) => {
      const y = fr.m.t + 28 + i * 50;
      body += text(16, y + 4, g.label + " (n=" + g.n + ")", p.axis, "start", present ? 14 : 12);
      body += line(sx(g.lo), y, sx(g.hi), y, p.ci, present ? 2.4 : 1.8);
      body += line(sx(g.lo), y - 4, sx(g.lo), y + 4, p.ci, present ? 2.4 : 1.8);
      body += line(sx(g.hi), y - 4, sx(g.hi), y + 4, p.ci, present ? 2.4 : 1.8);
      body += circle(sx(g.coef), y, present ? 6 : 4, p.pt, p.pt);
    });
    return wrap(W, H, body, view);
  }

  function trends(data, view) {
    const present = view === "present";
    const W = present ? 720 : 560, H = present ? 400 : 320;
    const fr = svgFrame(W, H, { l: 54, r: 110, t: 18, b: 42 });
    const p = pal(view), c = data.cols, T = data.T;
    const meanByGroupTime = (treated) => {
      const s = new Array(T).fill(0), n = new Array(T).fill(0);
      for (let i = 0; i < data.n; i++) if ((c.everTreated[i] === 1) === treated) { s[c.time[i]] += c.y[i]; n[c.time[i]]++; }
      return s.map((v, t) => (n[t] ? v / n[t] : null));
    };
    const tr = meanByGroupTime(true), co = meanByGroupTime(false);
    const allY = tr.concat(co).filter(v => v != null);
    const yd = padDomain(Math.min(...allY), Math.max(...allY), 0.1);
    const xd = [0, T - 1];
    const A = axes(fr, xd, yd, view, { fs: present ? 14 : 12, xticks: ticks(0, T - 1, Math.min(T, 7)) });
    let body = A.g;
    const polyline = (arr, color, dash) => {
      const pts = arr.map((v, t) => (v == null ? null : [A.sx(t), A.sy(v)])).filter(Boolean);
      let d = "M" + pts.map(pt => pt[0].toFixed(1) + "," + pt[1].toFixed(1)).join(" L");
      let s = '<path d="' + d + '" fill="none" stroke="' + color + '" stroke-width="' + (present ? 2.6 : 1.8) + '"' + (dash ? ' stroke-dasharray="' + dash + '"' : "") + "/>";
      for (const pt of pts) s += circle(pt[0], pt[1], present ? 3.5 : 2.5, color, color);
      return s;
    };
    body += polyline(tr, p.treated, "");
    body += polyline(co, p.control, "5 4");
    body += text(A.sx(T - 1) + 8, A.sy(tr[T - 1]), "Ever-treated", p.treated, "start", present ? 13 : 11);
    body += text(A.sx(T - 1) + 8, A.sy(co[T - 1]), "Never-treated", p.control, "start", present ? 13 : 11);
    body += text(fr.m.l + fr.iw / 2, fr.m.t + fr.ih + 36, "Period", p.axis, "middle", present ? 13 : 11);
    return wrap(W, H, body, view);
  }

  function firstStage(ivf, view) {
    const present = view === "present";
    const W = present ? 640 : 500, H = present ? 380 : 300;
    const fr = svgFrame(W, H, { l: 50, r: 18, t: 16, b: 42 });
    const p = pal(view);
    const xs = ivf.scatter.map(o => o.z), ys = ivf.scatter.map(o => o.d);
    const xd = padDomain(Math.min(...xs), Math.max(...xs), 0.05), yd = padDomain(Math.min(...ys), Math.max(...ys), 0.05);
    const A = axes(fr, xd, yd, view, { fs: present ? 14 : 12 });
    let body = A.g;
    for (const o of ivf.scatter) body += circle(A.sx(o.z), A.sy(o.d), present ? 2.6 : 1.8, p.control, "none", 0.5);
    body += line(A.sx(xd[0]), A.sy(ivf.fsLine.a + ivf.fsLine.b * xd[0]), A.sx(xd[1]), A.sy(ivf.fsLine.a + ivf.fsLine.b * xd[1]), p.pt, present ? 3 : 2);
    body += text(fr.m.l + 8, fr.m.t + 16, "First-stage F = " + f(ivf.firstStage.F, 1), p.pt, "start", present ? 15 : 12);
    body += text(fr.m.l + fr.iw / 2, fr.m.t + fr.ih + 36, "Instrument z", p.axis, "middle", present ? 13 : 11);
    return wrap(W, H, body, view);
  }

  function rdPlot(rd, view) {
    const present = view === "present";
    const W = present ? 680 : 540, H = present ? 400 : 320;
    const fr = svgFrame(W, H, { l: 50, r: 18, t: 16, b: 42 });
    const p = pal(view);
    const xs = rd.binned.map(b => b.x), ys = rd.binned.map(b => b.y);
    const xd = padDomain(Math.min(...xs), Math.max(...xs), 0.03), yd = padDomain(Math.min(...ys), Math.max(...ys), 0.1);
    const A = axes(fr, xd, yd, view, { fs: present ? 14 : 12 });
    let body = A.g;
    body += line(A.sx(rd.cutoff), fr.m.t, A.sx(rd.cutoff), fr.m.t + fr.ih, p.onset, 1.4, "5 4");
    for (const b of rd.binned) body += circle(A.sx(b.x), A.sy(b.y), present ? 4 : 2.8, b.side ? p.pt : p.control, b.side ? p.pt : p.control);
    body += text(A.sx(rd.cutoff) + 6, fr.m.t + 16, "Jump = " + f(rd.jump, 3), p.onset, "start", present ? 15 : 12);
    body += text(fr.m.l + fr.iw / 2, fr.m.t + fr.ih + 36, "Running variable (centered at cutoff)", p.axis, "middle", present ? 13 : 11);
    return wrap(W, H, body, view);
  }

  /* ---- SVG primitives -------------------------------------------------- */
  function wrap(W, H, body, view) {
    return '<svg viewBox="0 0 ' + W + " " + H + '" class="econ-fig ' + view + '" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg">' + body + "</svg>";
  }
  function line(x1, y1, x2, y2, color, w, dash) {
    return '<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) +
      '" stroke="' + color + '" stroke-width="' + (w || 1) + '"' + (dash ? ' stroke-dasharray="' + dash + '"' : "") + "/>";
  }
  function circle(cx, cy, r, fill, stroke, sw) {
    return '<circle cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + r + '" fill="' + fill + '"' +
      (stroke ? ' stroke="' + stroke + '" stroke-width="' + (sw || 1) + '"' : "") + "/>";
  }
  function rect(x, y, w, h, fill, op) { return '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + Math.max(w, 0).toFixed(1) + '" height="' + h.toFixed(1) + '" fill="' + fill + '" opacity="' + (op == null ? 1 : op) + '"/>'; }
  function text(x, y, s, color, anchor, fs) {
    return '<text x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" fill="' + color + '" text-anchor="' + (anchor || "start") +
      '" font-size="' + (fs || 12) + '" font-family="inherit">' + esc(s) + "</text>";
  }
  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }

  function ticks(lo, hi, n) {
    if (lo === hi) return [lo];
    const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const norm = step0 / mag, step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
    const out = []; let t = Math.ceil(lo / step) * step;
    for (; t <= hi + 1e-9; t += step) out.push(Math.round(t / step) * step);
    return out;
  }
  function padDomain(lo, hi, frac) { if (lo === hi) { lo -= 1; hi += 1; } const pad = (hi - lo) * frac; return [lo - pad, hi + pad]; }

  global.Render = {
    table1, mainTable, robustnessTable, balanceTable, resultsTable, attritionTable,
    latexMain, latexRobustness, latexResults,
    eventStudy, forest, trends, firstStage, rdPlot,
  };
})(typeof window !== "undefined" ? window : this);
