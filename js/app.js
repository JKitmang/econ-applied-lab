/* ============================================================================
   app.js — UI state + wiring. Generates data, runs models, draws exhibit cards.
   Depends on globals: Econ, DGP, Models, Render, TAXONOMY.
   ========================================================================== */
(function () {
  "use strict";
  const $ = (id) => document.getElementById(id);
  const T = window.TAXONOMY;

  const state = {
    method: "did",
    view: "paper",
    se: "cluster",
  };

  /* ---- read simulation controls ---------------------------------------- */
  function params() {
    const num = (id, dflt) => { const v = parseFloat($(id).value); return isNaN(v) ? dflt : v; };
    return {
      nClusters: Math.round(num("nClusters", 6)),
      unitsPerCluster: Math.round(num("unitsPerCluster", 8)),
      T: Math.round(num("T", 12)),
      effect: num("effect", 1.0),
      icc: num("icc", 0.2),
      preTrend: num("preTrend", 0.0),
      timing: $("timing").value,
      seed: Math.round(num("seed", 12345)),
      reps: Math.round(num("reps", 199)),
    };
  }

  /* ---- build dataset for the active methodology ------------------------ */
  function buildData(method, p) {
    if (method === "did") return DGP.did({ nClusters: p.nClusters, unitsPerCluster: p.unitsPerCluster, T: p.T, tau: p.effect, ramp: 3, preTrend: p.preTrend, icc: p.icc, timing: p.timing, seed: p.seed });
    if (method === "rct") return DGP.rct({ n: p.nClusters * p.unitsPerCluster * 6, tau: p.effect * 0.3, nClusters: Math.max(p.nClusters, 8), icc: p.icc, seed: p.seed });
    if (method === "iv") return DGP.iv({ n: Math.max(p.nClusters * p.unitsPerCluster * 40, 1600), beta: p.effect * 0.5, strength: 0.6, seed: p.seed });
    if (method === "rdd") return DGP.rdd({ n: Math.max(p.nClusters * p.unitsPerCluster * 40, 2000), jump: p.effect * 0.6, seed: p.seed });
  }

  /* ---- per-exhibit rendering ------------------------------------------- */
  function renderExhibit(methodKey, ex, data, p) {
    const rng = () => Econ.makeRNG(p.seed + 101); // fresh, reproducible RNG per call
    let inner = "", latex = null, extra = "";
    try {
      switch (ex.render) {
        case "summary": inner = Render.table1(Models.didSummary(data), state.view); break;
        case "trends": inner = Render.trends(data, state.view); break;
        case "eventstudy": {
          const es = Models.didEventStudy(data, state.se, rng(), p.reps);
          inner = Render.eventStudy(es, state.view);
          extra = badge("SE: " + (Models.METHOD_LABEL[es.method] || es.method)) + badge("true effect = " + p.effect.toFixed(2));
          break;
        }
        case "main": {
          const m = Models.didMain(data, state.se, rng(), p.reps);
          inner = Render.mainTable(m, state.view);
          latex = Render.latexMain(m);
          break;
        }
        case "robustness": {
          const rob = Models.didRobustnessSE(data, rng(), p.reps);
          inner = Render.robustnessTable(rob, state.view);
          latex = Render.latexRobustness(rob);
          extra = rob.nClusters <= 30 ? badge("⚠ few clusters (G=" + rob.nClusters + ") → trust the wild bootstrap", "warn") : "";
          break;
        }
        case "heterogeneity": inner = Render.forest(Models.didHeterogeneity(data, state.se, rng(), p.reps), state.view); break;
        case "balance": {
          const bal = Models.rctBalance(data);
          const itt = Models.rctITT(data, state.se, rng(), p.reps);
          inner = Render.balanceTable(bal, state.view) +
            "<p class='inline-est'>ITT estimate: <strong>" + itt.coef.toFixed(3) + "</strong> " +
            (itt.stars ? "<sup class='st'>" + itt.stars + "</sup> " : "") + "(SE " + (isNaN(itt.se) ? "—" : itt.se.toFixed(3)) + ", " + Models.METHOD_LABEL[state.se] + ")</p>";
          break;
        }
        case "firststage": {
          const ivf = Models.ivFit(data);
          inner = Render.firstStage(ivf, state.view);
          extra = badge("First-stage F = " + ivf.firstStage.F.toFixed(1)) + badge("2SLS β = " + ivf.tsls.toFixed(3)) + badge("true β = " + (p.effect * 0.5).toFixed(3));
          break;
        }
        case "rdplot": {
          const rd = Models.rddFit(data, 22);
          inner = Render.rdPlot(rd, state.view);
          extra = badge("estimated jump = " + rd.jump.toFixed(3)) + badge("true jump = " + rd.trueJump.toFixed(3));
          break;
        }
        default: inner = "<em>Not implemented.</em>";
      }
    } catch (err) {
      inner = "<div class='err'>Error: " + (err && err.message ? err.message : err) + "</div>";
    }
    return { inner, latex, extra };
  }

  /* ---- code recipe (maps each exhibit to the user's Stata/R skills) ----- */
  function recipe(methodKey, ex) {
    const S = {
      summary: ["estpost summarize y x, by(ever_treated)\nesttab using table1.tex, cells(\"mean sd\") booktabs", "datasummary(y + x ~ Mean + SD, data = df, output = \"table1.tex\")"],
      trends: ["collapse (mean) y, by(period ever_treated)\ntwoway (line y period if ever_treated) (line y period if !ever_treated)", "df |> group_by(period, ever_treated) |> summarise(y=mean(y)) |> ggplot(aes(period,y,color=ever_treated))+geom_line()"],
      eventstudy: ["eventdd y i.unit i.period, timevar(evt) cluster(group)\n// or csdid / eventstudyinteract (Sun-Abraham)", "feols(y ~ i(evt, ref=-1) | unit + period, df, cluster=~group) |> iplot()"],
      main: ["eststo: reghdfe y D, absorb(unit period) cluster(group)\nesttab using main.tex, se star(* .10 ** .05 *** .01) booktabs", "feols(y ~ D | unit + period, df, cluster=~group) |> modelsummary(output=\"main.tex\")"],
      robustness: ["reghdfe y D, absorb(unit period) cluster(group)\nboottest D, reps(999) cluster(group)   // wild bootstrap", "feols(y~D|unit+period, df) ; boottest::boottest(m, clustid=~group, B=999, param=\"D\")"],
      heterogeneity: ["reghdfe y c.D##i.subgroup, absorb(unit period) cluster(group)\ncoefplot, keep(*D*)", "feols(y ~ D | unit+period, split=~subgroup, df) |> coefplot()"],
      balance: ["iebaltab x1 x2 x3, grpvar(treat) save(balance) format(%9.3f)", "modelsummary(datasummary_balance(~treat, df))"],
      firststage: ["ivreg2 y (d = z), robust first   // reports first-stage F", "feols(y ~ 1 | d ~ z, df) |> fitstat(~ivf)"],
      rdplot: ["rdplot y run, c(0)\nrdrobust y run, c(0)   // CCT robust bias-corrected", "rdrobust::rdplot(df$y, df$run, c=0); rdrobust(df$y, df$run, c=0)"],
    };
    const r = S[ex.render] || ["", ""];
    return { stata: r[0], r: r[1] };
  }

  function badge(txt, cls) { return "<span class='badge " + (cls || "") + "'>" + esc(txt) + "</span>"; }
  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }

  /* ---- main render ----------------------------------------------------- */
  let pendingId = 0;
  function renderAll() {
    const p = params();
    const meth = T.methodologies[state.method];
    $("methBlurb").innerHTML = "<strong>" + meth.name + ".</strong> " + meth.blurb;
    $("methMeta").innerHTML =
      "<div><span class='ml'>Estimators</span> " + meth.estimators.join(" · ") + "</div>" +
      "<div><span class='ml'>Diagnostics</span> " + (meth.diagnostics || []).join(" · ") + "</div>" +
      "<div><span class='ml'>Reference packages</span> <a href='" + meth.ejdSearch + "' target='_blank' rel='noopener'>browse on EJD →</a></div>";

    const cards = $("cards");
    cards.innerHTML = "<div class='computing'>Computing simulations…</div>";
    const myId = ++pendingId;

    // defer heavy compute so the "computing" message paints first
    setTimeout(() => {
      if (myId !== pendingId) return;
      const data = buildData(state.method, p);
      let html = "";
      meth.exhibits.forEach((ex, i) => {
        const out = renderExhibit(state.method, ex, data, p);
        const rc = recipe(state.method, ex);
        const cid = "card-" + state.method + "-" + ex.id;
        html += "<section class='card' id='" + cid + "'>";
        html += "<header class='card-h'><div><span class='pos'>" + ex.n + "</span><h3>" + ex.title + "</h3>" +
          "<span class='typ " + ex.type + "'>" + ex.type + "</span></div></header>";
        html += "<div class='story'>" + ex.story + "</div>";
        if (out.extra) html += "<div class='badges'>" + out.extra + "</div>";
        html += "<div class='exhibit'>" + out.inner + "</div>";
        // EJD references
        if (ex.ejd && ex.ejd.length) {
          html += "<div class='ejd'>Reference packages: " +
            ex.ejd.map(e => "<a href='" + e.url + "' target='_blank' rel='noopener'>" + esc(e.label) + " →</a>").join(" · ") + "</div>";
        }
        // export / recipe
        html += "<details class='recipe'><summary>Replication recipe (Stata · R)" + (out.latex ? "  +  LaTeX export" : "") + "</summary>";
        html += "<div class='reci'><div><span class='ml'>Stata</span><pre>" + esc(rc.stata) + "</pre></div>" +
          "<div><span class='ml'>R</span><pre>" + esc(rc.r) + "</pre></div></div>";
        if (out.latex) {
          html += "<div class='ltx'><div class='ltx-h'><span class='ml'>LaTeX (booktabs / threeparttable)</span>" +
            "<button class='copy' data-target='ltx-" + cid + "'>Copy LaTeX</button></div>" +
            "<pre id='ltx-" + cid + "'>" + esc(out.latex) + "</pre></div>";
        }
        html += "</details>";
        html += "</section>";
      });
      if (myId !== pendingId) return;
      cards.innerHTML = html;
      bindCopy();
    }, 30);
  }

  function bindCopy() {
    document.querySelectorAll(".copy").forEach(btn => {
      btn.addEventListener("click", () => {
        const pre = document.getElementById(btn.dataset.target);
        const txt = pre.textContent;
        navigator.clipboard && navigator.clipboard.writeText(txt).then(
          () => { btn.textContent = "Copied ✓"; setTimeout(() => (btn.textContent = "Copy LaTeX"), 1500); },
          () => selectText(pre)
        ) || selectText(pre);
      });
    });
  }
  function selectText(el) { const r = document.createRange(); r.selectNodeContents(el); const s = getSelection(); s.removeAllRanges(); s.addRange(r); }

  /* ---- controls wiring ------------------------------------------------- */
  function init() {
    // methodology tabs
    document.querySelectorAll("[data-method]").forEach(b => b.addEventListener("click", () => {
      state.method = b.dataset.method;
      document.querySelectorAll("[data-method]").forEach(x => x.classList.toggle("active", x === b));
      document.body.classList.toggle("method-did", state.method === "did");
      renderAll();
    }));
    // view toggle
    document.querySelectorAll("[data-view]").forEach(b => b.addEventListener("click", () => {
      state.view = b.dataset.view;
      document.querySelectorAll("[data-view]").forEach(x => x.classList.toggle("active", x === b));
      document.body.classList.toggle("present", state.view === "present");
      renderAll();
    }));
    // SE selector
    $("seMethod").addEventListener("change", () => { state.se = $("seMethod").value; renderAll(); });
    // sim controls (debounced)
    let t;
    ["nClusters", "unitsPerCluster", "T", "effect", "icc", "preTrend", "timing", "seed", "reps"].forEach(id => {
      $(id).addEventListener("input", () => { clearTimeout(t); t = setTimeout(renderAll, 220); });
      $(id).addEventListener("change", () => { clearTimeout(t); renderAll(); });
    });
    // range value mirrors
    document.querySelectorAll("input[type=range]").forEach(r => {
      const out = document.querySelector("[data-out='" + r.id + "']");
      const upd = () => { if (out) out.textContent = r.value; };
      r.addEventListener("input", upd); upd();
    });
    $("reseed").addEventListener("click", () => { $("seed").value = (parseInt($("seed").value, 10) || 0) + 1; renderAll(); });
    renderAll();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
