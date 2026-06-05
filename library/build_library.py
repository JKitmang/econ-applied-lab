#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_library.py — Build the multi-paper Exhibit Library for econ-applied-lab.

A researcher browses tables/figures across papers, filters by methodology, and
gets the exact code that produces each exhibit (to replicate / adapt).

Each paper is captured from its published PDF (exhibits preserve paper format)
and paired with the package code that produces it. No proprietary data is
redistributed. Aggregates everything into library/library.json.

Sources:
  * Katrina  — imported from the already-built replication-katrina catalog.
  * Kapor    — captured from the paper PDF in the user's Zotero; code from the
               local replication package.
"""
import json, os, re, shutil, subprocess

LAB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # econ-applied-lab repo root
PAPERS = os.path.join(LAB, "library", "papers")
DPI = 120

# ---- global methodology vocabulary (filter chips, ordered) -----------------
TAGS = [
    "Summary statistics", "Descriptive figure (map)", "OLS + fixed effects",
    "Reduced form", "IV / 2SLS", "Quartile heterogeneity",
    "Linear-in-means peer effects", "Model specification tests",
    "Placebo / robustness", "School-level aggregate", "Non-test outcomes",
    "Structural estimation", "Counterfactual simulation", "Model fit",
    "Welfare / diversity",
]

def run(cmd): return subprocess.run(cmd, capture_output=True, text=True)

def capture(pdf, pages, dest_base):
    tmp = dest_base + "._raw.pdf"
    if len(pages) == 1:
        run(["pdfseparate", "-f", str(pages[0]), "-l", str(pages[0]), pdf, tmp])
    else:
        parts = []
        for i, p in enumerate(pages):
            pp = f"{dest_base}._p{i}.pdf"; run(["pdfseparate", "-f", str(p), "-l", str(p), pdf, pp]); parts.append(pp)
        run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", f"-sOutputFile={tmp}"] + parts)
        for pp in parts:
            if os.path.exists(pp): os.remove(pp)
    cropped = dest_base + ".pdf"
    run(["pdfcrop", "--margins", "6", tmp, cropped])
    if not os.path.exists(cropped): shutil.copy(tmp, cropped)
    run(["pdftoppm", "-png", "-r", str(DPI), "-singlefile", cropped, dest_base])
    if os.path.exists(tmp): os.remove(tmp)
    return os.path.exists(cropped)

def snippet(path, patterns, n=10):
    if not os.path.exists(path): return ""
    pat = re.compile(patterns, re.I); keep = []
    with open(path, errors="ignore") as f:
        for line in f:
            if pat.search(line):
                keep.append(line.rstrip()[:160])
                if len(keep) >= n: break
    return "\n".join(keep)

# ===========================================================================
#  PAPER 1 — Katrina (import from the existing replication-katrina catalog)
# ===========================================================================
KAT_SRC = "/Users/jostin.kitmang/Downloads/Imberman Replication"
KAT_CODE_URL = "https://github.com/JKitmang/replication-katrina/blob/main/programs/"

def import_katrina():
    cat = json.load(open(os.path.join(KAT_SRC, "exhibits", "catalog.json")))
    labels = cat["tagLabels"]
    pid = "katrina"
    dst = os.path.join(PAPERS, pid, "exhibits")
    os.makedirs(os.path.join(dst, "main"), exist_ok=True)
    os.makedirs(os.path.join(dst, "appendix"), exist_ok=True)
    out = []
    for e in cat["exhibits"]:
        # copy assets
        for k in ("png", "pdf"):
            src = os.path.join(KAT_SRC, e[k])
            rel = e[k].replace("exhibits/", "")
            dest = os.path.join(dst, rel)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            if os.path.exists(src): shutil.copy(src, dest)
        out.append({
            "paper": pid, "section": e["section"], "num": e["num"], "kind": e["kind"],
            "title": e["title"], "tags": [labels[t] for t in e["tags"]],
            "programs": [{"name": p, "url": KAT_CODE_URL + p} for p in e["programs"]],
            "lang": "Stata", "exports": e["exports"], "snippet": e["snippet"],
            "evaluation": e["evaluation"], "noteMap": e.get("noteMap", ""),
            "png": f"library/papers/{pid}/exhibits/{e['png'].replace('exhibits/','')}",
            "pdf": f"library/papers/{pid}/exhibits/{e['pdf'].replace('exhibits/','')}",
        })
    meta = {
        "id": pid, "short": "Imberman, Kugler & Sacerdote (2012)",
        "title": "Katrina's Children: Evidence on the Structure of Peer Effects from Hurricane Evacuees",
        "journal": "American Economic Review", "year": 2012, "doi": "10.1257/aer.102.5.2048",
        "methods": ["Reduced form", "IV / 2SLS", "OLS + fixed effects", "Quartile heterogeneity", "Linear-in-means peer effects"],
        "lang": "Stata", "dataNote": "Micro-data are proprietary (HISD + Louisiana DOE); exhibits captured from the published PDF.",
        "site": "https://jkitmang.github.io/replication-katrina/website/exhibits.html",
        "n_exhibits": len(out),
    }
    print(f"  katrina: imported {len(out)} exhibits")
    return meta, out

# ===========================================================================
#  PAPER 2 — Kapor, "Transparency and Percent Plans" (Econometrica 2025)
# ===========================================================================
KAP_PDF = "/Users/jostin.kitmang/Zotero/storage/LJI9X39V/Kapor - 2025 - Transparency and Percent Plans.pdf"
KAP_PKG = "/Users/jostin.kitmang/Downloads/replication_package"

KAP_EXHIBITS = [
    dict(id="t1", num="Table I", kind="table", pages=[25], section="main",
         title="Estimated parameters (preferences and information)",
         tags=["Structural estimation"], code=["code/analysis/describeParameters2024.jl", "code/analysis/estimate_2024.jl"],
         eval=["A structural-parameters table: the estimated primitives (preferences, information, costs) that drive every counterfactual — the foundation the rest of the paper rests on.",
               "Reporting structural estimates with standard errors (bootstrap) is the analog of a reduced-form coefficient table; transparency about identification lives in the surrounding text.",
               "Improvement: a one-line 'what identifies this' annotation per parameter block helps readers who skip the methods section."]),
    dict(id="t2", num="Table II", kind="table", pages=[24], section="main",
         title="Admissions and Outcome Parameters",
         tags=["Structural estimation"], code=["code/analysis/fitOutcomeModels2024.jl", "code/analysis/estimate_2024.jl"],
         eval=["Separates the admissions/outcome side of the model from preferences (Table I), keeping each estimation block legible.",
               "Good practice: outcome-model parameters are shown before they are used to simulate counterfactuals, so the reader can audit the inputs.",
               "Improvement: flag which parameters are estimated vs calibrated to avoid over-crediting the estimation."]),
    dict(id="t3", num="Table III", kind="table", pages=[28], section="main",
         title="Main Results",
         tags=["Counterfactual simulation", "Welfare / diversity"], code=["code/analysis/counterfactuals_2024.jl", "code/analysis/describeResults2024.jl"],
         eval=["The headline table: counterfactual outcomes (transparency vs opaque percent plans) — the paper's actual contribution in one exhibit.",
               "Columns as policy regimes is the structural analog of 'progressive specification columns' — the reader compares worlds, not controls.",
               "Reports bootstrap uncertainty on counterfactuals, which is easy to omit and important to include.",
               "Improvement: pairing each counterfactual column with a welfare/diversity summary row makes the trade-off explicit at a glance."]),
    dict(id="f1", num="Figure 1", kind="figure", pages=[23], section="main",
         title="Model Fit, Flagship Universities",
         tags=["Model fit"], code=["code/analysis/describeFit2024.jl"],
         eval=["A model-fit figure is the structural paper's identification figure: it shows the estimated model reproduces key moments (enrollment by type) before any counterfactual is trusted.",
               "Overlaying model vs data is exactly the right visual contract with the reader.",
               "Improvement: shading the moments that were targeted vs untargeted distinguishes fit from out-of-sample validation."]),
    dict(id="f2", num="Figure 2", kind="figure", pages=[25], section="main",
         title="Costs, Information, and Admission Parameters",
         tags=["Structural estimation"], code=["code/analysis/describeParameters2024.jl"],
         eval=["Visualizes estimated parameter gradients (costs/information by group) — far more legible than a dense parameter table for heterogeneous primitives.",
               "Panels keep distinct primitives separate while sharing a scale.",
               "Improvement: confidence bands on each curve would carry the uncertainty the companion table reports."]),
    dict(id="f3", num="Figure 3", kind="figure", pages=[26], section="main",
         title="Preferences and Awareness",
         tags=["Structural estimation"], code=["code/analysis/describeResults2024.jl"],
         eval=["Turns latent preference/awareness estimates into average choice probabilities — translating structure into quantities readers understand.",
               "Good storytelling: it answers 'who knows about the policy and who wants what' visually.",
               "Improvement: a brief axis note on units (probabilities vs utils) avoids misreading."]),
    dict(id="f4", num="Figure 4", kind="figure", pages=[31], section="main",
         title="Demographics and Diversity",
         tags=["Counterfactual simulation", "Welfare / diversity"], code=["code/analysis/counterfactuals_2024.jl", "code/analysis/describeResults2024.jl"],
         eval=["The payoff figure: decomposes the policy's impact on diversity by demographic group — the equity question the paper exists to answer.",
               "A decomposition figure beats a table here because the composition shift is the message.",
               "Improvement: anchoring to the status-quo baseline within the figure makes 'how much changed' unambiguous."]),
]

def build_kapor():
    pid = "kapor-percent-plans"
    dst = os.path.join(PAPERS, pid, "exhibits", "main")
    os.makedirs(dst, exist_ok=True)
    code_dst = os.path.join(PAPERS, pid, "code")
    os.makedirs(code_dst, exist_ok=True)
    out = []
    # copy referenced code files into the library (a few small source files)
    copied = set()
    for ex in KAP_EXHIBITS:
        for c in ex["code"]:
            src = os.path.join(KAP_PKG, c)
            if os.path.exists(src) and c not in copied:
                shutil.copy(src, os.path.join(code_dst, os.path.basename(c))); copied.add(c)
    for ex in KAP_EXHIBITS:
        base = os.path.join(dst, ex["id"])
        ok = capture(KAP_PDF, ex["pages"], base)
        snip_src = os.path.join(KAP_PKG, ex["code"][0])
        snip = snippet(snip_src, r"(savefig|plot\(|Plots\.|@df|StatsPlots|CSV\.write|open\(|println|latex|writedlm|DataFrame\()")
        progs = [{"name": os.path.basename(c), "url": f"library/papers/{pid}/code/{os.path.basename(c)}"} for c in ex["code"]]
        print(("  ok " if ok else " FAIL "), "kapor", ex["num"], ex["pages"])
        out.append({
            "paper": pid, "section": ex["section"], "num": ex["num"], "kind": ex["kind"],
            "title": ex["title"], "tags": ex["tags"], "programs": progs, "lang": "Julia",
            "exports": True, "snippet": snip, "evaluation": ex["eval"], "noteMap": "",
            "png": f"library/papers/{pid}/exhibits/main/{ex['id']}.png",
            "pdf": f"library/papers/{pid}/exhibits/main/{ex['id']}.pdf",
        })
    meta = {
        "id": pid, "short": "Kapor (2025)",
        "title": "Transparency and Percent Plans",
        "journal": "Econometrica", "year": 2025, "doi": "",
        "methods": ["Structural estimation", "Counterfactual simulation", "Model fit", "Welfare / diversity"],
        "lang": "Julia (+ Stata setup)", "dataNote": "Restricted-use THEOP data; exhibits captured from the published PDF. Code from the replication package.",
        "site": "", "n_exhibits": len(out),
    }
    return meta, out

# ===========================================================================
def main():
    os.makedirs(PAPERS, exist_ok=True)
    papers, exhibits = [], []
    # Focus: applied papers with QUASI-EXPERIMENTAL designs (RD / DiD / IV).
    # Katrina (reduced-form / IV / peer effects) qualifies; Kapor (structural) does not.
    for builder in (import_katrina,):
        meta, exs = builder()
        papers.append(meta); exhibits.extend(exs)
    library = {
        "title": "Applied-Economics Exhibit Library",
        "blurb": "Real tables and figures from applied papers with credible research designs — quasi-experimental (RD, difference-in-differences, IV) and randomized (RCT) — organised by methodology and each paired with the code that produces it, so you can replicate the exhibit and adapt it to your context.",
        "tags": TAGS, "papers": papers, "exhibits": exhibits,
        "source": "Replication packages indexed by EJD (ejd.econ.mathematik.uni-ulm.de). Exhibits are shown for replication and methodological commentary.",
    }
    with open(os.path.join(LAB, "library", "library.json"), "w") as f:
        json.dump(library, f, indent=1, ensure_ascii=False)
    print(f"\nLibrary: {len(papers)} papers, {len(exhibits)} exhibits → library/library.json")

if __name__ == "__main__":
    main()
