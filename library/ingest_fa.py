#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ingest_fa.py — Auto-ingest a set of papers (from Zotero PDFs) into the library.

For each PDF: detect Table/Figure exhibits + their pages, capture them as
cropped PDF + PNG, infer methodology tags from the caption + paper text, write
a (templated) good-practice note, and link to the replication package via DOI /
EObrowse on EJD. Appends to library/library.json (keeps katrina + kapor).

Inline code is only attached where a local replication package is available
(AEA packages live behind openICPSR login, so for those we link the package).
"""
import json, os, re, subprocess, sys

LAB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPERS = os.path.join(LAB, "library", "papers")
LIBJSON = os.path.join(LAB, "library", "library.json")
DPI = 120
EJD = "https://ejd.econ.mathematik.uni-ulm.de/?search="

def run(cmd): return subprocess.run(cmd, capture_output=True, text=True)
def pages_of(pdf):
    r = run(["pdfinfo", pdf]); m = re.search(r"Pages:\s+(\d+)", r.stdout); return int(m.group(1)) if m else 0
def ptext(pdf, p):
    return run(["pdftotext", "-f", str(p), "-l", str(p), "-layout", pdf, "-"]).stdout

# ---- methodology vocabulary (must be a subset of library TAGS) -------------
KW = [
    ("Regression discontinuity", r"regression discontinuit|discontinuit|running variable|bandwidth|\bcutoff\b|\bRD\b|forcing variable|rdrobust|McCrary"),
    ("IV / 2SLS", r"instrument|2sls|two-stage least|first stage|\bIV\b|exclusion restriction"),
    ("Difference-in-differences", r"difference-in-difference|diff-in-diff|\bDiD\b|event[- ]study|parallel trend|two-way fixed effect|staggered"),
    ("RCT / experiment", r"randomi[sz]ed|\bRCT\b|experiment|treatment group|control group|intent[- ]to[- ]treat|\bITT\b|lottery"),
    ("Structural estimation", r"structural|estimate the model|utility function|demand model|equilibrium|counterfactual|willingness to pay"),
    ("Summary statistics", r"summary statistic|descriptive statistic|sample characteristic|means and standard"),
    ("Balance / first stage", r"balance|covariate balance|orthogonal"),
    ("Robustness / placebo", r"robustness|placebo|sensitivity|specification check|falsification"),
    ("Heterogeneity", r"heterogen|subgroup|by (gender|race|income|quartile)|interaction"),
    ("Event study", r"event[- ]study|dynamic effect|leads and lags"),
]
def tags_from(text):
    t = text.lower(); out = []
    for label, pat in KW:
        if re.search(pat, t, re.I): out.append(label)
    return out

# ---- exhibit caption detector ---------------------------------------------
# matches caption headers across AEA / Econometrica / Elsevier / REStat styles
CAP = re.compile(r"^\s{0,12}(TABLE|Table|FIGURE|Figure|Fig\.?)\s+(A?\d{1,2}|[IVXL]{1,5})\b\s*([.:—–-]|\b)", )
BADPREFIX = re.compile(r"(see|in|from|of|and|panel|appendix\s+figure|appendix\s+table|columns?|rows?)\s*$", re.I)

FUNCWORD = re.compile(r"^(the|a|an|is|are|was|were|includes?|reports?|shows?|plots?|presents?|displays?|gives?|lists?|provides?|for|in|of|to|and|but|that|this|these|those|we|it|note|notes|panel|column|row|above|below|see)\b", re.I)

def detect_exhibits(pdf):
    np = pages_of(pdf)
    found = {}; order = []
    for p in range(1, np + 1):
        lines = ptext(pdf, p).splitlines()
        for i, line in enumerate(lines):
            m = CAP.match(line)
            if not m: continue
            kind = "table" if m.group(1).lower().startswith("t") else "figure"
            num = m.group(2)
            pre = line[:m.start()].strip()
            if BADPREFIX.search(pre): continue
            if pre and pre[-1:].islower(): continue            # mid-sentence mention
            after = line[m.end():].strip(" .—–-:")
            # AEA style: title may sit on the next non-empty line after "Table N—"
            if len(after) < 4:
                for nx in lines[i + 1:i + 3]:
                    if nx.strip(): after = nx.strip(" .—–-:"); break
            # a real caption's title starts with a capital and is not an in-text verb phrase
            if not after or not after[0].isupper(): continue
            if FUNCWORD.match(after): continue
            if re.sub(r"[^a-z]", "", after.lower()) in ("table", "figure", "tablenum", ""): continue
            sec = "appendix" if num.upper().startswith("A") else "main"
            key = (kind, num.upper())
            if key not in found:
                title = re.sub(r"\s+", " ", after)[:120]
                found[key] = dict(kind=kind, num=num, page=p, title=title, section=sec)
                order.append(key)
    return [found[k] for k in order]

def evaluation(ex, methods):
    kind, num, title = ex["kind"], ex["num"], ex["title"].lower()
    meth = ", ".join(methods[:2]) if methods else "the paper's design"
    if re.search(r"summary|descriptive|sample|characteristic", title) or (num in ("1", "I") and kind == "table"):
        return ["Orienting descriptive exhibit: establishes the sample and key variables before identification — the reader sizes the setting first."]
    if kind == "figure" and re.search(r"discontinuit|cutoff|bandwidth|rd\b", " ".join(methods).lower() + " " + title):
        return ["Identification figure (RD plot): binned means around the cutoff make the discontinuity visible — the emotional core of an RD paper. Check that bins and bandwidth aren't over-smoothed."]
    if kind == "figure" and re.search(r"event|dynamic|trend", title + " " + " ".join(methods).lower()):
        return ["Event-study figure: leads test parallel trends, lags trace the dynamic effect. Confidence intervals and a marked reference period are the things to verify."]
    if re.search(r"robustness|placebo|sensitivity|specification", title):
        return ["Robustness/placebo exhibit defending the design — shown for completeness; the value is that the estimate survives alternative choices."]
    if re.search(r"balance", title):
        return ["Balance exhibit: checks that treatment/comparison groups are comparable on baseline covariates before causal claims."]
    if re.search(r"heterogen|subgroup|by ", title):
        return [f"Heterogeneity exhibit: effects across subgroups. A forest/coefficient plot reads faster than a wide interaction table."]
    if kind == "table":
        return [f"Results table ({meth}): look for progressive specification columns, the outcome mean, and a clearly flagged preferred estimate."]
    return [f"Exhibit supporting the {meth} analysis."]

def capture(pdf, page, dest_base):
    tmp = dest_base + "._raw.pdf"
    run(["pdfseparate", "-f", str(page), "-l", str(page), pdf, tmp])
    cropped = dest_base + ".pdf"
    run(["pdfcrop", "--margins", "6", tmp, cropped])
    if not os.path.exists(cropped):
        if os.path.exists(tmp): os.rename(tmp, cropped)
        else: return False
    run(["pdftoppm", "-png", "-r", str(DPI), "-singlefile", cropped, dest_base])
    if os.path.exists(tmp): os.remove(tmp)
    return os.path.exists(cropped)

def slug(author, year, title):
    s = re.sub(r"[^a-z0-9]+", "-", (author + "-" + (year or "") + "-" + title[:24]).lower()).strip("-")
    return s[:48]

def ingest(meta, max_exhibits=24):
    pdf = meta["pdf"]
    pid = slug(meta["author"], meta.get("year", ""), meta["title"])
    dst = os.path.join(PAPERS, pid, "exhibits")
    exhibits = detect_exhibits(pdf)
    if not exhibits or len(exhibits) > max_exhibits + 12:
        # too few or suspiciously many (detector noise) — cap to first max
        exhibits = exhibits[:max_exhibits]
    methods = sorted(set(tags_from(meta["title"]) + tags_from(ptext(pdf, 1) + ptext(pdf, 2))))
    doi = meta.get("doi", "")
    pkg = f"https://doi.org/{doi}" if doi else ""
    ejdq = EJD + meta["title"][:40].replace(" ", "+")
    recs = []
    for ex in exhibits[:max_exhibits]:
        sec = ex["section"]
        base = os.path.join(dst, sec, f"{ex['kind']}-{ex['num']}")
        os.makedirs(os.path.dirname(base), exist_ok=True)
        if not capture(pdf, ex["page"], base): continue
        extags = sorted(set(methods + tags_from(ex["title"]))) or methods or ["Results"]
        recs.append({
            "paper": pid, "section": sec, "num": f"{'Table' if ex['kind']=='table' else 'Figure'} {ex['num']}",
            "kind": ex["kind"], "title": ex["title"] or f"{ex['kind'].title()} {ex['num']}",
            "tags": extags[:4], "programs": [], "packageUrl": pkg, "ejdUrl": ejdq,
            "lang": "—", "exports": False, "snippet": "",
            "evaluation": evaluation(ex, methods), "noteMap": "",
            "png": f"library/papers/{pid}/exhibits/{sec}/{ex['kind']}-{ex['num']}.png",
            "pdf": f"library/papers/{pid}/exhibits/{sec}/{ex['kind']}-{ex['num']}.pdf",
        })
    pmeta = {
        "id": pid, "short": f"{meta['author']} ({meta.get('year','')})", "title": meta["title"],
        "journal": meta.get("journal", ""), "year": meta.get("year", ""), "doi": doi,
        "methods": methods, "lang": "—",
        "dataNote": "Exhibits captured from the published PDF. Replication package linked via DOI / EJD.",
        "site": "", "package": pkg, "n_exhibits": len(recs),
    }
    print(f"  {pid}: {len(recs)} exhibits  [{', '.join(methods[:3])}]")
    return pmeta, recs

def main():
    fa = json.load(open("/tmp/fa_papers.json"))
    SKIP_DOIS = {"10.3982/ECTA18385"}  # Kapor already hand-curated
    lib = json.load(open(LIBJSON))
    existing_ids = {p["id"] for p in lib["papers"]}
    for r in fa:
        if r.get("doi") in SKIP_DOIS or not r.get("pdf"): continue
        try:
            pmeta, recs = ingest(r)
        except Exception as e:
            print(f"  FAIL {r['title'][:40]}: {e}"); continue
        if not recs: continue
        if pmeta["id"] in existing_ids:  # rebuild: replace
            lib["papers"] = [p for p in lib["papers"] if p["id"] != pmeta["id"]]
            lib["exhibits"] = [e for e in lib["exhibits"] if e["paper"] != pmeta["id"]]
        lib["papers"].append(pmeta); lib["exhibits"].extend(recs)
        existing_ids.add(pmeta["id"])
    # ensure FA methodology tags are in the global tag list
    for label, _ in KW:
        if label not in lib["tags"]: lib["tags"].append(label)
    json.dump(lib, open(LIBJSON, "w"), indent=1, ensure_ascii=False)
    print(f"\nLibrary now: {len(lib['papers'])} papers, {len(lib['exhibits'])} exhibits")

if __name__ == "__main__":
    main()
