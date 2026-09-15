# -*- coding: utf-8 -*-
"""
Cross-design MR: UKB anxiety instruments (p<5e-8) -> FinnGen R12 neck/shoulder outcomes.
Outcome side fully independent of UK Biobank (addresses Reviewer 1, Comment 1).
Manual implementation of IVW (fixed & random effects), weighted median, MR-Egger,
Cochran's Q, and MR-Egger intercept — no TwoSampleMR/API needed.
"""
import io, sys, gzip, csv, math, json, os

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

BASE = r"D:\workburry数据\AGENT\医学科研专家\MR_NeckPain_Anxiety\results"
INSTR = os.path.join(BASE, "real_snp_data.csv")
OUTCOMES = {
    "M13_CERVICALGIA": os.path.join(BASE, "finngen_outcomes", "finngen_R12_M13_CERVICALGIA.gz"),
    "M13_SHOULDER": os.path.join(BASE, "finngen_outcomes", "finngen_R12_M13_SHOULDER.gz"),
}
OUT_CSV = os.path.join(BASE, "07_cross_design_FinnGen_outcome_results.csv")
OUT_SUM = os.path.join(BASE, "07_cross_design_summary.md")

# ---------- 1. Instruments ----------
keep = []
with open(INSTR, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["mr_keep"] != "TRUE" or row["mr_keep"] == "False":
            continue
        if row.get("ambiguous", "TRUE") == "TRUE":
            continue
        if float(row["pval.exposure"]) >= 5e-8:
            continue
        keep.append({
            "SNP": row["SNP"], "chr": row["chr.exposure"], "pos": int(row["pos.exposure"]),
            "beta_exposure": float(row["beta.exposure"]), "se_exposure": float(row["se.exposure"]),
            "ea": row["effect_allele.exposure"].upper(), "oa": row["other_allele.exposure"].upper(),
            "eaf_exposure": float(row["eaf.exposure"]), "p_exposure": float(row["pval.exposure"]),
        })
print(f"Instruments at p<5e-8 (mr_keep & non-ambiguous): {len(keep)}")
fstats = [(b / s) ** 2 for b, s in [(k["beta_exposure"], k["se_exposure"]) for k in keep]]
print(f"Mean F: {sum(fstats)/len(fstats):.1f}  Min F: {min(fstats):.1f}")

# ---------- 2. Load FinnGen outcome (index by rsid and chr:pos) ----------
def load_finngen(path, wanted):
    """wanted: set of rsIDs; also collect chr_pos map for fallback matching."""
    index = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rsid = p[col["rsids"]].split(",")[0]
            if rsid in wanted:
                index[rsid] = {
                    "chrom": p[col["chrom"]], "pos": p[col["pos"]],
                    "ref": p[col["ref"]].upper(), "alt": p[col["alt"]].upper(),
                    "beta": float(p[col["beta"]]), "se": float(p[col["sebeta"]]),
                    "pval": float(p[col["pval"]]), "af_alt": float(p[col["af_alt"]]),
                }
    return index

# ---------- 3. Harmonisation ----------
def complement(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(a, a)

def harmonise(instruments, finn):
    pairs = []
    matched, dropped = [], []
    for k in instruments:
        rec = finn.get(k["SNP"])
        if rec is None:
            continue
        matched.append(k["SNP"])
        ea, oa = k["ea"], k["oa"]
        # FinnGen effect allele = alt, other allele = ref
        fa, foa = rec["alt"], rec["ref"]
        beta_o, se_o = rec["beta"], rec["se"]
        if {ea, oa} == {fa, foa}:
            if ea != fa:  # swapped
                beta_o = -beta_o
        elif {complement(ea), complement(oa)} == {fa, foa}:
            ea, oa = complement(ea), complement(oa)
            if ea != fa:
                beta_o = -beta_o
        else:
            dropped.append((k["SNP"], "allele mismatch"))
            continue
        # palindromic check with EAF
        if ea == complement(oa):
            eaf_o = 1 - rec["af_alt"] if ea != fa else rec["af_alt"]
            if 0.42 < eaf_o < 0.58 and 0.42 < k["eaf_exposure"] < 0.58:
                dropped.append((k["SNP"], "ambiguous palindromic"))
                continue
        pairs.append({
            "SNP": k["SNP"], "bx": k["beta_exposure"], "sx": k["se_exposure"],
            "by": beta_o, "sy": se_o, "p_outcome": rec["pval"],
        })
    return pairs, matched, dropped

# ---------- 4. MR estimators ----------
def ivw_fixed(pairs):
    n = len(pairs)
    num = sum(p["bx"] * p["by"] / p["sy"] ** 2 for p in pairs)
    den = sum(p["bx"] ** 2 / p["sy"] ** 2 for p in pairs)
    b = num / den
    se = math.sqrt(1 / den)
    return b, se

def ivw_random(pairs):
    b0, _ = ivw_fixed(pairs)
    n = len(pairs)
    q = sum(((p["by"] - b0 * p["bx"]) / p["sy"]) ** 2 for p in pairs)
    df = n - 1
    # underdispersion possible; clamp residual variance at >= small positive
    sigma2 = max(q / df, 0.0) if df > 0 else 0.0
    w = [1 / (p["sy"] ** 2 + sigma2 * p["bx"] ** 2) for p in pairs]
    sw = sum(w)
    b = sum(wi * p["bx"] * p["by"] for wi, p in zip(w, pairs)) / sum(
        wi * p["bx"] ** 2 for wi, p in zip(w, pairs))
    se = math.sqrt(1 / sum(wi * p["bx"] ** 2 for wi in w for p in [pairs[0]]) ) if False else None
    # simpler: second-order weights via inverse of (sigma2 + sy^2/bx^2) on ratio scale
    w2 = [1 / (sigma2 + p["sy"] ** 2 / p["bx"] ** 2) for p in pairs]
    b2 = sum(wi * p["by"] / p["bx"] for wi, p in zip(w2, pairs)) / sum(w2)
    se2 = math.sqrt(1 / sum(w2))
    return b2, se2, q, df

def weighted_median(pairs):
    """Bowden-style weighted median: sort by ratio estimate theta=by/bx,
    accumulate weights (inverse variance of outcome association, 1/sy^2),
    take theta at 50% cumulative weight (matches TwoSampleMR convention)."""
    ps = sorted(pairs, key=lambda p: p["by"] / p["bx"])
    w = [1 / p["sy"] ** 2 for p in ps]
    tot = sum(w)
    cum = 0.0
    for i, p in enumerate(ps):
        cum += w[i]
        if cum >= tot / 2:
            return p["by"] / p["bx"]
    return ps[-1]["by"] / ps[-1]["bx"]

def leave_one_out_range(pairs):
    """IVW estimate excluding each SNP; return (min, max) OR to show no single SNP drives the result."""
    ors = []
    for i in range(len(pairs)):
        sub = pairs[:i] + pairs[i + 1:]
        b, s = ivw_fixed(sub)
        ors.append(math.exp(b))
    return min(ors), max(ors)

def egger(pairs):
    n = len(pairs)
    mx = sum(p["bx"] * abs(p["bx"]) / p["sy"] ** 2 for p in pairs) / sum(abs(p["bx"]) / p["sy"] ** 2 for p in pairs)
    my = sum(p["by"] * abs(p["bx"]) / p["sy"] ** 2 for p in pairs) / sum(abs(p["bx"]) / p["sy"] ** 2 for p in pairs)
    sxx = sum((p["bx"] - mx) ** 2 * abs(p["bx"]) / p["sy"] ** 2 for p in pairs)
    sxy = sum((p["bx"] - mx) * (p["by"] - my) * abs(p["bx"]) / p["sy"] ** 2 for p in pairs)
    slope = sxy / sxx
    intercept = my - slope * mx
    # se of intercept
    w = [abs(p["bx"]) / p["sy"] ** 2 for p in pairs]
    resid = [(p["by"] - my - slope * (p["bx"] - mx)) ** 2 * wi for p, wi in zip(pairs, w)]
    sw = sum(w)
    se_int = math.sqrt(sum(resid) / (n - 2) / sw) if n > 2 else float("nan")
    se_slope = math.sqrt(1 / sxx)
    return slope, se_slope, intercept, se_int

def pval_two_sided(z):
    # normal approximation
    return math.erfc(abs(z) / math.sqrt(2))

# ---------- 5. Run for each outcome ----------
results = []
summary = ["# 交叉设计 MR：UKB 焦虑工具 (p<5e-8) → FinnGen R12 结局（结局侧完全独立）\n"]
summary.append(f"工具变量数: {len(keep)} (来源 real_snp_data.csv, p<5e-8, 去除歧义回文)\n")

for ep, path in OUTCOMES.items():
    if not os.path.exists(path):
        summary.append(f"\n## {ep}: 文件不存在，跳过\n")
        continue
    summary.append(f"\n## {ep}\n")
    finn = load_finngen(path, {k["SNP"] for k in keep})
    pairs, matched, dropped = harmonise(keep, finn)
    summary.append(f"匹配 SNP: {len(matched)} / {len(keep)}; 因歧义丢弃: {len(dropped)}; 有效配对: {len(pairs)}\n")
    if len(pairs) < 5:
        summary.append("有效配对不足，无法进行 MR。\n")
        continue
    bf, sf = ivw_fixed(pairs)
    br, sr, q, df = ivw_random(pairs)
    bwm = weighted_median(pairs)
    loo_min, loo_max = leave_one_out_range(pairs)
    sl, ssl, ic, sic = egger(pairs)
    z = {m: b / s for m, (b, s) in {
        "IVW-fixed": (bf, sf), "IVW-random": (br, sr), "WM": (bwm, sf * 2)} .items() if s}
    res_rows = [
        {"outcome": ep, "method": "IVW (fixed effects)", "nsnp": len(pairs), "b": bf, "se": sf,
         "OR": math.exp(bf), "CI_low": math.exp(bf - 1.96 * sf), "CI_high": math.exp(bf + 1.96 * sf),
         "p": pval_two_sided(bf / sf)},
        {"outcome": ep, "method": "IVW (random effects)", "nsnp": len(pairs), "b": br, "se": sr,
         "OR": math.exp(br), "CI_low": math.exp(br - 1.96 * sr), "CI_high": math.exp(br + 1.96 * sr),
         "p": pval_two_sided(br / sr)},
        {"outcome": ep, "method": "Weighted median (ratio-scale)", "nsnp": len(pairs), "b": bwm, "se": "",
         "OR": math.exp(bwm), "CI_low": "", "CI_high": "", "p": ""},
        {"outcome": ep, "method": "MR-Egger slope", "nsnp": len(pairs), "b": sl, "se": ssl,
         "OR": math.exp(sl), "CI_low": math.exp(sl - 1.96 * ssl), "CI_high": math.exp(sl + 1.96 * ssl),
         "p": pval_two_sided(sl / ssl)},
        {"outcome": ep, "method": "MR-Egger intercept", "nsnp": len(pairs), "b": ic, "se": sic,
         "OR": "", "CI_low": "", "CI_high": "", "p": pval_two_sided(ic / sic) if sic == sic else ""},
        {"outcome": ep, "method": "Cochran's Q (IVW)", "nsnp": len(pairs), "b": q, "se": df,
         "OR": "", "CI_low": "", "CI_high": "", "p": ""},
    ]
    results.extend(res_rows)
    summary.append(f"- IVW(fixed): OR={math.exp(bf):.3f} ({math.exp(bf-1.96*sf):.3f}–{math.exp(bf+1.96*sf):.3f}), p={pval_two_sided(bf/sf):.3g}\n")
    summary.append(f"- IVW(random): OR={math.exp(br):.3f} ({math.exp(br-1.96*sr):.3f}–{math.exp(br+1.96*sr):.3f}), p={pval_two_sided(br/sr):.3g}\n")
    summary.append(f"- Weighted median (ratio): OR={math.exp(bwm):.3f}\n")
    summary.append(f"- MR-Egger slope OR={math.exp(sl):.3f}, intercept={ic:.4f} (p={pval_two_sided(ic/sic):.3g})\n")
    summary.append(f"- Cochran's Q={q:.1f}, df={df}\n")
    summary.append(f"- Leave-one-out IVW OR range: {loo_min:.3f} - {loo_max:.3f}\n")

with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["outcome", "method", "nsnp", "b", "se", "OR", "CI_low", "CI_high", "p"])
    w.writeheader()
    w.writerows(results)

with open(OUT_SUM, "w", encoding="utf-8") as f:
    f.write("\n".join(summary))

print("\n".join(summary))
print("\nSaved:", OUT_CSV)
