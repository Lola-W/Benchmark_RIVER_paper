#!/usr/bin/env python3
"""
Batch-fit ZIP vs ZINB (NB2 / Gamma-Poisson) to coverage hist TSVs (depth,count).

Outputs:
  - fig1_per_tsv/<sample>.zip_vs_zinb.png   (per TSV: full depth range, observed + fitted pmf)
  - coverage_vs_alpha_colored_named.png     (all TSVs: coverage vs alpha, colored by family, labeled by sample)
  - all_fits.zip_zinb.tsv                   (combined table)

Headless safe: forces matplotlib Agg backend (no DISPLAY needed).
"""

import os
import re
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------
# Only the TSVs you specified
# ----------------------------
TSV_LIST = [
    "hidef_rerun_coverage_hist.merged.tsv",
    "illumina_100x_coverage_hist.merged.tsv",
    "illumina_200x_coverage_hist.merged.tsv",
    "illumina_30x_coverage_hist.merged.tsv",
    "illumina_6566_blood_coverage_hist.merged.tsv",
    "illumina_6566_sperm_coverage_hist.merged.tsv",
    "illumina_full_coverage_hist.merged.tsv",
    "nanoseq_clean_6566_coverage_hist.merged.tsv",
    "nanoseq_clean_7614_coverage_hist.merged.tsv",
    "ppmseq_100x_coverage_hist.merged.tsv",
    "ppmseq_200x_coverage_hist.merged.tsv",
    "ppmseq_30x_coverage_hist.merged.tsv",
    "ppmseq_6566_coverage_hist.merged.tsv",
    "ppmseq_7614_coverage_hist.merged.tsv",
    "PTA_01_coverage_hist.merged.tsv",
    "PTA_02_coverage_hist.merged.tsv",
    "PTA_03_coverage_hist.merged.tsv",
    "PTA_04_coverage_hist.merged.tsv",
    "udseq_clean_30x_coverage_hist.merged.tsv",
    "udseq_clean_6566_coverage_hist.merged.tsv",
    "udseq_clean_7614_coverage_hist.merged.tsv",
]

FIG1_DIR = "fig1_per_tsv"
OUT_TABLE = "all_fits.zip_zinb.tsv"
OUT_FIG2 = "coverage_vs_alpha_colored_named.png"


# ----------------------------
# Family mapping + colors
# ----------------------------
FAMILY_COLORS = {
    "illumina": "tab:blue",
    "nanoseq": "tab:orange",
    "ppmseq": "tab:green",
    "udseq": "tab:red",
    "pta": "tab:purple",
    "hidef": "black",
    "other": "gray",
}

def get_family(sample: str) -> str:
    s = sample.lower()
    if s.startswith("illumina"): return "illumina"
    if s.startswith("nanoseq"): return "nanoseq"
    if s.startswith("ppmseq"): return "ppmseq"
    if s.startswith("udseq"): return "udseq"
    if s.startswith("pta"): return "pta"
    if s.startswith("hidef"): return "hidef"
    return "other"

def get_sample_name(fn: str) -> str:
    # everything before "_coverage_hist"
    return re.sub(r"_coverage_hist.*", "", fn, flags=re.IGNORECASE)


# ----------------------------
# Distributions
# ----------------------------
def poisson_logpmf(k, lam):
    return k * np.log(lam) - lam - gammaln(k + 1.0)

def nb2_logpmf(k, mu, alpha):
    # NB2: Var = mu + alpha*mu^2
    size = 1.0 / alpha
    prob = size / (size + mu)
    return (
        gammaln(k + size)
        - gammaln(size)
        - gammaln(k + 1.0)
        + size * np.log(prob)
        + k * np.log(1.0 - prob)
    )

def zip_loglik(x, k, n, eps=1e-12):
    t_pi, t_lam = x
    pi = 1.0 / (1.0 + np.exp(-t_pi))
    lam = np.exp(t_lam)

    logp = poisson_logpmf(k, lam)
    p = np.where(k == 0,
                 pi + (1.0 - pi) * np.exp(logp),
                 (1.0 - pi) * np.exp(logp))
    p = np.clip(p, eps, 1.0)
    return np.sum(n * np.log(p))

def zinb_loglik(x, k, n, eps=1e-12):
    t_pi, t_mu, t_alpha = x
    pi = 1.0 / (1.0 + np.exp(-t_pi))
    mu = np.exp(t_mu)
    alpha = np.exp(t_alpha)

    logp = nb2_logpmf(k, mu, alpha)
    p = np.where(k == 0,
                 pi + (1.0 - pi) * np.exp(logp),
                 (1.0 - pi) * np.exp(logp))
    p = np.clip(p, eps, 1.0)
    return np.sum(n * np.log(p))

def fit_mle(loglik_fn, x0, k, n):
    res = minimize(lambda x: -loglik_fn(x, k, n), x0, method="L-BFGS-B")
    return res

def zip_pmf(k, pi, lam):
    k = np.asarray(k)
    p_pois = np.exp(poisson_logpmf(k, lam))
    p = (1.0 - pi) * p_pois
    p[k == 0] = pi + (1.0 - pi) * np.exp(poisson_logpmf(np.array([0]), lam))[0]
    return p

def zinb_pmf(k, pi, mu, alpha):
    k = np.asarray(k)
    p_nb = np.exp(nb2_logpmf(k, mu, alpha))
    p = (1.0 - pi) * p_nb
    p[k == 0] = pi + (1.0 - pi) * np.exp(nb2_logpmf(np.array([0]), mu, alpha))[0]
    return p


# ----------------------------
# Main
# ----------------------------
def main():
    os.makedirs(FIG1_DIR, exist_ok=True)
    rows = []

    for fn in TSV_LIST:
        if not os.path.exists(fn):
            print(f"[SKIP missing] {fn}")
            continue
        if fn.endswith(".fit_params.tsv"):
            print(f"[SKIP params file] {fn}")
            continue

        df = pd.read_csv(fn, sep="\t")
        df.columns = [c.strip().lower() for c in df.columns]
        if "depth" not in df.columns or "count" not in df.columns:
            print(f"[SKIP bad format] {fn} (need depth,count)")
            continue

        df = df[["depth", "count"]].copy()
        df["depth"] = df["depth"].astype(int)
        df["count"] = df["count"].astype(np.int64)
        df = df.sort_values("depth")

        k = df["depth"].to_numpy(dtype=int)
        n = df["count"].to_numpy(dtype=np.int64)

        total = n.sum()
        zero = n[k == 0].sum() if np.any(k == 0) else 0
        nonzero = total - zero
        coverage = nonzero / total if total > 0 else np.nan

        mean = np.sum(k * n) / total
        var = np.sum(((k - mean) ** 2) * n) / total

        frac0 = zero / total
        lam0 = max(mean / max(1e-6, (1.0 - frac0)), 1e-6)

        # ZIP init + fit
        x0_zip = np.array([np.log(frac0 / max(1e-12, 1.0 - frac0)), np.log(lam0)], dtype=float)
        res_zip = fit_mle(zip_loglik, x0_zip, k, n)
        pi_zip = 1.0 / (1.0 + np.exp(-res_zip.x[0]))
        lam_zip = np.exp(res_zip.x[1])
        ll_zip = -res_zip.fun

        # ZINB init + fit
        mu0 = lam0
        alpha0 = max((var - mu0) / (mu0 * mu0 + 1e-12), 1e-6)
        x0_zinb = np.array([
            np.log(frac0 / max(1e-12, 1.0 - frac0)),
            np.log(mu0),
            np.log(alpha0),
        ], dtype=float)
        res_zinb = fit_mle(zinb_loglik, x0_zinb, k, n)
        pi_zinb = 1.0 / (1.0 + np.exp(-res_zinb.x[0]))
        mu_zinb = np.exp(res_zinb.x[1])
        alpha_zinb = np.exp(res_zinb.x[2])
        ll_zinb = -res_zinb.fun

        # AIC/BIC
        zip_k = 2
        zinb_k = 3
        zip_aic = 2 * zip_k - 2 * ll_zip
        zinb_aic = 2 * zinb_k - 2 * ll_zinb
        zip_bic = zip_k * np.log(total) - 2 * ll_zip
        zinb_bic = zinb_k * np.log(total) - 2 * ll_zinb

        sample = get_sample_name(fn)
        family = get_family(sample)

        rows.append({
            "tsv": fn,
            "sample": sample,
            "family": family,
            "total_loci": int(total),
            "zero": int(zero),
            "nonzero": int(nonzero),
            "coverage_pct": coverage * 100.0,
            "ZIP_pi": pi_zip,
            "ZIP_lambda": lam_zip,
            "ZIP_loglik": ll_zip,
            "ZIP_AIC": zip_aic,
            "ZIP_BIC": zip_bic,
            "ZIP_converged": bool(res_zip.success),
            "ZINB_pi": pi_zinb,
            "ZINB_mu": mu_zinb,
            "ZINB_alpha": alpha_zinb,
            "ZINB_loglik": ll_zinb,
            "ZINB_AIC": zinb_aic,
            "ZINB_BIC": zinb_bic,
            "ZINB_converged": bool(res_zinb.success),
        })

        # -------- Figure 1 (full depth range) --------
        frac_obs = n / n.sum()

        p_zip = zip_pmf(k, pi_zip, lam_zip)
        p_zinb = zinb_pmf(k, pi_zinb, mu_zinb, alpha_zinb)

        # Renormalize to the observed support (k present in file)
        p_zip = p_zip / p_zip.sum()
        p_zinb = p_zinb / p_zinb.sum()

        plt.figure(figsize=(10, 4))
        plt.bar(k, frac_obs, width=0.9, alpha=0.4, label="observed")
        plt.plot(k, p_zip, marker="o", linewidth=1.1,
                 label=f"ZIP pi={pi_zip:.3f}, lambda={lam_zip:.3f}")
        plt.plot(k, p_zinb, marker="o", linewidth=1.1,
                 label=f"ZINB pi={pi_zinb:.3f}, mu={mu_zinb:.3f}, alpha={alpha_zinb:.3f}")
        plt.xlabel("Depth")
        plt.ylabel("Fraction of genomic loci")
        plt.title(f"{sample} (coverage={coverage*100:.2f}%)")
        plt.legend()
        plt.tight_layout()
        out_fig1 = os.path.join(FIG1_DIR, f"{sample}.zip_vs_zinb.png")
        plt.savefig(out_fig1, dpi=250)
        plt.close()

        print(f"[OK] {sample:25s} cov={coverage*100:6.2f}%  alpha={alpha_zinb:.4g}  -> {out_fig1}")

    # Save combined table
    if len(rows) == 0:
        print("[DONE] no TSV processed.")
        return

    res_df = pd.DataFrame(rows)
    res_df.to_csv(OUT_TABLE, sep="\t", index=False)

    # -------- Figure 2 (coverage vs alpha) --------
    plt.figure(figsize=(7, 6))
    for _, r in res_df.iterrows():
        c = FAMILY_COLORS.get(r["family"], "gray")
        plt.scatter(r["coverage_pct"], r["ZINB_alpha"], color=c)
        plt.text(r["coverage_pct"], r["ZINB_alpha"], r["sample"], fontsize=7)

    plt.xlabel("Theoretical genomic coverage (%) = non-zero / (zero + non-zero)")
    plt.ylabel("ZINB overdispersion alpha")
    plt.title("Coverage vs ZINB alpha (colored by sequencing family)")
    plt.tight_layout()
    plt.savefig(OUT_FIG2, dpi=300)
    plt.close()

    print("\n[DONE]")
    print(f"Table: {OUT_TABLE}")
    print(f"Fig2 : {OUT_FIG2}")
    print(f"Fig1 : {FIG1_DIR}/")


if __name__ == "__main__":
    main()
