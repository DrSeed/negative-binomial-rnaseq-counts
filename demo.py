#!/usr/bin/env python3
# Self-contained demo: why RNA-seq counts are negative-binomial, not Poisson.
# Generates its own synthetic count matrix, compares Poisson vs negative-binomial
# behaviour, and saves a 4-panel figure plus a per-gene summary table.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(1234)

N_GENES = 2000
N_SAMPLES = 40


def simulate_counts():
    # Give every gene its own "true" mean expression spanning several orders
    # of magnitude, the way a real RNA-seq experiment looks.
    true_mean = 10 ** RNG.uniform(-0.5, 3.0, size=N_GENES)

    # Poisson: variance is forced to equal the mean (no biological noise).
    poisson = RNG.poisson(lam=true_mean[:, None], size=(N_GENES, N_SAMPLES))

    # Negative binomial: add biological variability via a per-gene dispersion.
    # Variance = mean + dispersion * mean^2. We let dispersion drift down with
    # expression, mirroring the DESeq2/edgeR mean-dispersion trend.
    dispersion = 0.02 + 2.0 / (true_mean + 1.0)
    # numpy parameterises NB by number of successes n and success prob p, where
    # mean = n(1-p)/p and n = 1/dispersion.
    n = 1.0 / dispersion
    p = n / (n + true_mean)
    nb = RNG.negative_binomial(n[:, None], p[:, None], size=(N_GENES, N_SAMPLES))
    return true_mean, dispersion, poisson, nb


def per_gene_stats(mat):
    mean = mat.mean(axis=1)
    var = mat.var(axis=1, ddof=1)
    zeros = (mat == 0).mean(axis=1)
    return mean, var, zeros


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    true_mean, dispersion, poisson, nb = simulate_counts()
    p_mean, p_var, p_zero = per_gene_stats(poisson)
    nb_mean, nb_var, nb_zero = per_gene_stats(nb)

    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    # Panel 1: mean-variance relationship. Poisson hugs variance = mean;
    # NB sits well above it -> overdispersion.
    a = ax[0, 0]
    a.scatter(p_mean, p_var, s=6, alpha=0.3, label="Poisson", color="#4C72B0")
    a.scatter(nb_mean, nb_var, s=6, alpha=0.3, label="Negative binomial", color="#C44E52")
    lim = np.array([1, max(nb_var.max(), p_var.max())])
    a.plot(lim, lim, "k--", lw=1, label="variance = mean")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("gene mean count"); a.set_ylabel("gene variance")
    a.set_title("Mean-variance: NB overdisperses"); a.legend(fontsize=8)

    # Panel 2: count distribution for one moderately expressed gene.
    a = ax[0, 1]
    idx = int(np.argmin(np.abs(true_mean - 30)))
    bins = np.arange(0, max(poisson[idx].max(), nb[idx].max()) + 2)
    a.hist(poisson[idx], bins=bins, alpha=0.6, label="Poisson", color="#4C72B0", density=True)
    a.hist(nb[idx], bins=bins, alpha=0.6, label="Negative binomial", color="#C44E52", density=True)
    a.set_xlabel("count"); a.set_ylabel("density")
    a.set_title(f"One gene (mean~30): NB has fatter tail"); a.legend(fontsize=8)

    # Panel 3: estimated dispersion vs mean (method of moments).
    a = ax[1, 0]
    est_disp = np.clip((nb_var - nb_mean) / (nb_mean ** 2 + 1e-9), 1e-4, None)
    order = np.argsort(nb_mean)
    a.scatter(nb_mean, est_disp, s=6, alpha=0.25, color="#C44E52", label="estimated")
    a.plot(true_mean[order], dispersion[order], "k-", lw=1.5, label="true trend")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("gene mean count"); a.set_ylabel("dispersion")
    a.set_title("Dispersion shrinks with expression"); a.legend(fontsize=8)

    # Panel 4: fraction of zeros vs mean. The NB curve already explains the
    # zeros -> you rarely need a separate zero-inflation term.
    a = ax[1, 1]
    a.scatter(nb_mean, nb_zero, s=6, alpha=0.3, color="#C44E52", label="NB observed zeros")
    a.scatter(p_mean, p_zero, s=6, alpha=0.3, color="#4C72B0", label="Poisson zeros")
    grid = np.logspace(-0.5, 3, 100)
    a.plot(grid, np.exp(-grid), "k--", lw=1, label="Poisson P(0)=e^-mean")
    a.set_xscale("log")
    a.set_xlabel("gene mean count"); a.set_ylabel("fraction of zero counts")
    a.set_title("NB explains the zeros on its own"); a.legend(fontsize=8)

    fig.suptitle("Why RNA-seq counts are negative-binomial (synthetic data)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig("figures/demo.png", dpi=120)

    summary = pd.DataFrame({
        "gene": [f"gene_{i}" for i in range(N_GENES)],
        "true_mean": true_mean,
        "poisson_mean": p_mean, "poisson_var": p_var, "poisson_frac_zero": p_zero,
        "nb_mean": nb_mean, "nb_var": nb_var, "nb_frac_zero": nb_zero,
        "nb_est_dispersion": est_disp,
    })
    summary.to_csv("results/summary.csv", index=False)

    print(f"Genes: {N_GENES}, samples: {N_SAMPLES}")
    print(f"Median Poisson var/mean ratio: {np.median(p_var / p_mean):.2f} (Poisson ~ 1.0)")
    print(f"Median NB var/mean ratio:      {np.median(nb_var / nb_mean):.2f} (overdispersed > 1)")
    print("Wrote figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()
