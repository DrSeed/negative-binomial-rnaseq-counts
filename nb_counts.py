#!/usr/bin/env python3
# Reusable helpers for simulating and characterising negative-binomial counts.
# Import these in your own analysis, or run this file for a tiny sanity check.
import numpy as np


def simulate_nb(true_mean, dispersion, n_samples, rng=None):
    # Draw an (n_genes x n_samples) NB count matrix with per-gene mean and
    # dispersion. Variance = mean + dispersion * mean^2.
    rng = rng or np.random.default_rng()
    true_mean = np.asarray(true_mean, dtype=float)
    dispersion = np.broadcast_to(dispersion, true_mean.shape).astype(float)
    n = 1.0 / dispersion
    p = n / (n + true_mean)
    return rng.negative_binomial(n[:, None], p[:, None], size=(true_mean.size, n_samples))


def estimate_dispersion(mat):
    # Method-of-moments dispersion per gene from a count matrix.
    mean = mat.mean(axis=1)
    var = mat.var(axis=1, ddof=1)
    return np.clip((var - mean) / (mean ** 2 + 1e-9), 1e-4, None)


def overdispersion_ratio(mat):
    # Variance / mean per gene. ~1 means Poisson-like; >1 means overdispersed.
    mean = mat.mean(axis=1)
    var = mat.var(axis=1, ddof=1)
    return var / (mean + 1e-9)


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    mat = simulate_nb([5, 50, 500], dispersion=0.1, n_samples=100, rng=rng)
    print("shape:", mat.shape)
    print("dispersion est:", np.round(estimate_dispersion(mat), 3))
    print("var/mean:", np.round(overdispersion_ratio(mat), 2))
