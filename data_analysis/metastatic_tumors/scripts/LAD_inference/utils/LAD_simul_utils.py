from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Iterable, Mapping


def get_autosome_length_weights(chr_lengths: Mapping[str, int] | None = None) -> np.ndarray:
    """
    Return normalized chromosome-length weights for chr1..chr22.

    Parameters
    ----------
    chr_lengths : dict or None
        Mapping like {"1": length1, ..., "22": length22}
        or {"chr1": length1, ..., "chr22": length22}.
        If None, a built-in autosome length dictionary is used.

    Returns
    -------
    weights : np.ndarray
        Shape (22,), sums to 1.
    """
    if chr_lengths is None:
        # Replace with your preferred assembly-specific lengths if needed.
        # These are common GRCh38 autosome lengths in bp.
        chr_lengths = {
            "1": 248956422,
            "2": 242193529,
            "3": 198295559,
            "4": 190214555,
            "5": 181538259,
            "6": 170805979,
            "7": 159345973,
            "8": 145138636,
            "9": 138394717,
            "10": 133797422,
            "11": 135086622,
            "12": 133275309,
            "13": 114364328,
            "14": 107043718,
            "15": 101991189,
            "16": 90338345,
            "17": 83257441,
            "18": 80373285,
            "19": 58617616,
            "20": 64444167,
            "21": 46709983,
            "22": 50818468,
        }

    vals = []
    for i in range(1, 23):
        k1 = str(i)
        k2 = f"chr{i}"
        if k1 in chr_lengths:
            vals.append(chr_lengths[k1])
        elif k2 in chr_lengths:
            vals.append(chr_lengths[k2])
        else:
            raise ValueError(f"Missing chromosome length for chromosome {i}")

    weights = np.asarray(vals, dtype=float)
    if np.any(weights <= 0):
        raise ValueError("All chromosome lengths must be positive.")
    weights /= weights.sum()
    return weights

def simulate_divisions(input_vector: Iterable[int], N: int, rng: np.random.Generator | None = None) -> np.ndarray:
    """
    Simulation of segregation of lesions in subsequent divisions.

    Starts from a vector of chromosome homolog counts (typically 22 entries, all 2).
    Meaning that both chromosomes in each pair contain lesions in the beginning.
    Applies N divisions with random segregation and update the number of homolog chromosomes in each pair
    using these rules:
      - if current value == 2 (both chroms with multi sites): 
                sample new value from {2,1,0} with probs {1/4,1/2,1/4}
      - if current value == 1 (one chrom with multi sites): 
                sample new value from {1,0} with probs {1/2,1/2}
      - if current value == 0: stays 0

    After all iterations, if all entries are zero, force the first entry to 1.
    Parameters
    ----------
    input_vector : iterable of int
        Initial state vector, e.g. [2] * 22
    N : int
        Number of iterations / divisions
    rng : np.random.Generator or None
        Random number generator

    Returns
    -------
    np.ndarray
        Updated vector.
    """
    if rng is None:
        rng = np.random.default_rng()

    vec = np.asarray(list(input_vector), dtype=int).copy()

    if N == 0:
        return vec

    for _ in range(N):
        for i in range(len(vec)):
            current_value = vec[i]

            if current_value == 2:
                vec[i] = rng.choice([2, 1, 0], p=[1/4, 1/2, 1/4])
            elif current_value == 1:
                vec[i] = rng.choice([1, 0], p=[1/2, 1/2])
            else:
                vec[i] = 0

    if vec.sum() == 0:
        vec[0] = 1

    return vec


def simulate_one_C(
    M: int,
    D: int,
    base_weights: np.ndarray,
    iv: np.ndarray | None = None,
    rng: np.random.Generator | None = None,
) -> int:
    """
    Simulate one realization of C given M and D.

    Parameters
    ----------
    M : int
        Number of multiallelic sites
    D : int
        LAD state
    base_weights : np.ndarray
        Chromosome baseline probabilities (e.g. from chromosome lengths), shape (22,)
    iv : np.ndarray or None
        Initial chromosome state vector, default = [2] * 22
    rng : np.random.Generator or None
        Random number generator

    Returns
    -------
    int
        Simulated number of chromosomes carrying multiallelic sites (C)
    """
    if rng is None:
        rng = np.random.default_rng()

    if iv is None:
        iv = np.repeat(2, 22)

    if len(base_weights) != 22:
        raise ValueError("base_weights must have length 22.")
    if M < 0:
        raise ValueError("M must be >= 0.")
    if D < 1:
        raise ValueError("D must be >= 1.")

    if M == 0:
        return 0

    div = D - 1 #because we start from cell that is LAD=1
    result_vector = simulate_divisions(iv, div, rng=rng)

    # Same structure as your R code:
    # probabilities ∝ chromosome_length_weight * homolog_count
    active_weights = base_weights * result_vector.astype(float)

    if active_weights.sum() == 0:
        return 0

    probabilities = active_weights / active_weights.sum()

    selected_chr = rng.choice(
        np.arange(22),   # 0-based internally
        size=M,
        replace=True,
        p=probabilities
    )

    C_sim = np.unique(selected_chr).size
    return int(C_sim)


def simulate_pc_given_m_d(
    M: int,
    D: int,
    n_iter: int = 10000,
    base_weights: np.ndarray | None = None,
    iv: np.ndarray | None = None,
    seed: int | None = None,
    pseudocount: float = 0.0,
    complete_grid: bool = True,
) -> pd.DataFrame:
    """
    Estimate P(C | M, D) by Monte Carlo.

    Parameters
    ----------
    M : int
        Number of multiallelic sites
    D : int
        LAD state
    n_iter : int
        Number of simulations
    base_weights : np.ndarray or None
        Baseline chromosome probabilities. If None, uses chromosome-length weights.
    iv : np.ndarray or None
        Initial homolog-count vector, default [2] * 22
    seed : int or None
        Random seed
    pseudocount : float
        Optional smoothing added to counts before normalization
    complete_grid : bool
        If True, return rows for all C in 0..min(M,22), filling missing probs with 0
        (or with pseudocount-adjusted small values if pseudocount > 0)

    Returns
    -------
    pd.DataFrame
        Columns: M, C, D, prob
    """
    rng = np.random.default_rng(seed)

    if base_weights is None:
        base_weights = get_autosome_length_weights()

    if iv is None:
        iv = np.repeat(2, 22)

    sim_C = np.array(
        [simulate_one_C(M=M, D=D, base_weights=base_weights, iv=iv, rng=rng) for _ in range(n_iter)],
        dtype=int
    )

    counts = pd.Series(sim_C).value_counts().sort_index()

    if complete_grid:
        c_min = 0 if M == 0 else 1
        c_max = min(M, 22)
        all_c = pd.Index(range(c_min, c_max + 1))
        counts = counts.reindex(all_c, fill_value=0)

    counts = counts.astype(float)
    if pseudocount > 0:
        counts = counts + pseudocount

    probs = counts / counts.sum()

    out = pd.DataFrame({
        "M": M,
        "C": probs.index.astype(int),
        "LAD": D,
        "prob": probs.values.astype(float),
    })

    return out


def build_likelihood_df(
    M_values: Iterable[int],
    D_values: Iterable[int] = (1, 2, 3, 4, 5),
    n_iter: int = 10000,
    chr_lengths: Mapping[str, int] | None = None,
    seed: int | None = 1,
    pseudocount: float = 0.0,
    complete_grid: bool = True,
) -> pd.DataFrame:
    """
    Build the full lookup table likelihood_df with probabilities P(C | M, D).

    Parameters
    ----------
    M_values : iterable of int
        Values of the number of multiallelic sites(M) to simulate
    D_values : iterable of int
        LAD states
    n_iter : int
        Number of simulations per (M, D)
    chr_lengths : dict or None
        Chromosome length mapping; if None, uses built-in defaults
    seed : int or None
        Random seed for reproducibility
    pseudocount : float
        Optional smoothing count
    complete_grid : bool
        Whether to include zero-probability rows for unseen C values

    Returns
    -------
    pd.DataFrame
        Columns: M, C, D, prob
    """
    base_weights = get_autosome_length_weights(chr_lengths)
    rng = np.random.default_rng(seed)

    out = []

    for M in M_values:
        for D in D_values:
            print(f"{M}-{D}")
            # Create a fresh seed per (M, D) for reproducibility but independence
            local_seed = int(rng.integers(0, 2**32 - 1))
            df_md = simulate_pc_given_m_d(
                M=M,
                D=D,
                n_iter=n_iter,
                base_weights=base_weights,
                iv=np.repeat(2, 22),
                seed=local_seed,
                pseudocount=pseudocount,
                complete_grid=complete_grid,
            )
            out.append(df_md)

    likelihood_df = pd.concat(out, ignore_index=True)

    # safety check
    check = likelihood_df.groupby(["M", "LAD"])["prob"].sum().reset_index()
    if not np.allclose(check["prob"].values, 1.0):
        raise RuntimeError("Probabilities do not sum to 1 within each (M, D) block.")

    return likelihood_df

def load_chr_oe(
    dataall_path: str | Path,
    chrom_path: str | Path,
) -> pd.DataFrame:
    """
    Load and prepare CHR_OE equivalent from the two Hartwig tables.

    Mirrors this R logic:
        dataAll <- read.table(...)
        dataAll$OE <- ...
        dataAll <- dataAll[!duplicated(dataAll$sample),]

        chrom <- read.table(...)
        chrom$sample <- substr(chrom$sample, 1, nchar(chrom$sample) - 1)

        cbind(chrom, rowSums(chrom[,c(7:14)]), rowSums(chrom[,c(3:6)])) -> chrom
        names(chrom)[c(15,16)] <- c('MultA', 'BiA')

        merge(dataAll[,c(1,2,9,10)], chrom, by="sample") -> CHR_OE

    Returns
    -------
    pd.DataFrame
        Dataframe analogous to CHR_OE in your R code.
    """
    data_all = pd.read_csv(dataall_path, sep="\t")
    chrom = pd.read_csv(chrom_path, sep="\t")

    # OE
    data_all = data_all.copy()
    data_all["OE"] = data_all["n_multi"] / (
        ((data_all["n_bi"] + data_all["n_multi"]) / 2500000000.0)
        * (data_all["n_bi"] + data_all["n_multi"])
    )

    # remove duplicated samples
    data_all = data_all.drop_duplicates(subset=["sample"]).copy()

    # trim last character from sample in chrom
    chrom = chrom.copy()
    chrom["sample"] = chrom["sample"].str[:-1]

    # R columns 7:14 -> Python iloc[:, 6:14]
    # R columns 3:6  -> Python iloc[:, 2:6]
    chrom["MultA"] = chrom.iloc[:, 6:14].sum(axis=1)
    chrom["BiA"] = chrom.iloc[:, 2:6].sum(axis=1)

    # R dataAll[, c(1,2,9,10)] -> Python positions [0,1,8,9]
    keep_cols = data_all.columns[[0, 1, 8, 9]].tolist()

    chr_oe = data_all[keep_cols].merge(chrom, on="sample", how="inner")
    return chr_oe   
