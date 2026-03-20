import glob
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from hmmlearn.hmm import MultinomialHMM
from typing import List, Dict, Optional, Tuple, Any

from hmm_models import run_HMM


def fit_vaf_mixture(
    x: np.ndarray,
    random_state: int = 0,
    bic_delta: float = 0.0,
) -> Tuple[GaussianMixture, np.ndarray, np.ndarray, bool, Dict[str, Any]]:
    """
    Fit VAF distribution with both 1- and 2-component Gaussian mixtures and
    select the better model using BIC.

    Parameters
    ----------
    x : np.ndarray
        1D array of VAF values.
    random_state : int
        Random seed for reproducibility.
    bic_delta : float
        Require BIC(1) - BIC(2) > bic_delta to accept k=2.
        Use 0 for "pick the smallest BIC". Use e.g. 10 for stronger evidence.

    Returns
    -------
    gmm : GaussianMixture
        The chosen model (k=2 if preferred, else k=1).
    post_p : np.ndarray
        Posterior probabilities. Shape:
          - (n, 2) if k=2 chosen
          - (n, 1) if k=1 chosen (all ones)
    mu : np.ndarray
        Component means. Shape:
          - (2,) if k=2 chosen
          - (1,) if k=1 chosen
    use_k2 : bool
        True if k=2 selected by BIC, else False (treat as 1-clone downstream).
    info : dict
        Diagnostics: BIC/AIC values, weights, covariances, etc.
    """
    x = np.asarray(x, dtype=float).reshape(-1, 1)

    g1 = GaussianMixture(n_components=1, covariance_type="full", random_state=random_state)
    g2 = GaussianMixture(n_components=2, covariance_type="full", random_state=random_state)

    g1.fit(x)
    g2.fit(x)

    bic1, bic2 = g1.bic(x), g2.bic(x)
    aic1, aic2 = g1.aic(x), g2.aic(x)

    # Smaller BIC is better. Accept k=2 only if it improves BIC by > bic_delta.
    use_k2 = (bic1 - bic2) > bic_delta

    if use_k2:
        gmm = g2
        post_p = gmm.predict_proba(x)          # (n, 2)
        mu = gmm.means_.ravel()               # (2,)
    else:
        gmm = g1
        post_p = np.ones((x.shape[0], 1))     # (n, 1) "all in one component"
        mu = gmm.means_.ravel()               # (1,)

    info = {
        "bic1": float(bic1),
        "bic2": float(bic2),
        "aic1": float(aic1),
        "aic2": float(aic2),
        "bic_delta_used": float(bic_delta),
        "use_k2": bool(use_k2),
        "weights": gmm.weights_.ravel().tolist(),
        "means": mu.tolist(),
        "covariances": np.asarray(gmm.covariances_).ravel().tolist(),
    }
    return gmm, post_p, mu, use_k2, info

def encode_obs(code_as: pd.Series):
    """
    Map emissions to integers for MultinomialHMM:
      "T>N" -> 0
      "A>N" -> 1
    Returns obs as shape (n,1) integers and mapping dict.
    """
    mapping = {"T>N": 0, "A>N": 1}
    obs = code_as.map(mapping)
    if obs.isna().any():
        bad = code_as[obs.isna()].unique()
        raise ValueError(f"Unexpected code_as values: {bad}")
    return obs.to_numpy(dtype=int).reshape(-1, 1), mapping

def map_clonal_to_subclonal_next(
    subclon_coords: np.ndarray,
    clon_coords: np.ndarray
) -> np.ndarray:
    """
    Map each clonal coordinate to the "next" subclonal coordinate,
    reproducing the R behavior:

    This function is used to pair clonal mutations with the nearest
    subsequent subclonal mutation along the genome.

    Parameters
    ----------
    subclon_coords : np.ndarray
        Must be sorted in ascending order.

    clon_coords : np.ndarray
        Does not need to be sorted, but should correspond to the
        same coordinate system as subclon_coords.

    Returns
    -------
    idx : np.ndarray
        Integer array of indices (same length as clon_coords).
        idx[i] is the index in `subclon_coords` that the i-th clonal
        mutation is mapped to.

    """
    # Ensure numeric arrays
    subclon_coords = np.asarray(subclon_coords, dtype=float)
    clon_coords = np.asarray(clon_coords, dtype=float)

    # Find index of rightmost subclonal <= clonal (R's findInterval)
    j = np.searchsorted(subclon_coords, clon_coords, side="right") - 1

    # Move to the next subclonal mutation (+1)
    idx = j + 1

    # Clip to valid range [0, len(subclon_coords)-1]
    idx = np.clip(idx, 0, len(subclon_coords) - 1)

    return idx


def safe_div(a, b):
    return a / b if b != 0 else 0.0


# ----------------------------
# Per-sample processing
# ----------------------------
def process_one_file(file_path: Path, out_dir: Path, random_state=0):
    """
    Returns a dict of summary metrics for the final master table.
    Also writes intermediate files.
    """
    sample = file_path.stem
    print(f"Processing {sample}")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(file_path, sep=",", low_memory=False)
    df["chr"] = df["chr"].astype(str)

    #select bi-allelic autosomal A>N/T>N
    mask = (df["code_multi"] == "B") & (df["code_as"].isin(["A>N", "T>N"])) & (df["chr"] != "X")
    SM = df.loc[mask].copy()

    # Ensure numeric
    SM["chr_num"] = pd.to_numeric(SM["chr"], errors="coerce")
    SM = SM.dropna(subset=["chr_num"])
    SM["chr_num"] = SM["chr_num"].astype(int)
    SM["altCount"] = pd.to_numeric(SM["altCount"], errors="coerce")
    SM["refCount"] = pd.to_numeric(SM["refCount"], errors="coerce")
    SM = SM.dropna(subset=["altCount", "refCount"])

    # calculate VAF and assign coord along through chromosomes
    SM["vaf"] = SM["altCount"] / (SM["altCount"] + SM["refCount"])
    SM["coord"] = SM["chr_num"] * 1e9 + pd.to_numeric(SM["pos"], errors="coerce")
    SM = SM.dropna(subset=["coord"])
    SM = SM.sort_values("coord").reset_index(drop=True)

    # fit mixture on VAF
    gmm, post_p, mu, use_k2, mix_info = fit_vaf_mixture(SM["vaf"].to_numpy(), bic_delta=10)
    print(use_k2)
    print(mix_info)

    Ncl = 2
    if not use_k2:
        Ncl = 1
        final_row = write_minimal_outputs(sample, SM, out_dir)
        subsets_fit = np.nan

    else:
        # Determine which component is subclonal vs clonal by mean values
        subclon_idx, clon_idx = (0, 1) if mu[0] < mu[1] else (1, 0)

        n_subclonal = int((post_p[:, subclon_idx] > 0.8).sum())
        n_clonal = int((post_p[:, clon_idx] > 0.8).sum())

        # Boundaries MAX/MIN of vaf for clonal and subclonal (only among confidently assigned)
        vaf_subclon_conf = SM["vaf"].to_numpy()[post_p[:, subclon_idx] > 0.8]
        vaf_clon_conf = SM["vaf"].to_numpy()[post_p[:, clon_idx] > 0.8]
        if len(vaf_subclon_conf) == 0 or len(vaf_clon_conf) == 0:
            # fallback: mixture didn't separate well
            MAX = np.nan
            MIN = np.nan
        else:
            MAX = float(min(vaf_subclon_conf.max(), vaf_clon_conf.max()))
            MIN = float(max(vaf_subclon_conf.min(), vaf_clon_conf.min()))

        # medians at posterior > 0.5 
        median_subclon = float(np.median(SM["vaf"].to_numpy()[post_p[:, subclon_idx] > 0.5])) if (post_p[:, subclon_idx] > 0.5).any() else np.nan
        median_clon = float(np.median(SM["vaf"].to_numpy()[post_p[:, clon_idx] > 0.5])) if (post_p[:, clon_idx] > 0.5).any() else np.nan

        # subset clonal/subclonal observations for HMM using posterior > 0.8
        obs_subclon = SM.loc[post_p[:, subclon_idx] > 0.8].copy()
        obs_clon = SM.loc[post_p[:, clon_idx] > 0.8].copy()

        if (len(obs_subclon) < 800) or (len(obs_clon) < 800):
            Ncl = 1
            final_row = write_minimal_outputs(sample, SM, out_dir)
            subsets_fit = np.nan

        else:            
            state_names = ["A1", "A2", "A3"]

            # Train/Decode HMMs (only if enough data)
            emissions_subclon_T = [np.nan, np.nan, np.nan]
            emissions_clon_T = [np.nan, np.nan, np.nan]
            counts_subclon = {"A1": 0, "A2": 0, "A3": 0}
            counts_clon = {"A1": 0, "A2": 0, "A3": 0}
            subsets_fit = "Good"

            _, emissions_clon, _, Z_as_clon, cats_clon = run_HMM(obs_clon[["code_as"]].to_numpy(), "HMM_ploidy")
            _, emissions_subclon, _, Z_as_subclon, cats_subclon = run_HMM(obs_subclon[["code_as"]].to_numpy(), "HMM_ploidy")

            t_idx_clon = int(np.where(cats_clon == "T>N")[0][0])
            t_idx_subclon = int(np.where(cats_subclon == "T>N")[0][0])

            obs_subclon["state_subclonal"] = Z_as_subclon
            obs_clon["state_clonal"] = Z_as_clon

            emissions_subclon_T = emissions_subclon[:, t_idx_subclon].tolist()
            emissions_clon_T = emissions_clon[:, t_idx_clon].tolist()

            counts_subclon = pd.Series(Z_as_subclon).value_counts().to_dict()
            counts_clon = pd.Series(Z_as_clon).value_counts().to_dict()

            # mark subset fit "bad" if emissions deviate too strong from the expectation
            print(emissions_subclon)
            subclon_A1, subclon_A2, subclon_A3 = emissions_subclon[:, 0]
            clon_A1, clon_A2, clon_A3 = emissions_clon[:, 0]

            if (
                subclon_A1 < 0.67 or clon_A1 < 0.67 or
                subclon_A3 > 0.33 or clon_A3 > 0.33 or
                not (0.4 <= subclon_A2 <= 0.6) or
                not (0.4 <= clon_A2 <= 0.6)
            ):
                subsets_fit = "Bad"

            # Also mark mixture fit "bad" if any state has <100 sites in either subset
            for s in state_names:
                if counts_subclon.get(s, 0) < 100 or counts_clon.get(s, 0) < 100:
                    subsets_fit = "Bad"

            # ----------------------------
            # PAir subclonal states with clonal states
            # ----------------------------
            cloneas_df = pd.DataFrame(columns=["subclonal", "clonal", "count", "prop"])
            o_e = {"11": 0.0, "33": 0.0, "31": 0.0, "13": 0.0}

            if ("state_subclonal" in obs_subclon.columns) and ("state_clonal" in obs_clon.columns) and (len(obs_subclon) > 0) and (len(obs_clon) > 0):
                obs_subclon = obs_subclon.sort_values("coord").reset_index(drop=True)
                obs_clon = obs_clon.sort_values("coord").reset_index(drop=True)

                idx = map_clonal_to_subclonal_next(obs_subclon["coord"].to_numpy(), obs_clon["coord"].to_numpy())

                AB = pd.DataFrame({
                    "subclonal": obs_subclon.loc[idx, "state_subclonal"].to_numpy(),
                    "clonal": obs_clon["state_clonal"].to_numpy(),
                })

                tab = AB.value_counts().reset_index(name="count")
                tab["prop"] = tab["count"] / tab["count"].sum()
                cloneas_df = tab.rename(columns={"level_0": "subclonal", "level_1": "clonal"})

                # Observed/expected ratios 
                # S1 = subclonal props by A1/A2/A3, S2 = clonal props
                all1 = sum(counts_subclon.get(s, 0) for s in state_names)
                all2 = sum(counts_clon.get(s, 0) for s in state_names)
                S1 = {s: safe_div(counts_subclon.get(s, 0), all1) for s in state_names}
                S2 = {s: safe_div(counts_clon.get(s, 0), all2) for s in state_names}

                # observed combined props
                obs_comb = {(r.subclonal + r.clonal): r.prop for r in tab.itertuples(index=False)}

                def oe(subclon_state, clon_state):
                    key = subclon_state + clon_state  # e.g. "A1A3"
                    exp = S1[subclon_state] * S2[clon_state]
                    return safe_div(obs_comb.get(key, 0.0), exp)

                # Map to your "11/33/31/13" notion (A1=1, A3=3)
                o_e["11"] = oe("A1", "A1")
                o_e["33"] = oe("A3", "A3")
                o_e["31"] = oe("A3", "A1")
                o_e["13"] = oe("A1", "A3")

                Mixture_tumours = "mixture_tumours"
                if o_e["11"] > 1.5 and o_e["33"] > 1.5 and o_e["31"] < 0.6 and o_e["13"] < 0.6:
                    Mixture_tumours = "one_tumour*"
                if o_e["11"] > 1.8 and o_e["33"] > 1.8 and o_e["31"] < 0.6 and o_e["13"] < 0.6:
                    Mixture_tumours = "one_tumour"
                if o_e["31"] > 2 and o_e["13"] > 2:
                    Mixture_tumours = "Symmetric"

                # ----------------------------
                # Write intermediate outputs 
                # ----------------------------

                # clonesize
                clonesize = pd.DataFrame([{
                    "sample": sample,
                    "median_subclon": median_subclon,
                    "median_clon": median_clon,
                    "vaf_max_overlap": MAX,
                    "vaf_min_overlap": MIN,
                    "n_subclon": n_subclonal,
                    "n_clon": n_clonal,
                    }])
                clonesize.to_csv(out_dir / f"{sample}.cloneSize.tsv", sep="\t", index=False)

                # BW
                bw = pd.DataFrame([{
                    "sample": sample,
                    "emission_TtoN_subclon_A1": emissions_subclon_T[0],
                    "emission_TtoN_subclon_A2": emissions_subclon_T[1],
                    "emission_TtoN_subclon_A3": emissions_subclon_T[2],
                    "emission_TtoN_clon_A1": emissions_clon_T[0],
                    "emission_TtoN_clon_A2": emissions_clon_T[1],
                    "emission_TtoN_clon_A3": emissions_clon_T[2],
                    "count_subclon_A1": counts_subclon.get("A1", 0),
                    "count_subclon_A2": counts_subclon.get("A2", 0),
                    "count_subclon_A3": counts_subclon.get("A3", 0),
                    "count_clon_A1": counts_clon.get("A1", 0),
                    "count_clon_A2": counts_clon.get("A2", 0),
                    "count_clon_A3": counts_clon.get("A3", 0),
                    }])
                bw.to_csv(out_dir / f"{sample}.bw.tsv", sep="\t", index=False)

                # CloneAS
                cloneas_df.to_csv(out_dir / f"{sample}.cloneAs.tsv", sep="\t", index=False)

                # ----------------------------
                # Return final per-sample summary row
                # ----------------------------
                final_row = {
                    "Sample": sample,
                    "Ncl": Ncl,
                    "subsets_fit": subsets_fit,
                    "o_e_11": o_e["11"],
                    "o_e_33": o_e["33"],
                    "o_e_31": o_e["31"],
                    "o_e_13": o_e["13"],
                    "Mixture_tumours": Mixture_tumours,
                    }
    return final_row


def write_minimal_outputs(sample: str, SM: pd.DataFrame, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([{
        "sample": sample,
        "median_subclonclon": np.nan,
        "median_clon": np.nan,
        "vaf_max_overlap": np.nan,
        "vaf_min_overlap": np.nan,
        "n_subclon": 0,
        "n_clonal": 0,
    }]).to_csv(out_dir / f"{sample}.clonesize.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "sample": sample,
        "emission_TtoN_subclon_A1": np.nan, "emission_TtoN_subclon_A2": np.nan, "emission_TtoN_subclon_A3": np.nan,
        "emission_TtoN_clon_A1": np.nan, "emission_TtoN_clon_A2": np.nan, "emission_TtoN_clon_A3": np.nan,
        "count_subclon_A1": 0, "count_subclon_A2": 0, "count_subclon_A3": 0,
        "count_clon_A1": 0, "count_clon_A2": 0, "count_clon_A3": 0,
    }]).to_csv(out_dir / f"{sample}.bw.tsv", sep="\t", index=False)

    pd.DataFrame(columns=["subclonal", "clonal", "count", "prop"]).to_csv(out_dir / f"{sample}.cloneas.tsv", sep="\t", index=False)

    return {
        "Sample": sample,
        "Ncl": 1,
        "subsets_fit": "1_clone",
        "o_e_11": "1_clone",
        "o_e_33": "1_clone",
        "o_e_31": "1_clone",
        "o_e_13": "1_clone",
        "Mixture_tumours": "1_clone",
    }


# ----------------------------
# Batch runner
# ----------------------------
def main():
    in_files = glob.glob("../../data/encoded_data_all/*.hmm")
    out_dir = Path("../results/HMM_ploidy")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for fn in in_files:
        row = process_one_file(Path(fn), out_dir, random_state=0)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "../Ploidy_summary.tsv", sep="\t", index=False,na_rep="NA")
    print(summary.head())


if __name__ == "__main__":
    main()