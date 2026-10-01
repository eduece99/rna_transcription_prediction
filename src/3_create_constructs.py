


import pandas as pd
import numpy as np
import xgboost as xgb
from pathlib import Path

from utils import extract_rna_features, motif_hits

# PROMPT  what about some sort of generative neural network which uses my xgboost model as a scoring function to suggest new UTR sequences?
# PROMPT  please provide me with code instead for a genetic algorithm which uses my xgboost model as a scoring function to suggest new UTR sequences


def valid_utr5(seq: str) -> bool:
    """Check whether a 5' UTR is valid length and alphabet.

    Args:
        seq: Candidate 5' UTR sequence.

    Returns:
        True if the sequence length is 50 nt and uses only A, C, G, U.
    """
    return len(seq) == 50 and set(seq) <= set("ACGU")

def valid_utr3(seq: str) -> bool:
    """Check whether a 3' UTR is valid length and alphabet.

    Args:
        seq: Candidate 3' UTR sequence.

    Returns:
        True if the sequence length is 120 nt and uses only A, C, G, U.
    """
    return len(seq) == 120 and set(seq) <= set("ACGU")

def score_candidate(utr5: str, utr3: str, model: xgb.XGBRegressor) -> float:
    """Predict day-28 liver RNA for a single candidate UTR pair.

    Args:
        utr5: 5' UTR sequence to evaluate.
        utr3: 3' UTR sequence to evaluate.
        model: Trained XGBoost  model.

    Returns:
        Model-predicted day-28 liver RNA value for the candidate pair.
    """
    X = pd.DataFrame([extract_rna_features(utr5, utr3)])
    return float(model.predict(X)[0])

def biological_penalty(utr5: str, utr3: str) -> float:
    """Penalize biologically unfavorable sequence patterns.

    This keeps the search near credible UTR structures by downweighting motifs
    associated with reduced translation or transcript durability.

    Args:
        utr5: Candidate 5' UTR sequence.
        utr3: Candidate 3' UTR sequence.

    Returns:
        A penalty score that should be subtracted from model output.
    """
    penalty = 0.0

    # 5'UTR: upstream AUGs are bad for translation
    penalty += 500.0 * sum(1 for i in range(len(utr5) - 2) if utr5[i:i+3] == "AUG")

    # 3'UTR: AU-rich decay elements and miR-122 seed reduce durability
    penalty += 250.0 * motif_hits(utr3, "AUUUA")
    penalty += 350.0 * motif_hits(utr3, "ACACUCC")

    # Simple extra penalty for very AU-rich / weakly structured tails
    if (utr3.count("A") + utr3.count("U")) / len(utr3) > 0.72:
        penalty += 80.0

    return penalty

def mutate_seq(seq: str, rng: np.random.Generator, p: float = 0.03) -> str:
    """Apply random point mutations to a sequence with a fixed probability.

    Args:
        seq: RNA sequence to mutate.
        rng: NumPy random generator used to sample mutation positions.
        p: Probability of mutating each base.

    Returns:
        Mutated sequence of the same length as the input sequence.
    """
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if rng.random() < p:
            choices = [b for b in "ACGU" if b != seq_list[i]]
            seq_list[i] = rng.choice(choices)
    return "".join(seq_list)

def generate_new_utr_pairs_from_model_med(
    model_med: xgb.XGBRegressor,
    model_low: xgb.XGBRegressor,
    model_high: xgb.XGBRegressor,
    utr_library_df: pd.DataFrame,
    n_designs: int = 10,
    n_rounds_mutation : int = 4,
    seed: int = 42,
) -> pd.DataFrame:
    """Generate candidate UTR pairs by mutating library seeds under the XGBoost model.

    This is a model-guided search: it starts from known good sequences, iteratively
    mutates them, and keeps the changes that improve the median prediction while
    respecting biological constraints.

    Args:
        model_med: Trained XGBoost model for day-28 liver RNA prediction.
        model_low: Trained XGBoost model to compute lower quantile predictions.
        model_high: Trained XGBoost model to compute upper quantile predictions.
        utr_library_df: DataFrame containing library sequences and metadata.
        n_designs: Number of candidate designs to return.
        seed: Random seed for reproducible mutation behavior.

    Returns:
        DataFrame with candidate UTR pairs and predicted scores.
    """
    rng = np.random.default_rng(seed)

    # Seed candidates from the library with strong baseline biology.
    # We deliberately avoid obvious bad patterns: upstream AUGs, AUUUA-heavy tails, miR-122 motif.

    print("Filtering library for valid seed constructs...")
    good_seeds = []
    for _, row in utr_library_df.iterrows():
        u5 = row["utr5"]
        u3 = row["utr3"]
        if not valid_utr5(u5) or not valid_utr3(u3):
            continue
        if sum(1 for i in range(len(u5) - 2) if u5[i:i+3] == "AUG") > 0:
            continue
        if motif_hits(u3, "AUUUA") > 1:
            continue
        if motif_hits(u3, "ACACUCC") > 0:
            continue
        good_seeds.append((u5, u3))

    if not good_seeds:
        raise ValueError("No valid seed constructs found in utr_library.csv")

    print(f"Found {len(good_seeds)} valid seed constructs for mutation search.")
    print(f"Starting mutation search with {n_rounds_mutation} rounds per seed...")
    candidates = []
    for u5_start, u3_start in good_seeds:
        current_u5 = u5_start
        current_u3 = u3_start

        current_score = score_candidate(current_u5, current_u3, model_med) - biological_penalty(current_u5, current_u3)

        print(f"Starting seed: {current_u5} | {current_u3} | score: {current_score:.2f}")

        for _ in range(n_rounds_mutation):
            trial_u5 = mutate_seq(current_u5, rng, p=0.04)
            trial_u3 = mutate_seq(current_u3, rng, p=0.025)

            if not valid_utr5(trial_u5) or not valid_utr3(trial_u3):
                continue

            trial_score = score_candidate(trial_u5, trial_u3, model_med) - biological_penalty(trial_u5, trial_u3)

            if trial_score > current_score:
                current_u5 = trial_u5
                current_u3 = trial_u3
                current_score = trial_score

        print(f"Final candidate: {current_u5} | {current_u3} | score: {current_score:.2f}")
        candidates.append((current_u5, current_u3, current_score))

    # Keep the best unique candidate pairs
    ranked = sorted(candidates, key=lambda x: x[2], reverse=True)
    out = []
    seen = set()
    for u5, u3, score in ranked:
        key = (u5, u3)
        if key in seen:
            continue
        seen.add(key)
        # design_id,utr5,utr3,predicted_day28,lower80,upper80,rationale
        out.append({
            "design_id": f"design_{len(out)+1}",
            "utr5": u5,
            "utr3": u3,
            "predicted_day28": score_candidate(u5, u3, model_med),
            "lower80": score_candidate(u5, u3, model_low),
            "upper80": score_candidate(u5, u3, model_high),
            "rationale": "Model-guided mutation from library seed"
        })
        if len(out) >= n_designs:
            break

    return pd.DataFrame(out)


if __name__ == "__main__":

    utr_df = pd.read_csv("utr_library.csv")

    model_med = xgb.XGBRegressor()
    model_med.load_model("models/model_med.json")
    model_low = xgb.XGBRegressor()
    model_low.load_model("models/model_low.json")
    model_high = xgb.XGBRegressor()
    model_high.load_model("models/model_high.json")

    designs_df =generate_new_utr_pairs_from_model_med(
        model_med, model_low, model_high, 
        utr_df, n_designs=10, seed=42,
        n_rounds_mutation=20
    )

    # output for existing constructs as a baseline reference for comparison with new designs
    baseline_list = []
    for u5, u3, construct_id in utr_df[["utr5", "utr3", "construct_id"]].values:
            key = (u5, u3)
            # design_id,utr5,utr3,predicted_day28,lower80,upper80,rationale
            baseline_list.append({
                "design_id": construct_id,
                "utr5": u5,
                "utr3": u3,
                "predicted_day28": score_candidate(u5, u3, model_med),
                "lower80": score_candidate(u5, u3, model_low),
                "upper80": score_candidate(u5, u3, model_high),
                "rationale": "Existing library construct with predictions attached as a baseline reference"
            })
    baseline_df = pd.DataFrame(baseline_list)

    all_designs_df = pd.concat([designs_df, baseline_df], ignore_index=True)

    all_designs_df.to_csv("designs.csv", index=False)
