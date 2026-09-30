"""
Common functions used by multiple scripts in the ScionXBio UTR design pipeline.
"""

import re
import pandas as pd


def motif_hits(seq: str, motif: str) -> int:
    """Count non-overlapping occurrences of a motif in a sequence.

    Args:
        seq: RNA sequence to search.
        motif: Motif string, such as "AUUUA" or "ACACUCC".

    Returns:
        Number of matches found using regex-based searching.
    """
    return len(re.findall(re.escape(motif), seq))


def extract_rna_features(utr5: str, utr3: str) -> dict:
    """Build a compact feature table describing 5' and 3' UTR biology.

    The feature set captures basic composition, cap-proximal GC content,
    upstream AUG counts, Kozak-like context, AU-rich decay motifs, and the
    miR-122 target seed that is relevant in liver.

    Args:
        utr5: 50-nt 5' untranslated region sequence.
        utr3: 120-nt 3' untranslated region sequence.

    Returns:
        Dictionary mapping feature names to numeric values.
    """
    feats = {}

    # --- general helpers ---
    def gc_content(seq: str) -> float:
        if len(seq) == 0:
            return 0.0
        return (seq.count("G") + seq.count("C")) / len(seq)

    def au_content(seq: str) -> float:
        if len(seq) == 0:
            return 0.0
        return (seq.count("A") + seq.count("U")) / len(seq)

    # --- 5'UTR features ---
    feats["utr5_gc"] = gc_content(utr5)
    feats["utr5_au"] = au_content(utr5)
    feats["utr5_a_count"] = utr5.count("A")
    feats["utr5_c_count"] = utr5.count("C")
    feats["utr5_g_count"] = utr5.count("G")
    feats["utr5_u_count"] = utr5.count("U")

    # cap-proximal GC / structure
    feats["utr5_capprox_gc"] = gc_content(utr5[:12])
    feats["utr5_capprox_au"] = au_content(utr5[:12])

    # upstream AUGs (can inhibit translation)
    upstream_aug_count = 0
    for i in range(len(utr5) - 2):
        if utr5[i:i+3] == "AUG":
            upstream_aug_count += 1
    feats["utr5_upstream_aug_count"] = upstream_aug_count

    # Kozak-like context: purine at -3 relative to AUG
    kozak_good = 0
    for i in range(len(utr5) - 2):
        if utr5[i:i+3] == "AUG" and i >= 3 and utr5[i-3] in ["A", "G"]:
            kozak_good += 1
    feats["utr5_kozak_purine_minus3"] = kozak_good

    # --- 3'UTR features ---
    feats["utr3_gc"] = gc_content(utr3)
    feats["utr3_au"] = au_content(utr3)
    feats["utr3_a_count"] = utr3.count("A")
    feats["utr3_c_count"] = utr3.count("C")
    feats["utr3_g_count"] = utr3.count("G")
    feats["utr3_u_count"] = utr3.count("U")

    # AU-rich decay motifs
    feats["utr3_auuua_count"] = motif_hits(utr3, "AUUUA")
    feats["utr3_auuu_count"] = motif_hits(utr3, "AUUU")  # not explicitly in primer.md
    feats["utr3_au_rich_blocks"] = len(re.findall(r"AUUUA|AUUUA", utr3))  # not explicitly in primer.md

    # miR-122 seed motif in liver
    feats["utr3_mir122_seed_acacucc"] = motif_hits(utr3, "ACACUCC")

    # simple compositional decay proxies
    feats["utr3_poly_a_run_max"] = max((len(m) for m in re.findall(r"A+", utr3)), default=0)
    feats["utr3_poly_u_run_max"] = max((len(m) for m in re.findall(r"U+", utr3)), default=0)

    return feats



# Example: apply to a dataframe
def generate_features_from_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the RNA feature extractor row-wise across a construct table.

    Args:
        df: DataFrame containing at least "utr5" and "utr3" columns.

    Returns:
        DataFrame with one row per construct and one column per engineered RNA feature.
    """
    feature_df = df.apply(
        lambda row: pd.Series(extract_rna_features(row["utr5"], row["utr3"])),
        axis=1
    )
    #return pd.concat([df, feature_df], axis=1)
    return feature_df