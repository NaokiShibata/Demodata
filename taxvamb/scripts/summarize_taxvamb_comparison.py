from __future__ import annotations

from pathlib import Path
from collections import defaultdict, deque

import pandas as pd


ROOT = Path("results")
OUT = ROOT / "summary"
OUT.mkdir(parents=True, exist_ok=True)

CHECKM2 = ROOT / "checkm2_all" / "quality_report.tsv"
GTDBTK_DIR = ROOT / "gtdbtk_all"
RNA_SUMMARY = ROOT / "rna_genes" / "rna_summary.tsv"
FASTANI = ROOT / "ani" / "fastani_all_vs_all.tsv"


def strip_ext(name: str) -> str:
    for ext in [".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"]:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def parse_binner(bin_id: str) -> str:
    if "__" not in bin_id:
        return "unknown"
    return bin_id.split("__", 1)[0]


def read_checkm2() -> pd.DataFrame:
    df = pd.read_csv(CHECKM2, sep="\t")

    if "Name" not in df.columns:
        raise ValueError(f"'Name' column not found in {CHECKM2}")

    df["bin_id"] = df["Name"].astype(str).map(strip_ext)
    df["binner"] = df["bin_id"].map(parse_binner)

    return df


def read_rna_summary() -> pd.DataFrame:
    df = pd.read_csv(RNA_SUMMARY, sep="\t")

    required = {
        "bin_id",
        "has_full_rrna_set",
        "has_mimag_trna_set",
        "trna_aa_types",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {RNA_SUMMARY}: {sorted(missing)}")

    return df


def read_gtdbtk() -> pd.DataFrame:
    files = list(GTDBTK_DIR.rglob("*.summary.tsv"))

    if not files:
        return pd.DataFrame(columns=["bin_id", "classification"])

    dfs = []

    for f in files:
        x = pd.read_csv(f, sep="\t")

        if "user_genome" not in x.columns or "classification" not in x.columns:
            continue

        x = x[["user_genome", "classification"]].copy()
        x["bin_id"] = x["user_genome"].astype(str).map(strip_ext)
        dfs.append(x[["bin_id", "classification"]])

    if not dfs:
        return pd.DataFrame(columns=["bin_id", "classification"])

    return pd.concat(dfs, ignore_index=True).drop_duplicates("bin_id")


def split_gtdb_taxonomy(classification: str) -> dict[str, str]:
    ranks = {
        "d__": "domain",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }

    out = {v: "" for v in ranks.values()}

    if not isinstance(classification, str):
        return out

    for item in classification.split(";"):
        item = item.strip()

        for prefix, rank in ranks.items():
            if item.startswith(prefix):
                out[rank] = item[len(prefix):].strip()

    return out


def classify_mimag(row: pd.Series) -> str:
    comp = float(row["Completeness"])
    cont = float(row["Contamination"])

    has_full_rrna = bool(row.get("has_full_rrna_set", False))
    has_trna = bool(row.get("has_mimag_trna_set", False))

    if comp > 90 and cont < 5 and has_full_rrna and has_trna:
        return "HQ"

    if comp >= 50 and cont < 10:
        return "MQ"

    if comp < 50 and cont < 10:
        return "LQ"

    return "contaminated"


def read_fastani_clusters(valid_bins: set[str]) -> pd.DataFrame:
    if not FASTANI.exists():
        return pd.DataFrame({
            "bin_id": sorted(valid_bins),
            "same_mag_group": [
                f"group_{i:05d}" for i, _ in enumerate(sorted(valid_bins), 1)
            ],
        })

    edges = defaultdict(set)

    with FASTANI.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 5:
                continue

            q = strip_ext(Path(parts[0]).name)
            r = strip_ext(Path(parts[1]).name)

            if q == r:
                edges[q].add(r)
                continue

            try:
                ani = float(parts[2])
                mapped = float(parts[3])
                total = float(parts[4])
            except ValueError:
                continue

            matched_fragment_ratio = mapped / total if total else 0.0

            if ani >= 95.0 and matched_fragment_ratio >= 0.30:
                edges[q].add(r)
                edges[r].add(q)

    for b in valid_bins:
        edges[b].add(b)

    seen = set()
    records = []
    group_index = 1

    for b in sorted(valid_bins):
        if b in seen:
            continue

        group_id = f"group_{group_index:05d}"
        group_index += 1

        queue = deque([b])
        seen.add(b)

        while queue:
            x = queue.popleft()
            records.append({"bin_id": x, "same_mag_group": group_id})

            for y in edges[x]:
                if y not in seen:
                    seen.add(y)
                    queue.append(y)

    return pd.DataFrame(records)


def choose_best_bin(df: pd.DataFrame) -> pd.Series:
    rank = {
        "HQ": 0,
        "MQ": 1,
        "LQ": 2,
        "contaminated": 3,
    }

    x = df.copy()
    x["_quality_rank"] = x["quality"].map(rank).fillna(9)

    x = x.sort_values(
        [
            "_quality_rank",
            "Completeness",
            "Contamination",
            "trna_aa_types",
        ],
        ascending=[True, False, True, False],
    )

    return x.iloc[0]


def main() -> None:
    checkm2 = read_checkm2()
    rna = read_rna_summary()
    gtdb = read_gtdbtk()

    merged = checkm2.merge(rna, on="bin_id", how="left")

    for col in ["has_full_rrna_set", "has_mimag_trna_set"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False)

    if "trna_aa_types" in merged.columns:
        merged["trna_aa_types"] = merged["trna_aa_types"].fillna(0)

    merged["quality"] = merged.apply(classify_mimag, axis=1)

    merged = merged.merge(gtdb, on="bin_id", how="left")

    tax = merged["classification"].apply(split_gtdb_taxonomy).apply(pd.Series)
    merged = pd.concat([merged, tax], axis=1)

    groups = read_fastani_clusters(set(merged["bin_id"]))
    merged = merged.merge(groups, on="bin_id", how="left")

    merged.to_csv(OUT / "bin_quality_taxonomy.tsv", sep="\t", index=False)

    table1 = (
        merged.groupby("binner")
        .agg(
            HQ=("quality", lambda s: (s == "HQ").sum()),
            MQ=("quality", lambda s: (s == "MQ").sum()),
            LQ=("quality", lambda s: (s == "LQ").sum()),
            contaminated=("quality", lambda s: (s == "contaminated").sum()),
            total_bins=("bin_id", "count"),
            median_completeness=("Completeness", "median"),
            median_contamination=("Contamination", "median"),
            bins_with_full_rrna_set=("has_full_rrna_set", "sum"),
            bins_with_mimag_trna_set=("has_mimag_trna_set", "sum"),
            median_trna_aa_types=("trna_aa_types", "median"),
        )
        .reset_index()
    )

    table1.to_csv(OUT / "table1_quality_by_binner.tsv", sep="\t", index=False)

    target = merged[merged["quality"].isin(["HQ", "MQ"])].copy()

    rows = []

    for binner, sub in target.groupby("binner"):
        rows.append({
            "binner": binner,
            "phyla": sub["phylum"].replace("", pd.NA).dropna().nunique(),
            "classes": sub["class"].replace("", pd.NA).dropna().nunique(),
            "orders": sub["order"].replace("", pd.NA).dropna().nunique(),
            "families": sub["family"].replace("", pd.NA).dropna().nunique(),
            "genera": sub["genus"].replace("", pd.NA).dropna().nunique(),
            "species": sub["species"].replace("", pd.NA).dropna().nunique(),
        })

    table2 = pd.DataFrame(rows)
    table2.to_csv(
        OUT / "table2_taxonomic_diversity_by_binner.tsv",
        sep="\t",
        index=False,
    )

    group_rows = []

    for group_id, sub in merged.groupby("same_mag_group"):
        best = choose_best_bin(sub)

        binners = sorted(set(sub["binner"]))
        good_binners = sorted(set(sub.loc[sub["quality"].isin(["HQ", "MQ"]), "binner"]))

        group_rows.append({
            "same_mag_group": group_id,
            "binners": ",".join(binners),
            "good_binners": ",".join(good_binners),
            "best_bin": best["bin_id"],
            "best_binner": best["binner"],
            "quality": best["quality"],
            "Completeness": best["Completeness"],
            "Contamination": best["Contamination"],
            "has_full_rrna_set": best.get("has_full_rrna_set", False),
            "has_mimag_trna_set": best.get("has_mimag_trna_set", False),
            "trna_aa_types": best.get("trna_aa_types", 0),
            "classification": best.get("classification", ""),
            "phylum": best.get("phylum", ""),
            "family": best.get("family", ""),
            "genus": best.get("genus", ""),
            "species": best.get("species", ""),
        })

    group_df = pd.DataFrame(group_rows)
    group_df.to_csv(OUT / "genome_clusters.tsv", sep="\t", index=False)

    taxvamb_specific = group_df[
        group_df["good_binners"].apply(
            lambda x: "taxvamb_taxometer" in x.split(",") and "vamb" not in x.split(",")
        )
    ].copy()

    taxvamb_specific.to_csv(
        OUT / "table3_taxvamb_specific_clusters.tsv",
        sep="\t",
        index=False,
    )

    vamb_only = group_df[
        group_df["good_binners"].apply(
            lambda x: "vamb" in x.split(",") and "taxvamb_taxometer" not in x.split(",")
        )
    ].copy()

    vamb_only.to_csv(
        OUT / "table4_vamb_only_clusters.tsv",
        sep="\t",
        index=False,
    )

    print("[OK] wrote summary tables to", OUT)


if __name__ == "__main__":
    main()
