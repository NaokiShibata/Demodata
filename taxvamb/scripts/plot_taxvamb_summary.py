from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


SUMMARY = Path("results/summary")
FIGDIR = Path("results/figures")
FIGDIR.mkdir(parents=True, exist_ok=True)


def savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def plot_quality_by_binner() -> None:
    df = pd.read_csv(SUMMARY / "table1_quality_by_binner.tsv", sep="\t")
    df = df.set_index("binner")

    cols = ["HQ", "MQ", "LQ", "contaminated"]
    ax = df[cols].plot(kind="bar", stacked=True, figsize=(8, 5))

    ax.set_xlabel("Binner")
    ax.set_ylabel("Number of bins")
    ax.set_title("MAG quality by binner")
    ax.legend(title="Quality", bbox_to_anchor=(1.02, 1), loc="upper left")

    savefig(FIGDIR / "fig1_quality_by_binner.png")


def plot_taxonomic_diversity() -> None:
    df = pd.read_csv(SUMMARY / "table2_taxonomic_diversity_by_binner.tsv", sep="\t")
    df = df.set_index("binner")

    cols = ["phyla", "classes", "orders", "families", "genera", "species"]
    ax = df[cols].plot(kind="bar", figsize=(9, 5))

    ax.set_xlabel("Binner")
    ax.set_ylabel("Number of unique taxa")
    ax.set_title("Taxonomic diversity of HQ/MQ MAGs")
    ax.legend(title="Rank", bbox_to_anchor=(1.02, 1), loc="upper left")

    savefig(FIGDIR / "fig2_taxonomic_diversity.png")


def plot_rrna_trna_support() -> None:
    df = pd.read_csv(SUMMARY / "table1_quality_by_binner.tsv", sep="\t")
    df = df.set_index("binner")

    cols = ["bins_with_full_rrna_set", "bins_with_mimag_trna_set"]
    ax = df[cols].plot(kind="bar", figsize=(8, 5))

    ax.set_xlabel("Binner")
    ax.set_ylabel("Number of bins")
    ax.set_title("rRNA/tRNA support by binner")
    ax.legend(
        ["5S/16S/23S rRNA detected", ">=18 tRNA isotypes detected"],
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    savefig(FIGDIR / "fig3_rrna_trna_support.png")


def plot_same_mag_group_patterns() -> None:
    df = pd.read_csv(SUMMARY / "genome_clusters.tsv", sep="\t")

    def simplify_pattern(x) -> str:
        if pd.isna(x) or str(x).strip() == "":
            return "No HQ/MQ recovery"

        items = sorted([i for i in str(x).split(",") if i and i != "nan"])
        return "+".join(items) if items else "No HQ/MQ recovery"

    df["pattern"] = df["good_binners"].apply(simplify_pattern)
    counts = df["pattern"].value_counts().sort_values(ascending=False).head(15)

    ax = counts.plot(kind="bar", figsize=(10, 6))

    ax.set_xlabel("Binners with HQ/MQ bin in the same MAG group")
    ax.set_ylabel("Number of same-MAG groups")
    ax.set_title("Patterns of HQ/MQ recovery across binners")

    savefig(FIGDIR / "fig4_same_mag_group_patterns.png")


def main() -> None:
    plot_quality_by_binner()
    plot_taxonomic_diversity()
    plot_rrna_trna_support()
    plot_same_mag_group_patterns()

    print("[OK] wrote figures to", FIGDIR)


if __name__ == "__main__":
    main()
