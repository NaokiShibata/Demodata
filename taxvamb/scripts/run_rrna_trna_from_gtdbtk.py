from __future__ import annotations

import subprocess
from pathlib import Path
import pandas as pd


GENOME_DIR = Path("results/all_bins")
GTDBTK_DIR = Path("results/gtdbtk_all")
OUTDIR = Path("results/rna_genes")

BARRNAP_DIR = OUTDIR / "barrnap"
TRNASCAN_DIR = OUTDIR / "trnascan"

BARRNAP_DIR.mkdir(parents=True, exist_ok=True)
TRNASCAN_DIR.mkdir(parents=True, exist_ok=True)

STANDARD_AA = {
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Gln",
    "Glu",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val",
}


def read_gtdbtk_domains() -> pd.DataFrame:
    summary_files = list(GTDBTK_DIR.rglob("*.summary.tsv"))

    if not summary_files:
        raise SystemExit(f"No GTDB-Tk summary files found in {GTDBTK_DIR}")

    rows = []

    for path in summary_files:
        df = pd.read_csv(path, sep="\t")

        if "user_genome" not in df.columns or "classification" not in df.columns:
            continue

        for _, row in df.iterrows():
            bin_id = str(row["user_genome"])
            classification = str(row["classification"])

            if classification.startswith("d__Bacteria"):
                domain = "Bacteria"
                barrnap_kingdom = "bac"
                trnascan_mode = "-B"
            elif classification.startswith("d__Archaea"):
                domain = "Archaea"
                barrnap_kingdom = "arc"
                trnascan_mode = "-A"
            else:
                domain = "unknown"
                barrnap_kingdom = ""
                trnascan_mode = ""

            rows.append(
                {
                    "bin_id": bin_id,
                    "classification": classification,
                    "domain": domain,
                    "barrnap_kingdom": barrnap_kingdom,
                    "trnascan_mode": trnascan_mode,
                }
            )

    if not rows:
        raise SystemExit("No usable GTDB-Tk records found")

    return pd.DataFrame(rows).drop_duplicates("bin_id")


def parse_barrnap_gff(path: Path) -> dict[str, int]:
    counts = {
        "rrna_5s": 0,
        "rrna_16s": 0,
        "rrna_23s": 0,
    }

    if not path.exists():
        return counts

    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attrs = parts[8].lower()

            if feature_type != "rRNA":
                continue

            if "16s" in attrs:
                counts["rrna_16s"] += 1
            elif "23s" in attrs:
                counts["rrna_23s"] += 1
            elif "5s" in attrs and "5.8s" not in attrs:
                counts["rrna_5s"] += 1

    return counts


def parse_trnascan(path: Path) -> dict[str, object]:
    aa_types = set()
    trna_count = 0

    if not path.exists():
        return {
            "trna_count": 0,
            "trna_aa_types": 0,
            "trna_aa_list": "",
        }

    with path.open() as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith("Sequence") or line.startswith("Name"):
                continue

            if line.startswith("-"):
                continue

            parts = line.split()

            # tRNAscan-SE tabular output:
            # Name tRNA# Begin End Type Codon ...
            if len(parts) < 5:
                continue

            trna_type = parts[4]

            if trna_type in {"Pseudo", "Undet", "Sup", "SeC"}:
                continue

            trna_count += 1

            if trna_type in STANDARD_AA:
                aa_types.add(trna_type)

    return {
        "trna_count": trna_count,
        "trna_aa_types": len(aa_types),
        "trna_aa_list": ",".join(sorted(aa_types)),
    }


def run_for_bin(row: pd.Series) -> dict[str, object]:
    bin_id = row["bin_id"]
    domain = row["domain"]
    barrnap_kingdom = row["barrnap_kingdom"]
    trnascan_mode = row["trnascan_mode"]

    genome = GENOME_DIR / f"{bin_id}.fna"

    if not genome.exists():
        print(f"[WARN] genome not found: {genome}")
        return {
            "bin_id": bin_id,
            "domain": domain,
            "rna_mode": "missing_genome",
            "rrna_5s": 0,
            "rrna_16s": 0,
            "rrna_23s": 0,
            "trna_count": 0,
            "trna_aa_types": 0,
            "trna_aa_list": "",
        }

    if domain not in {"Bacteria", "Archaea"}:
        print(f"[WARN] unsupported or unknown domain for {bin_id}: {domain}")
        return {
            "bin_id": bin_id,
            "domain": domain,
            "rna_mode": "skipped_unknown_domain",
            "rrna_5s": 0,
            "rrna_16s": 0,
            "rrna_23s": 0,
            "trna_count": 0,
            "trna_aa_types": 0,
            "trna_aa_list": "",
        }

    print(f"[INFO] {bin_id}: {domain}")

    barrnap_out = BARRNAP_DIR / f"{bin_id}.gff"
    trnascan_out = TRNASCAN_DIR / f"{bin_id}.tsv"
    trnascan_stdout = TRNASCAN_DIR / f"{bin_id}.stdout"
    trnascan_stderr = TRNASCAN_DIR / f"{bin_id}.stderr"

    with barrnap_out.open("w") as out:
        subprocess.run(
            [
                "pixi",
                "run",
                "barrnap",
                "--kingdom",
                barrnap_kingdom,
                str(genome),
            ],
            check=True,
            stdout=out,
        )

    with trnascan_stdout.open("w") as out, trnascan_stderr.open("w") as err:
        subprocess.run(
            [
                "pixi",
                "run",
                "tRNAscan-SE",
                trnascan_mode,
                "-o",
                str(trnascan_out),
                str(genome),
            ],
            check=True,
            stdout=out,
            stderr=err,
        )

    rrna = parse_barrnap_gff(barrnap_out)
    trna = parse_trnascan(trnascan_out)

    return {
        "bin_id": bin_id,
        "domain": domain,
        "rna_mode": f"{barrnap_kingdom}/{trnascan_mode}",
        **rrna,
        **trna,
    }


def main() -> None:
    domains = read_gtdbtk_domains()
    domains.to_csv(OUTDIR / "bin_domain_from_gtdbtk.tsv", sep="\t", index=False)

    records = []

    for _, row in domains.iterrows():
        records.append(run_for_bin(row))

    summary = pd.DataFrame(records)

    summary["has_5s"] = summary["rrna_5s"] > 0
    summary["has_16s"] = summary["rrna_16s"] > 0
    summary["has_23s"] = summary["rrna_23s"] > 0
    summary["has_full_rrna_set"] = (
        summary["has_5s"] & summary["has_16s"] & summary["has_23s"]
    )
    summary["has_mimag_trna_set"] = summary["trna_aa_types"] >= 18

    summary.to_csv(OUTDIR / "rna_summary.tsv", sep="\t", index=False)

    print("[OK] wrote", OUTDIR / "rna_summary.tsv")


if __name__ == "__main__":
    main()
