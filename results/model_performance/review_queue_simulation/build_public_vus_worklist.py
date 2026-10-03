#!/usr/bin/env python3
"""Build a public-source, research-only VUS evidence-review worklist."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
INPUT = (
    ROOT
    / "results/model_performance/final_cardioqueue_three_zone_all_rows"
    / "all_registry_scores.tsv"
)
OUT = Path(__file__).resolve().parent


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(INPUT, sep="\t", low_memory=False)
    public_source = frame[["in_clinvar", "in_emerge", "in_cardioboost"]].any(axis=1)
    keep = (
        frame["model_label_3class"].eq("VUS")
        & frame["primary_model_inclusion"].eq("include")
        & frame["in_hiro"].ne(True)
        & public_source
    )
    worklist = frame.loc[
        keep,
        [
            "variant_id",
            "primary_gene",
            "modeling_scope",
            "cardioqueue_score",
            "review_zone",
            "in_clinvar",
            "in_emerge",
            "in_cardioboost",
        ],
    ].copy()
    worklist = worklist.sort_values(
        ["cardioqueue_score", "variant_id"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    worklist.insert(0, "research_review_rank", worklist.index + 1)
    worklist.insert(1, "queue_fraction", (worklist.index + 1) / len(worklist))
    worklist["within_gene_rank"] = (
        worklist.groupby("primary_gene")["cardioqueue_score"]
        .rank(method="first", ascending=False)
        .astype("Int64")
    )
    worklist["use_constraint"] = "Research evidence-review priority; not a clinical classification"
    if worklist["variant_id"].duplicated().any():
        raise ValueError("Public VUS worklist must contain one row per variant_id")
    if worklist.empty:
        raise ValueError("Public-source VUS worklist is empty")
    worklist.to_csv(OUT / "public_source_vus_research_worklist.tsv", sep="\t", index=False)

    summary = {
        "input": str(INPUT.relative_to(ROOT)),
        "rows": len(worklist),
        "genes": int(worklist["primary_gene"].nunique()),
        "accelerated_review": int(worklist["review_zone"].eq("accelerated").sum()),
        "deferred_review": int(worklist["review_zone"].eq("deferred").sum()),
        "lower_priority_review": int(worklist["review_zone"].eq("lower_priority").sum()),
        "excludes_hiro_linked_rows": True,
        "required_use_statement": "Research evidence-review priority; not a clinical classification",
    }
    (OUT / "public_source_vus_research_worklist.summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
