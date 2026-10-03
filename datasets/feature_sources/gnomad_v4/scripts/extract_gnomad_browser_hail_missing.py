#!/usr/bin/env python3
"""Extract missing variants from the gnomAD v4.1 browser Hail Table.

This script is intentionally written as a Hail job wrapper. It requires:

- Java runtime
- Python `hail`
- Google Cloud credentials with a billing project, because the official browser
  table is in a requester-pays bucket

It exports one row per requested variant and stores the matched browser row as
JSON so downstream parsing can be schema-aware after a pilot run.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


GNOMAD_BROWSER_HT = "gs://gcp-public-data--gnomad/release/4.1/ht/browser/gnomad.browser.v4.1.sites.ht"
DEFAULT_GCS_CONNECTOR_JAR = (
    "datasets/feature_sources/gnomad_v4/raw/browser_hail/jars/"
    "gcs-connector-hadoop3-2.2.31-shaded.jar"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp-dir", default="tmp/hail/gnomad_browser_missing")
    parser.add_argument("--gcp-project", help="Google Cloud billing project for requester-pays access.")
    parser.add_argument("--hail-log", default="datasets/feature_sources/gnomad_v4/interim/gnomad_browser_hail.log")
    parser.add_argument("--browser-ht", default=GNOMAD_BROWSER_HT)
    parser.add_argument(
        "--gcs-connector-jar",
        default=DEFAULT_GCS_CONNECTOR_JAR,
        type=Path,
        help="Local shaded Google Cloud Storage connector jar for Spark gs:// access.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gcs_connector_jar = args.gcs_connector_jar.resolve()
    if not gcs_connector_jar.exists():
        raise FileNotFoundError(f"Missing GCS connector jar: {gcs_connector_jar}")
    adc_path = Path.home() / ".config/gcloud/application_default_credentials.json"
    service_account_key_path = Path.home() / ".config/gcloud/gnomad-hail-local-key.json"
    if adc_path.exists():
        os.environ.setdefault("GOOGLE_APPLICATION_CREDENTIALS", str(adc_path))

    import hail as hl

    spark_conf = {
        "spark.jars": str(gcs_connector_jar),
        "spark.driver.extraClassPath": str(gcs_connector_jar),
        "spark.executor.extraClassPath": str(gcs_connector_jar),
        "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
        "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",
        "spark.hadoop.fs.gs.auth.type": "APPLICATION_DEFAULT",
        "spark.hadoop.google.cloud.auth.type": "APPLICATION_DEFAULT",
        "spark.hadoop.fs.gs.client.type": "HTTP_API_CLIENT",
    }
    if service_account_key_path.exists():
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(service_account_key_path)
        spark_conf.update(
            {
                "spark.driverEnv.GOOGLE_APPLICATION_CREDENTIALS": str(service_account_key_path),
                "spark.executorEnv.GOOGLE_APPLICATION_CREDENTIALS": str(service_account_key_path),
                "spark.hadoop.fs.gs.auth.type": "SERVICE_ACCOUNT_JSON_KEYFILE",
                "spark.hadoop.fs.gs.auth.service.account.enable": "true",
                "spark.hadoop.fs.gs.auth.service.account.json.keyfile": str(service_account_key_path),
                "spark.hadoop.google.cloud.auth.service.account.enable": "true",
                "spark.hadoop.google.cloud.auth.service.account.json.keyfile": str(service_account_key_path),
            }
        )
    elif adc_path.exists():
        spark_conf.update(
            {
                "spark.driverEnv.GOOGLE_APPLICATION_CREDENTIALS": str(adc_path),
                "spark.executorEnv.GOOGLE_APPLICATION_CREDENTIALS": str(adc_path),
                "spark.hadoop.fs.gs.auth.service.account.enable": "false",
                "spark.hadoop.google.cloud.auth.service.account.enable": "false",
            }
        )
        with adc_path.open() as handle:
            adc = json.load(handle)
        if {"client_id", "client_secret", "refresh_token"}.issubset(adc):
            spark_conf.update(
                {
                    "spark.hadoop.google.cloud.auth.client.id": adc["client_id"],
                    "spark.hadoop.google.cloud.auth.client.secret": adc["client_secret"],
                    "spark.hadoop.google.cloud.auth.refresh.token": adc["refresh_token"],
                    "spark.hadoop.fs.gs.auth.client.id": adc["client_id"],
                    "spark.hadoop.fs.gs.auth.client.secret": adc["client_secret"],
                    "spark.hadoop.fs.gs.auth.refresh.token": adc["refresh_token"],
                }
            )
    if args.gcp_project:
        spark_conf.update(
            {
                "spark.hadoop.fs.gs.project.id": args.gcp_project,
                "spark.hadoop.fs.gs.requester.pays.mode": "ENABLED",
                "spark.hadoop.fs.gs.requester.pays.project.id": args.gcp_project,
            }
        )

    init_kwargs = {
        "tmp_dir": args.tmp_dir,
        "log": args.hail_log,
        "default_reference": "GRCh38",
        "spark_conf": spark_conf,
    }
    if args.gcp_project:
        init_kwargs["gcs_requester_pays_configuration"] = args.gcp_project
    hl.init(**init_kwargs)

    requested = hl.import_table(
        str(args.input),
        delimiter="\t",
        impute=False,
        types={"pos": hl.tint32},
        force_bgz=False,
    )
    chrom = requested.chrom
    chrom_hail = hl.if_else(
        chrom.startswith("chr"),
        chrom,
        hl.if_else((chrom == "MT") | (chrom == "M"), "chrM", "chr" + chrom),
    )
    requested = requested.annotate(
        locus=hl.locus(chrom_hail, requested.pos, reference_genome="GRCh38"),
        alleles=[requested.ref, requested.alt],
    ).key_by("locus", "alleles")

    browser = hl.read_table(args.browser_ht)
    browser = browser.key_by("locus", "alleles")

    annotated = requested.annotate(gnomad_browser_row=browser[requested.key])
    annotated = annotated.select(
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "genes",
        "sources",
        "labels_3class",
        gnomad_browser_status=hl.if_else(hl.is_defined(annotated.gnomad_browser_row), "ok", "not_found"),
        gnomad_browser_missing_reason=hl.if_else(
            hl.is_defined(annotated.gnomad_browser_row),
            "",
            "variant absent from gnomAD browser v4.1 table for exact GRCh38 locus+alleles",
        ),
        gnomad_browser_row_json=hl.if_else(
            hl.is_defined(annotated.gnomad_browser_row),
            hl.json(annotated.gnomad_browser_row),
            "",
        ),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    annotated.export(str(args.output))
    hl.stop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
