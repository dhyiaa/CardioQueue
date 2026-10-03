# Gene 001 RYR2 AF3 Result Review

Date: 2026-07-06

Batch reviewed:

```text
af3_job_results/folds_2026_07_06_02_03/
```

This review uses the same support thresholds used in the AF3Draft model work:

| Support class | Criteria |
|---|---|
| Strong | pair ipTM >= 0.70, pair min PAE <= 5, max contact probability >= 0.60 |
| Moderate | pair ipTM >= 0.45, pair min PAE <= 10, max contact probability >= 0.35 |
| Unsupported | Does not pass moderate threshold |

For each AF3 job, the best-ranked model was selected by `ranking_score`.

## Main Result

| Job | Main pair | pair ipTM | min PAE | max contact | Support | Decision |
|---|---|---:|---:|---:|---|---|
| `RYR2_CaMBD1_1940_1965_CALM1_CA4` | RYR2 fragment-CALM1 | 0.50 | 3.97 | 0.78 | Moderate | Use with sensitivity flag |
| `RYR2_CaMBD2_3580_3611_CALM1_CA4` | RYR2 fragment-CALM1 | 0.84 | 1.00 | 0.99 | Strong | Promote to trusted AF3Draft feature |
| `RYR2_CaMBD3_4246_4275_CALM1_CA5` | RYR2 fragment-CALM1 | 0.44 | 4.27 | 0.61 | Unsupported by threshold | Borderline, review/redesign |
| `RYR2_LL1_4521_4573_CASQ2` | RYR2 LL1-CASQ2 | 0.16 | 13.99 | 0.26 | Unsupported | Do not promote |
| `CASQ2_TRDN_FULL_LENGTH` | CASQ2-TRDN | 0.21 | 15.87 | 0.16 | Unsupported | Do not promote |

The strongest result is:

```text
RYR2 CaMBD2 3580-3611 + CALM1 + Ca2+
```

This is the only RYR2 AF3 job from this batch that should be promoted as a trusted interaction feature without caveat.

## Job-Level QC

| Job | Ranking | ipTM | pTM | Fraction disordered | Best protein-pair support |
|---|---:|---:|---:|---:|---|
| `CASQ2_TRDN_FULL_LENGTH` | 0.58 | 0.21 | 0.38 | 0.68 | Unsupported |
| `RYR2_CaMBD1_1940_1965_CALM1_CA4` | 0.66 | 0.52 | 0.72 | 0.19 | Moderate |
| `RYR2_CaMBD2_3580_3611_CALM1_CA4` | 0.94 | 0.85 | 0.85 | 0.18 | Strong |
| `RYR2_CaMBD3_4246_4275_CALM1_CA5` | 0.61 | 0.45 | 0.71 | 0.21 | Unsupported |
| `RYR2_LL1_4521_4573_CASQ2` | 0.40 | 0.16 | 0.73 | 0.25 | Unsupported |

## Interpretation

### CaMBD2-CALM1

This is a clean positive result.

```text
pair ipTM = 0.84
min PAE = 1.00
max contact probability = 0.99
```

It supports using RYR2 CaMBD2/CALM1 interface features in the AF3Draft layer after variant-position mapping is improved. This is also biologically aligned with the UniProt CALM-interaction annotation around RYR2 residues 3581-3610.

### CaMBD1-CALM1

This is a moderate result.

```text
pair ipTM = 0.50
min PAE = 3.97
max contact probability = 0.78
```

It passes the moderate threshold, but not the strong threshold. It can be used as a secondary AF3Draft feature with a sensitivity flag.

### CaMBD3-CALM1

This result is biologically interesting but fails the formal threshold.

```text
pair ipTM = 0.44
min PAE = 4.27
max contact probability = 0.61
```

The PAE and contact probability look acceptable, but pair ipTM is just below the moderate cutoff of 0.45. Because this job was motivated by RCSB `7KL5`, it should be manually reviewed or redesigned before excluding the biology.

Recommended next action:

```text
Compare the AF3 CaMBD3 output against 7KL5 and consider rerunning with 4 Ca2+ instead of 5 Ca2+ or adjusted flanking residues.
```

### RYR2 LL1-CASQ2

This does not support a trusted interface.

```text
pair ipTM = 0.16
min PAE = 13.99
max contact probability = 0.26
```

The isolated RYR2 luminal loop may not be enough to recreate the CASQ2 interaction, or the interaction may require additional luminal/SR context. Do not use this output as an interface feature.

### CASQ2-TRDN Full-Length

This does not support a trusted interface.

```text
pair ipTM = 0.21
min PAE = 15.87
max contact probability = 0.16
fraction disordered = 0.68
```

The high disorder and weak pair metrics suggest the full-length protein design is not suitable. Redesign with a TRDN luminal/KEKE fragment and/or calcium context if this biology remains important.

## Feature Promotion Decision

| Job | Feature decision |
|---|---|
| `RYR2_CaMBD2_3580_3611_CALM1_CA4` | Promote to trusted AF3Draft interface feature |
| `RYR2_CaMBD1_1940_1965_CALM1_CA4` | Promote as moderate/sensitivity feature |
| `RYR2_CaMBD3_4246_4275_CALM1_CA5` | Do not promote yet; manual review/redesign |
| `RYR2_LL1_4521_4573_CASQ2` | Do not promote; redesign needed |
| `CASQ2_TRDN_FULL_LENGTH` | Do not promote; redesign needed |

## Output Tables

```text
gene_001_RYR2_af3_job_qc_summary.tsv
gene_001_RYR2_af3_pair_qc.tsv
gene_001_RYR2_af3_chain_qc.tsv
gene_001_RYR2_af3_feature_decisions.tsv
```

## Bottom Line

The RYR2 AF3 systematic run produced one strong and one moderate calmodulin-domain result. That is useful and biologically coherent.

It did not support the first CASQ2/TRDN luminal designs. Those should be redesigned rather than promoted.
