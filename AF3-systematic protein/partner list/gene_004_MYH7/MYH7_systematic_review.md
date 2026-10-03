# Gene 004: MYH7 AF3 Systematic Review

Date: 2026-07-06

## Current Decision

MYH7 is gene #4 in the systematic AF3 queue. It is one of the most important cardiomyopathy genes in the project, but it should not be handled as one full-length thick-filament AF3 job.

Reason:

```text
MYH7 canonical protein length = 1,935 aa
```

MYH7 is a modular sarcomeric motor protein. The top-3 AF3 policy is especially useful here because many sarcomere partners are biologically plausible. For this gene packet, the default upload set is limited to three nonredundant mechanisms:

1. Actomyosin motor interface.
2. Lever-arm/light-chain context.
3. MYBPC3/S2 regulatory context.

## Local Dataset Context

| Metric | Value |
|---|---:|
| Included MYH7 rows from priority file | 5228 |
| Current matrix MYH7 rows | 5228 |
| Binary rows | 2333 |
| VUS rows | 2869 |
| Pathogenic rows | 423 |
| Benign rows | 1910 |
| Current rows with `protein_position` | 81 |
| Current rows with VEP protein position-like value | 4208 |
| Current rows with dbNSFP protein position raw | 2971 |

Important caveat: MYH7 has many transcript/protein-position representations in dbNSFP, including `-1` values and long semicolon-separated position lists. The AF3 mapper must reconcile these to UniProt P12883 before feature claims.

## Sequence Audit

| Field | Value |
|---|---|
| UniProt accession | P12883 |
| Reviewed entry | yes |
| Canonical protein length | 1,935 aa |
| Baseline sequence policy | reviewed UniProt canonical MYH7 amino-acid sequence |
| PTM/isoform policy | no PTMs or alternative isoforms in baseline; separate sensitivity jobs only if justified |
| Full-length AF3 | not recommended as first systematic design |

## Candidate Partner Summary

| Partner | Priority | Why |
|---|---|---|
| ACTC1/cardiac actin | top-3 | MYH7 is an actin-based ATPase motor and UniProt maps actin-binding regions within the motor domain |
| MYL2 + MYL3 | top-3 | MYH7 myosin is a hexamer with regulatory and essential light chains; IQ/neck region is compact |
| MYBPC3 | top-3 | Cardiac myosin-binding protein C binds MHC/F-actin/thin filaments and regulates the crossbridge region |
| MYH7 homodimer/rod | hold | Important for thick-filament assembly but less specific for top-3 interaction mapping |
| TPM1/thin filament expansion | hold | Biologically relevant but too large/redundant for the MYH7 top-3 packet |

Detailed partner evidence table:

```text
03_candidate_partner_evidence.tsv
```

## Proposed AF3 Designs

| Design | Residues | Partner | Status |
|---|---:|---|---|
| Motor-ACTC1 | 85-778 | ACTC1 | top-3 draft JSON ready |
| Neck-light-chains | 760-840 | MYL2 + MYL3 | top-3 draft JSON ready |
| S2-MYBPC3-Nterm | 838-963 | MYBPC3 1-452 | top-3 draft JSON ready |
| Rod homodimer | 838-1935/subfragment | MYH7 copy | held, not generated |
| Expanded thin filament | motor domain | ACTC1/TPM1/TNNI3/TNNT2 | held, not generated |

Detailed design table:

```text
04_af3_design_decisions.tsv
```

## Source Links For Manual Review

- UniProt MYH7: https://www.uniprot.org/uniprotkb/P12883/entry
- UniProt ACTC1: https://www.uniprot.org/uniprotkb/P68032/entry
- UniProt MYL2: https://www.uniprot.org/uniprotkb/P10916/entry
- UniProt MYL3: https://www.uniprot.org/uniprotkb/P08590/entry
- UniProt MYBPC3: https://www.uniprot.org/uniprotkb/Q14896/entry
- Human beta-cardiac myosin motor paper: https://pubmed.ncbi.nlm.nih.gov/26246073/
- Human cardiac beta-myosin S2-delta structure paper: https://pubmed.ncbi.nlm.nih.gov/17095604/

Specific source/citation table for Ctrl-F review:

```text
09_specific_citation_sources.tsv
```

## AF3 Draft Job JSONs

Top-3 upload candidate:

```text
af3_job_json/AF3_SYS_gene_004_MYH7_TOP3_upload_alphafoldserver.json
```

Manifest:

```text
af3_job_json/AF3_SYS_gene_004_MYH7_json_manifest.tsv
```

Jobs intentionally held/not generated:

```text
af3_job_json/AF3_SYS_gene_004_MYH7_not_generated_jobs.tsv
```

| AF3 Job | Entities | Upload Status | Evidence IDs |
|---|---|---|---|
| `AF3_SYS_G004_MYH7_MOTOR_85_778_ACTC1_DRAFT` | MYH7 85-778 + ACTC1 | top-3 | `SRC_UNIPROT_MYH7_MOTOR`, `SRC_UNIPROT_MYH7_ACTIN_REGIONS` |
| `AF3_SYS_G004_MYH7_NECK_760_840_MYL2_MYL3_DRAFT` | MYH7 760-840 + MYL2 + MYL3 | top-3 | `SRC_UNIPROT_MYH7_SUBUNIT`, `SRC_UNIPROT_MYH7_IQ` |
| `AF3_SYS_G004_MYH7_S2_838_963_MYBPC3_NTERM_1_452_DRAFT` | MYH7 838-963 + MYBPC3 1-452 | top-3 | `SRC_UNIPROT_MYBPC3_FUNCTION`, `SRC_MYH7_S2_PMID17095604` |

## What To Watch When AF3 Results Come Back

For MYH7, useful AF3 results should be interpreted by mechanism. The motor-ACTC1 job should be treated as a simplified actomyosin interface because native actin is filamentous. The light-chain job should show stable contacts around the MYH7 IQ/neck region. The MYBPC3 job is exploratory but important; it should only be promoted if pair confidence and interface PAE support a plausible MYH7 S2/MYBPC3 N-terminal interaction.

## Next Action For MYH7

1. Upload only the top-3 JSON unless you explicitly want to spend quota on held jobs.
2. Place AlphaFold Server output folders under:

```text
af3_job_results/
```

3. After output is added, run pair-level ipTM/PAE/contact QC and then map interface features to UniProt-reconciled MYH7 variant positions.
