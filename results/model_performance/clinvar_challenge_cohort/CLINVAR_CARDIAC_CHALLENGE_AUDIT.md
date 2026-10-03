# ClinVar Cardiac Challenge Eligibility Audit

This audit found **12,341** binary, at-least-two-star ClinVar rows with either direct cardiac-VCEP verification or a cardiac phenotype keyword without an obvious RASopathy keyword, no HiRO/eMERGE/CardioBoost overlap, and an existing development assignment.

- Three-star expert-panel rows: **113** (53 B/LB; 60 P/LP)
- Cardiac-VCEP-verified rows with concordant current labels: **92** (52 B/LB; 40 P/LP)
- Two-star multiple-submitter/no-conflict rows: **12,228**
- Existing source-held-out cohort: **430** (175 B/LB; 255 P/LP)

## Defensible design

Keep the 430-row HiRO/eMERGE/CardioBoost result as the primary source-held-out analysis. Build a separate 570-row ClinVar challenge cohort only after assertion-level disease-context checks. Use the locally verified cardiac-VCEP rows as a seed, then draw the remainder from condition-filtered high-confidence records using a prespecified gene-by-consequence design. A provisional target of 325 B/LB and 245 P/LP would produce a descriptive pooled benchmark of exactly 500 B/LB and 500 P/LP, but the two cohorts must also be reported separately.

## Important boundary

ClinVar review stars describe review status, not clinical specialty. The phenotype keyword filter improves relevance but cannot prove that the reviewing panel was a cardiogenetics panel. The manifest therefore remains an eligibility audit and does not change model splits.
