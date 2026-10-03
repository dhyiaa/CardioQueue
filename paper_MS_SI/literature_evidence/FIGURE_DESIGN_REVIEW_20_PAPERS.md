# Figure design review: 20 related studies

Twenty primary papers and their main or supplementary figures were reviewed on 27 September 2026. The audit focused on study-flow diagrams, score distributions, abstention or threshold displays, calibration and workload plots, structure visualizations, and the boundary between computation and clinical interpretation.

| Papers | Figure pattern retained for CardioQueue |
|---|---|
| [CardioBoost](https://doi.org/10.1038/s41436-020-00972-3), [CardioClassifier](https://pmc.ncbi.nlm.nih.gov/articles/PMC6558251/), [CardioVAI](https://doi.org/10.1002/humu.23665), [CVD-PP](https://doi.org/10.1161/CIRCGEN.123.004464) | Give thresholds a visual identity, show the automated-to-expert boundary, and connect performance to a clinical workload rather than AUROC alone. |
| [GENESIS](https://pmc.ncbi.nlm.nih.gov/articles/PMC9018586/), [CardioVar](https://doi.org/10.1093/bioadv/vbag135), [AlphaMissense](https://doi.org/10.1126/science.adg7492), [EVE](https://doi.org/10.1038/s41586-021-04043-8) | Separate development, evaluation, and VUS deployment; show abstention explicitly; retain exact cohort sizes; avoid overcrowded grids of nearly identical curves. |
| [REVEL](https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/), [M-CAP](https://doi.org/10.1038/ng.3703), [ClinPred](https://pmc.ncbi.nlm.nih.gov/articles/PMC6174354/), [MVP](https://pmc.ncbi.nlm.nih.gov/articles/PMC7820281/) | Pair score distributions with operating thresholds, report candidate-list reduction or review burden, and show ablation effects around a zero-reference line when incremental value is the question. |
| [VARITY](https://pmc.ncbi.nlm.nih.gov/articles/PMC8546039/), [AlphScore](https://doi.org/10.1093/bioinformatics/btad280), [SIGMA](https://pubmed.ncbi.nlm.nih.gov/38211594/), [GeneTerpret](https://doi.org/10.1186/s12920-022-01166-3) | Use decision-relevant operating points, group structure features by family, report negative structure ablations honestly, and end with a clinician-controlled review step. |
| [TGex](https://pmc.ncbi.nlm.nih.gov/articles/PMC6937949/), [InterVar](https://pmc.ncbi.nlm.nih.gov/articles/PMC5294755/), [Kroncke et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6383132/), [Suay-Corredera et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8260873/) | Keep workflows staged and compact, never draw score-to-classification arrows, and use structure graphics only when they answer a measured mechanistic or performance question. |

## Changes implemented

1. Figure 1 now follows six numbered stages: registry assembly, pre-fit quarantine, model fitting, held-out evaluation, robustness tests, and review-order translation. Exact current counts replace the obsolete 430-row headline.
2. The clinical boundary is sequential: score and review priority, then expert ACMG/AMP adjudication, then amended classification or management only if independently supported. A clinician override bypasses score order for urgent or management-sensitive cases.
3. Figure 4 now uses the CardioBoost threshold precedent but avoids CardioBoost's clinical class wording. It shows the 775 held-out B/LB and P/LP score distributions, complete fixed-zone outcomes, and VUS workload by source.
4. The x-axis is a CardioQueue score, not a patient-level or ACMG pathogenicity probability. The outer zones are lower-priority and accelerated *review* zones; the center is explicit deferral.
5. FoldX remains a quantitative negative ablation. No decorative protein rendering was promoted to a main figure because numeric DDG did not independently improve performance over mapping eligibility.

## Deferred improvements

The remaining legacy plots should eventually adopt the same typography, colorblind-safe palette, line-style redundancy, and vector export. A source-specific forest plot and a threshold coverage-versus-performance curve are useful future additions, but were not inserted merely to increase panel count.
