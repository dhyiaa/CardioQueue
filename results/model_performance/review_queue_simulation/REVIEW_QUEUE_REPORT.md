# Held-out variant review-queue simulation

This retrospective simulation ordered 430 source-held-out variants by the frozen CardioQueue score. The comparison queue used 10,000 random permutations as a proxy for arrival-order review when arrival order is unrelated to pathogenicity. No observed timestamps, staffing data, turnaround times, or patient outcomes were available.

At a 20% review budget, the score-ranked queue reviewed 86 variants and recovered 86 P/LP variants, corresponding to 33.7% recall and 100.0% precision. At a 50% budget, it recovered 211 of 255 P/LP variants (82.7%) with 98.1% precision.

Recovering 80% of P/LP variants required 207 score-ranked reviews, compared with a median 344 reviews under random ordering. The observed workload reduction was 39.8%. At 95% recovery, the corresponding counts were 255 and 410, a 37.8% reduction.

On matched observations, CardioQueue recovered 72.8% of P/LP observations in the first half of the queue, compared with 66.3% for official CardioBoost. Reaching 90% recovery required 155 versus 185 observations.

These are worklist-efficiency estimates on a consensus-filtered retrospective cohort. They do not measure emergency triage, diagnostic yield, clinical action, or turnaround time.
