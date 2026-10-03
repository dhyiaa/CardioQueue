# gnomAD v4 External Inputs

This folder stores small inputs for targeted gnomAD queries.

Expected files:

- Reviewed gene panel, one gene symbol per line.
- Variant query tables with either:
  - `variant_id`, formatted as `chrom-pos-ref-alt`, or
  - `chrom`, `pos`, `ref`, `alt` columns.

The example gene list is only a seed for testing scripts. It is not yet the
final study panel.

