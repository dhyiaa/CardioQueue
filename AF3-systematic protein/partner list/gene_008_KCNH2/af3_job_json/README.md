# KCNH2 AF3 Job JSONs

Default upload file:

```text
AF3_SYS_gene_008_KCNH2_TOP3_upload_alphafoldserver.json
```

This file contains three jobs:

1. KCNH2 full-length homotetramer.
2. KCNH2 PAS/PAC 1-150 plus KCNH2 CNBHD 720-870 regulatory-domain pair.
3. KCNH2 full length plus KCNE2 full length comparator.

Upload order matters. The KCNE2 job is biologically supported but previously weak in AF3, so it is a comparator, not the lead feature source. Held jobs are documented in `AF3_SYS_gene_008_KCNH2_not_generated_jobs.tsv`.
