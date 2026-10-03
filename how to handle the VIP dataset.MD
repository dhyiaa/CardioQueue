Yes, **Neptune is worth downloading and testing**, but I would not treat it as a clean new independent 100k cardiogenetics dataset.

My read: **Neptune/VIP probably does contain cardiogenetics genes with ACMG/AMP-based clinical classifications**, including variants from the eMERGE reportable-gene panel and BCM HeartCare cardiac panel. But it will likely have **substantial overlap with your current ClinVar and eMERGE sources**, so its main value is not simply “more N.” Its value is **source concordance, discordance, and possibly some additional HeartCare/HGSC/internal-curated cardiac variants**.

## What Neptune/VIP actually is

The Neptune paper describes two things:

1. **Neptune software**, which manages reporting, EHR integration, reanalysis, FHIR/JSON/HTML reports, etc.
2. **VIP database**, which is the part you care about.

The VIP database contains variant position, allele, frequency, transcript data, gene annotations, disease association, inheritance, PubMed IDs, clinical-site comments, and categories. The paper reports **381,564 variants** in the VIP database. It was seeded by the two clinical reporting labs for **eMERGE III**, then updated with novel HGSC clinical lab variants and public variant resources. Importantly, it draws from **ClinVar, OMIM, literature review, and internal datasets**. 

The manuscript explicitly says variant filtering and interpretation for eMERGE and HeartCare were performed outside Neptune using **ACMG/AMP guidelines**, with ClinGen recommendations adopted when available, including MYH7 and CNV guidance. 

So: **yes, it has ACMG/AMP-based classifications**, but likely as final category labels such as B/LB, VUS, LP/P, risk allele, etc., not necessarily structured ACMG criteria like PM2, PP3, PS4 for every variant.

## Does it contain cardiogenetics genes?

Yes, almost certainly.

There are two relevant panels:

**eMERGE III:** Neptune was applied to a targeted gene panel with **68 consensus reportable genes**, and the paper says variants in these 68 genes were replayed through the VIP database as part of the review-burden analysis. 

**BCM HeartCare:** the HeartCare study added small variants and CNVs from **168 genes related to cardiac disease**, plus cardiovascular pharmacogenomics, Lp(a) risk alleles, and coronary artery disease PRS. 

That means Neptune is not just general Mendelian genetics. It has direct relevance to cardiogenetics, especially inherited cardiovascular disease and ACMG-actionable/reportable genes.

## Is it downloadable?

The paper says the VIP database is available for download through the Neptune GitLab repository and that Neptune interacts with a VIP snapshot in VCF format. 

I checked the current GitLab tree. The public repo still exists and includes a `data/vip` folder with two VIP snapshot files:

* `data/vip/2020-09-03.vip.txt`
* `data/vip/2020-10-13.vip.txt`

The public README also describes Neptune as taking genotype/coverage input, comparing it against a VIP database of known genetic variation, marking previously curated variants, selecting novel variants for review, and enabling reanalysis over time. ([GitLab][1]) ([GitLab][2])

So this is not theoretical. It is a real file source you should download and parse.

## Will it duplicate ClinVar/eMERGE?

Yes. I would assume **major duplication** until proven otherwise.

Reasons:

Neptune’s VIP was seeded by the eMERGE III clinical reporting laboratories, and your project already includes an eMERGE/Glazer arrhythmia source with 2,754 unique linked variants.  

The VIP also explicitly draws from **ClinVar**, and your model already uses ClinVar as the main public label source, with 83,915 primary post-QC ClinVar rows.  

Your current registry is already designed around one row per `chrom-pos-ref-alt` with source flags and source-membership tables, which is exactly the right structure for handling Neptune overlap. 

So if you add Neptune naïvely, you risk:

* duplicate variants under a new source name
* label leakage into train/test splits
* reusing eMERGE-derived labels as if they were independent
* inflating external validation because the same ClinVar/eMERGE labels appear through VIP
* mixing final clinical classifications with source-specific reporting categories

## But could Neptune still add useful new variants?

Yes, likely. The valuable part is the **non-overlapping subset**, especially from:

* BCM HeartCare’s 168 cardiac genes
* HGSC clinical lab variants added after eMERGE seeding
* internal curated variants not already in ClinVar
* VUS and likely benign variants from cardiac genes
* variants with literature/clinical-site comments/categories

But I would expect the useful new cardiogenetics subset to be **much smaller than 381,564**. The full VIP database is broad and includes many non-cardiac/reportable/pharmacogenomic/risk-allele records.

My realistic expectation:

| Category                                      |  Expected overlap with your current registry |
| --------------------------------------------- | -------------------------------------------: |
| ClinVar-derived VIP variants                  |                                    Very high |
| eMERGE 68-gene variants                       |                                    Very high |
| P/LP cardiac variants                         |                                 High overlap |
| B/LB/VUS cardiac variants                     |                     Moderate to high overlap |
| HeartCare-only/internal HGSC cardiac variants |               Potentially useful non-overlap |
| Risk alleles / PGx / PRS SNPs                 | Mostly exclude from your pathogenicity model |

## How I would use Neptune in your study

I would **not** immediately add Neptune into training.

I would first create a separate source audit:

### Step 1: Parse VIP into a raw source-record table

Keep it as:

```text
datasets/neptune/raw/neptune_vip_2020_10_13_source_records.tsv
```

Do not collapse immediately.

Fields to extract:

```text
chrom
pos
ref
alt
gene
transcript
classification/category
disease/condition
inheritance
frequency
pubmed_ids
comments
source_category
vip_snapshot_date
```

### Step 2: Normalize to your registry key

Use the same rules as your existing project:

```text
chrom-pos-ref-alt
GRCh37/GRCh38 status
liftOver status
REF allele match status
normalization status
```

Be careful because at least one later study using Neptune reportedly downloaded the VIP database in VCF format and lifted it from **GRCh37 to GRCh38**, so you need to check assembly before joining. 

### Step 3: Cardiac gene-panel filter

Apply your current primary/sensitivity/TTN logic:

* primary cardiogenetics genes
* sensitivity genes
* TTN separate
* quarantine/remove genes
* phenocopy/syndromic-overlap genes
* PGx/risk allele exclusion

This matters because Neptune includes Mendelian, PGx, CNV, PRS, and risk allele reporting, not only monogenic cardiac variant pathogenicity. 

### Step 4: Overlap analysis before modeling

Make a table like:

| Neptune category                 | Count |
| -------------------------------- | ----: |
| Raw VIP records                  |       |
| Coordinate-normalized records    |       |
| In your cardiogenetics panel     |       |
| In primary gene panel            |       |
| Already in ClinVar slice         |       |
| Already in eMERGE                |       |
| Already in HiRO                  |       |
| Already in CardioBoost           |       |
| New unique variants              |       |
| New unique B/LB                  |       |
| New unique VUS                   |       |
| New unique P/LP                  |       |
| Discordant with current label    |       |
| Same coordinate, same label      |       |
| Same coordinate, different label |       |

This will tell you if Neptune is actually useful.

### Step 5: Use it as one of three things

I would classify each VIP variant as:

| Neptune result                                   | How to use                                                           |
| ------------------------------------------------ | -------------------------------------------------------------------- |
| Same coordinate and same label as ClinVar/eMERGE | Concordance feature or confidence support                            |
| Same coordinate but discordant label             | Discordance audit, exclude from primary training or sensitivity only |
| New coordinate, cardiac gene, clean B/LB/P/LP    | Possible secondary training source                                   |
| New coordinate, cardiac gene, VUS                | VUS scoring/evaluation pool                                          |
| Risk allele / PGx / PRS SNP                      | Exclude from primary pathogenicity model                             |
| CNV                                              | Separate CNV analysis only, not your SNV/indel CatBoost model        |

## How this affects your current model results

Your current binary model is already very strong: corrected matrix with 85,677 variant-level rows, 42,990 binary rows, 41,731 VUS scoring rows, SpliceAI available for 77,141 rows, and source-held-out validation across HiRO, eMERGE, and CardioBoost. The source-held-out binary model had AUROC 0.9877 on HiRO, 0.9910 on eMERGE, and 0.9520 on CardioBoost. 

Because performance is already this strong, adding Neptune directly to training could create more problems than benefit unless you first prove it contributes **new, independent, non-overlapping cardiac labels**.

For the paper, Neptune is most valuable if it can support one of these claims:

1. **External validation:** “After removing all variants overlapping ClinVar, eMERGE, HiRO, and CardioBoost, the model generalized to an independent Neptune/VIP cardiac subset.”
2. **Source concordance:** “Neptune/VIP labels were used to quantify cross-source concordance and identify discordant variants.”
3. **Sparse-gene enrichment:** “Neptune/VIP contributed additional curated variants in sparse cardiac genes.”
4. **Reanalysis relevance:** “A model like ours could prioritize variants for reanalysis, similar to the Neptune/eMERGE reanalysis workflow.”

The last point is especially nice because Neptune’s manuscript already emphasizes reanalysis. They compared ClinVar snapshots, reviewed upgrade/downgrade candidates, and issued corrected reports after variant reinterpretation. 

## My recommendation

Download Neptune/VIP and run the overlap audit, but keep it **out of the primary model** until the audit is done.

My best prediction:

**It will contain cardiogenetics genes and ACMG/AMP-based final classifications. It will duplicate a lot of ClinVar/eMERGE. It may still add a useful smaller HeartCare/HGSC/internal-curated cardiac subset. Its best role is external/concordance/sparse-gene evidence, not simply increasing your training N.**

For your current study, the strongest design would be:

```text
Primary model:
ClinVar + non-Neptune curated sources, with HiRO/eMERGE/CardioBoost source-held-out validation.

Neptune analysis:
Independent audit after removing all overlapping coordinates.
Report:
1) overlap with existing registry
2) label concordance
3) new cardiac variants added
4) performance on non-overlapping Neptune cardiac variants
5) discordant-label examples for reanalysis relevance
```

That would make Neptune a strength rather than a leakage risk.

[1]: https://gitlab.com/api/v4/projects/bcm-hgsc%2Fneptune/repository/tree?path=data%2Fvip&per_page=100 "gitlab.com"
[2]: https://gitlab.com/bcm-hgsc/neptune/-/raw/master/README.md "gitlab.com"



Here is the dataset location.

The Neptune repo is on GitLab here: **BCM-HGSC / neptune**. The VIP dataset files are inside:

```text
data/vip/
```

GitLab’s API confirms two downloadable VIP snapshot files in that folder:

```text
data/vip/2020-09-03.vip.txt
data/vip/2020-10-13.vip.txt
```

([GitLab][1])

Use the newer one first:

```bash
mkdir -p datasets/neptune/raw

curl -L \
  "https://gitlab.com/api/v4/projects/bcm-hgsc%2Fneptune/repository/files/data%2Fvip%2F2020-10-13.vip.txt/raw?ref=master" \
  -o datasets/neptune/raw/2020-10-13.vip.txt
```

Also download the older snapshot for comparison/reanalysis logic:

```bash
curl -L \
  "https://gitlab.com/api/v4/projects/bcm-hgsc%2Fneptune/repository/files/data%2Fvip%2F2020-09-03.vip.txt/raw?ref=master" \
  -o datasets/neptune/raw/2020-09-03.vip.txt
```

Then check the file:

```bash
ls -lh datasets/neptune/raw/
head -n 5 datasets/neptune/raw/2020-10-13.vip.txt
wc -l datasets/neptune/raw/2020-10-13.vip.txt
```

For browsing manually, open the GitLab repo and go to:

```text
neptune → data → vip
```

The repo itself is publicly visible and was created by BCM-HGSC. ([GitLab][2])

After download, the first thing to do is **not modeling**. First run an overlap audit against your existing `chrom-pos-ref-alt` registry to see how much is truly new versus ClinVar/eMERGE duplication.

[1]: https://gitlab.com/api/v4/projects/bcm-hgsc%2Fneptune/repository/tree?path=data%2Fvip&per_page=100&ref=master "gitlab.com"
[2]: https://gitlab.com/bcm-hgsc/neptune "BCM-HGSC / neptune · GitLab"
