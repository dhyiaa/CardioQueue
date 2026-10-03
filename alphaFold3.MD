This is actually one of the coolest parts of AlphaFold 3, and the answer is it depends on what you mean by "protein."

A single job can contain multiple biological molecules, not just one protein. That's one of the biggest improvements over AlphaFold 2.

For example, one job could include:

Protein A (wild type)
Protein B (binding partner)
DNA
ATP
Mg²⁺
RNA

all in the same prediction if they're expected to interact.

However, you cannot upload 10 unrelated proteins and ask AlphaFold 3 to analyze them independently in one job. A job predicts one molecular complex.

For your project

Suppose you're studying MYH7.

One job could be:

MYH7 wild type

Another job:

MYH7 p.Arg403Gln

Then you compare the outputs.

Or, if MYH7 binds another protein:

Job 1

MYH7 WT
MYBPC3

Job 2

MYH7 Arg403Gln
MYBPC3

Then compare whether the interface changes.

Your 30 jobs/day

That means approximately:

WT + mutant = 2 jobs per variant (if comparing individually)
Around 15 variants/day

unless you become clever.

Here's something interesting...

For many proteins, the wild-type structure is already known (either experimentally or in the AlphaFold Protein Structure Database).

That means you often don't need to rerun the WT.

Instead:

Use the existing WT structure.
Run only the mutant.

Now each variant only costs 1 job instead of 2.

Your 30 jobs/day becomes roughly 30 variants/day.

Even better...

Since your project is focused on 80,000 variants, I wouldn't even think about running all of them.

I'd instead create something like:

50 pathogenic variants
50 benign variants
100 interesting VUS

(~200 variants total)

Those could become a mechanistic validation experiment in your paper, rather than a feature for every row.

One thing I want to check for you

I'm actually curious whether the AlphaFold 3 Server allows multiple sequence variants of the same protein in one prediction (for example, WT and mutant chains simultaneously) or whether each mutated sequence requires a separate job. That detail isn't obvious from the interface, and it could potentially cut your workload in half if supported. I can look into the documentation carefully and tell you exactly how to optimize your 30 jobs/day quota.

but in my case I need to upload one protein right? or there is a way to put the other protein/DNA/RNA that interacts with the proteins? for my study? uk my data strucures

This is actually where I think your project could become very interesting. I've been thinking about your data structure (85k variants, one row per variant, hundreds of engineered features, CatBoost classifier), and there are two very different ways you could use AlphaFold 3.

Option 1 (easy, probably not that useful)

For each variant:

Protein only
WT sequence

and

Protein only
Mutant sequence

Compare the structures.

This is basically an improved version of what you're already doing with AlphaFold + FoldX. You'll see local conformational changes, confidence changes, etc.

I don't think this adds much beyond your current structural features.

Option 2 (this is where AF3 shines)

Suppose your protein normally binds something.

Examples:

another protein
ATP
calmodulin
DNA
RNA
calcium
sodium
potassium

Instead of submitting only

Protein

you submit

Protein
+
binding partner

Now AF3 predicts the entire complex.

Then you mutate the protein.

Mutant protein
+
same binding partner

Now you compare

interface
binding geometry
contacts
confidence
predicted interface quality

This is something FoldX cannot really do.

Now comes the important question...

You asked:

"I don't know the interacting proteins."

That depends entirely on the gene.

For example,

RYR2

You could include

calmodulin
FKBP12.6

because they're known interactors.

KCNQ1

Could include

KCNE1
SCN5A

Could include

calmodulin
ankyrin-G
LMNA

Different interactors.

MYH7

Different interactors again.

So you do not have to discover the binding partner yourself. Many cardiac proteins already have well-characterized interaction partners in the literature.

Here's where I think your project could become genuinely novel

Right now your model has things like

REVEL
AlphaMissense
SpliceAI
gnomAD
FoldX ΔΔG
AlphaFold pLDDT

Imagine adding another feature family:

Interaction features

Predicted interface confidence change

Predicted interface distance

Interface residue disruption

Loss of hydrogen bonds

Loss of salt bridges

Predicted binding affinity change

Predicted interface area

Those features are completely different from ΔΔG.

ΔΔG asks:

Does this mutation destabilize the folded protein?

Interaction features ask:

The protein folds normally...

...but can it still bind its biological partner?

Those are different mechanisms of pathogenicity.

The challenge for your 85K dataset

This is the catch.

Every gene has different partners.

You can't realistically automate this for 120+ genes without building an enormous interaction database.

For a paper, I would not try to generate AF3 interaction features for all 85,677 variants.

Instead, I'd do it for a focused set—for example:

the 20–30 highest-impact cardiac genes,
or 100–200 representative variants,
or variants where your CatBoost and AlphaMissense predictions disagree.

That would make an excellent mechanistic validation section.

I actually have another idea that I think is even stronger than this.

Knowing your entire pipeline, I think there is a way to use AlphaFold 3 to generate high-throughput structural features that are much more compatible with an 85k-row machine learning dataset than running one AF3 prediction per variant with manually chosen binding partners. It would require building the feature extraction around known protein structures and interaction annotations rather than treating AF3 as a per-variant prediction tool. I think that approach has a much higher chance of becoming a publishable extension to your current work.

what do u mean tell me more
Thought for a couple of seconds

I mean: do not use AlphaFold 3 as “run every mutation one by one.” That is too slow and not scalable.

Instead, use AlphaFold/structure data to create variant-level features that fit your 85K CatBoost table.

The better idea

For each variant, ask:

Is this amino acid located in an important structural or interaction region?

Not:

Let me rerun AF3 for this mutation.

So for each variant you add features like:

Is residue in known protein interface? yes/no
Distance to ligand-binding site
Distance to ion-binding site
Distance to DNA/RNA-binding site
Distance to protein-protein interface
AlphaFold confidence at residue
Solvent accessibility
Buried vs exposed
In functional domain?
In transmembrane helix?
In ATP-binding pocket?
In calcium-binding region?
In myosin/actin interaction region?

This gives broad coverage because you only need the wild-type protein structure and annotations, not a new AF3 job for every mutation.

Then where does AF3 fit?

AF3 becomes useful for selected complexes, not all variants.

Example:

For KCNQ1:

KCNQ1 + KCNE1

AF3 predicts the complex. Then you map your KCNQ1 variants onto the structure and ask:

Is this variant at the KCNQ1-KCNE1 interface?
How close is it to the interface?
Is it buried in the interaction surface?

Now 1 AF3 job can help annotate many variants in that protein, because once you have the complex, you can map all KCNQ1 variants onto it.

That is the key.

Why this is better than ΔΔG

ΔΔG asks:

Does the mutation destabilize the protein?

This new AF3-style feature asks:

Is the mutation located where the protein binds something important?

A variant can have normal ΔΔG but still be pathogenic because it disrupts binding.

Practical strategy for your paper

I would do this:

Pick 5–10 high-value genes, like MYH7, MYBPC3, TNNT2, TNNI3, RYR2, KCNQ1, SCN5A, LMNA.
Identify known binding partners or ligands.
Run AF3 for the wild-type complex, not every variant.
Map all your variants onto that complex.
Add features such as “interface residue,” “distance to interface,” “distance to ligand/ion,” and “local confidence.”
Test whether adding these features improves CatBoost performance or helps explain VUS/pathogenic-like discordant cases.

So the novel angle is:

AF3 is not used as a per-variant mutation engine. It is used as a complex-mapping engine to create interaction-aware features for many variants at once.

eading my MS/SI/and draft docs, 1) give me a list of all the genes I looked at (how many and what are they?) , I am trying to see how many runs I would need, 2) make a detailed plan and description, 1) why are we doing this, 2) how are we doing this, 3) step by step, for my CODEX to understand this well, including how we will get the protein parterns/other parterns (ions and others), and how we will get their sequence ? I only I guess I have the gene as the input, , make a detailed plan fully

I’ll first extract the gene panel from your uploaded manuscript/SI/proposal files, then I’ll verify interaction-partner sources online and turn this into a Codex-ready implementation plan.

I found a key limitation in the uploaded docs: they mention a 53-gene panel, but the full 53-symbol list is not printed in the manuscript/SI text you uploaded. I can still give you the explicit gene lists that are present, and a plan for Codex to recover the full list from your local matrix.

Thought for 2m 33s

Yes, this is the right way to think about it. But I need to separate what your uploaded files prove from what Codex should recover from your actual local TSV matrix.

1) What genes did your uploaded MS/SI/draft explicitly mention?

Your files confirm the project has an analysis-ready matrix of 85,677 unique GRCh38 variant rows, with 42,990 binary training rows, and the protein-structure family currently includes AlphaFold/UniProt/DSSP/FreeSASA/FoldX features, but FoldX ΔΔG is sparse: 313 rows in the final MS/SI table, or 0.37% coverage.

Your proposal/status file also mentions a 53-gene panel, but the full 53-gene symbol list is not printed in the uploaded manuscript/SI text I could see. It says the gnomAD gene constraint table is available locally for the 53-gene panel, so Codex should recover the actual full list directly from your local matrix rather than from the pasted prose.

Formal gene-QC table in the SI: 27 explicitly listed genes
Bucket	Genes	Count	Meaning
Excluded or quarantined	AKAP9, DMD, FPGT, FPGT-TNNI3K, KCNE1B, KNCH2	6	Removed from primary panel
Sensitivity/separate handling	ANK2, SOS1, TNNI3K, TRPM4, TTN	5	Held out of primary ClinVar training or handled separately
Conditional phenocopy / syndromic overlap	EMD, GLA, HRAS, KRAS, LAMP2, MAP2K2, NRAS, PTPN11, RAF1, RIT1, SLC22A5	11	Retained with conditional-overlap flags
Retain, monitor	CACNB2, FHOD3, MIB1, SCN1B, SNTA1	5	Retained but monitored

That formal QC table totals 27 explicitly named genes/symbols.

One important QC warning: “KNCH2” appears in the uploaded text, but this may be a typo for KCNH2. Codex should not silently treat that as correct. It should check whether KNCH2 exists in the matrix and whether KCNH2 is present, because accidentally excluding or misspelling KCNH2 would be a major cardiogenetics error.

Additional genes mentioned elsewhere in the draft/proposal

Your proposal also names sparse-gene examples: PLN, ACTC1, MYL2, MYL3, CALM1, CALM2, CALM3, TRDN, SMAD3, TGFBR1, TGFBR2.

Your AF3 plan/draft also names candidate complex genes: KCNQ1, KCNE1, SCN5A, RYR2, MYBPC3, MYH7, TNNI3, TNNT2, TNNC1, TPM1, ACTC1, LMNA, EMD.

Across the uploaded prose, the unique gene/protein symbols explicitly mentioned for the project/AF3 discussion are therefore 49 unique symbols:

ACTC1, AKAP9, ANK2, CACNB2, CALM1, CALM2, CALM3, DMD, EMD,
FHOD3, FPGT, FPGT-TNNI3K, GLA, HRAS, KCNE1, KCNE1B, KCNQ1,
KNCH2, KRAS, LAMP2, LMNA, MAP2K2, MIB1, MYBPC3, MYH7, MYL2,
MYL3, NRAS, PLN, PTPN11, RAF1, RIT1, RYR2, SCN1B, SCN5A,
SLC22A5, SMAD3, SNTA1, SOS1, TGFBR1, TGFBR2, TNNC1, TNNI3,
TNNI3K, TNNT2, TPM1, TRDN, TRPM4, TTN

But again: this is not the same as the final 53-gene panel. It is the set of genes explicitly printed in your uploaded text. Codex should recover the exact 53 from your local matrix.

2) How many AlphaFold 3 runs would this need?

The answer depends on strategy.

Bad/scalable-unfriendly strategy

If you run:

WT protein
Mutant protein

for every variant, this becomes impossible.

You have 85,677 variant rows. Even if only missense variants are used, your VEP table reports ~40,280 missense variants in an intermediate matrix. That is not practical by AlphaFold Server quota or manuscript defensibility.

Better strategy

Use AF3 to predict wild-type complexes, then map many variants onto each complex.

That means:

1 AF3 job = 1 biological complex

not

1 AF3 job = 1 variant

Example:

KCNQ1 + KCNE1 + CALM + Ca2+

One job can then annotate many KCNQ1 and KCNE1 variants by asking whether each variant lies at the interface or ligand/ion-binding region.

AlphaFold Server can model biomolecular complexes containing combinations of proteins, DNA, RNA, ligands, ions, and chemical modifications. It also uses confidence outputs such as pLDDT and PAE, and complex-level outputs include metrics such as pTM/ipTM/pairwise interface metrics.

Practical run estimates
Plan	Jobs	Time at 30 jobs/day	My recommendation
3 pilot complexes	3	Same day	Best starting point
10 high-value complexes	10	Same day	Strong pilot substudy
30 selected complexes	30	1 day	Good maximum first batch
53 genes, protein-only	53	~2 days	Not worth it because AF2 already covers this better
53 genes × up to 5 partners	Up to 265	~9 days	Too broad for first pass; do only after automated prioritization
20 mutation-specific WT/mutant tests	~21 to 40	1 to 2 days	Good secondary validation only

I would not use multiple Google accounts to multiply quota. It looks like quota circumvention and would be hard to defend in Methods. Also, AlphaFold Server output terms restrict using outputs to train biomolecular structure-prediction models similar to AlphaFold; your CatBoost pathogenicity model is different, but you should still frame this as non-commercial, hypothesis-generating structural annotation, not as scaling AF3 into a derived structure model.

3) Why are we doing this?

Your current ΔΔG/FoldX layer asks:

Does this amino-acid substitution destabilize the folded protein?

That is useful, but it misses another mechanism:

The protein may fold okay, but the variant may disrupt binding to a partner,
ion, ligand, DNA/RNA, or complex interface.

That is where AF3 is useful. AF3 is not mainly replacing FoldX. It is adding a new feature family:

interaction-aware structural context

Your current manuscript already says protein-structure features are valuable but sparse, and FoldX ΔΔG should be interpreted as informative for mapped missense variants rather than a universal predictor.

So the scientific claim becomes:

FoldX ΔΔG captures predicted stability disruption. AlphaFold 3 complex mapping captures whether variants fall at biologically important molecular interfaces. These mechanisms are complementary.

That is a much stronger and cleaner story than “we used AF3 because it is new.”

4) Does this require only the protein, or known interacting partners?

For your goal, it requires known interacting partners or biologically justified ligands/ions.

Protein-only AF3 input:

KCNQ1

This mostly gives you another structure prediction. That overlaps with AlphaFold 2, pLDDT, DSSP, FreeSASA, UniProt domains, and FoldX.

Interaction-aware AF3 input:

KCNQ1 + KCNE1 + CALM + Ca2+

This gives you a complex. Then you can map variants to:

KCNQ1-KCNE1 interface
KCNQ1-CALM interface
distance to Ca2+
local confidence
interface PAE

That is the useful feature set.

5) Priority AF3 targets and up to 5 known partners/entities

This should be automated for the full panel, but these are the best first targets.

Target	Candidate partners/entities	Why this is worth modeling
KCNQ1	KCNE1, KCNE2/KCNE3 if desired, CALM, AKAP9, Ca²⁺	KCNQ1 associates with KCNE1 to form the cardiac IKs current, and UniProt reports KCNQ1 interactions including AKAP9; recent structural work also supports KCNQ1-KCNE1-CaM complex modeling.
SCN5A / Nav1.5	CALM1/2/3, ANK3 domain, PKP2, FGF13, SNTA1	Nav1.5 is regulated by intracellular partners including calmodulin and ankyrin-related binding, and UniProt/search results list interactions including ANK3 and PKP2.
RYR2	FKBP1B, FKBP1A, CALM, CASQ2, TRDN, Ca²⁺	UniProt reports RYR2 interaction with FKBP1A/FKBP1B, and NCBI describes the cardiac channel complex with FKBP1B; CASQ2/TRDN are relevant SR calcium-release complex partners.
MYBPC3 / MYH7	MYH7/myosin fragment, ACTC1/F-actin, thin-filament fragment, TTN region if supported, ATP/Mg²⁺ for myosin-focused jobs	MYBPC3 binds myosin heavy chain, F-actin, and native thin filaments, so interface mapping could capture sarcomere biology missed by ΔΔG.
Troponin/thin filament genes: TNNI3, TNNT2, TNNC1, TPM1, ACTC1	TNNI3, TNNT2, TNNC1, TPM1, ACTC1, Ca²⁺	Troponin I/T/C form the thin-filament regulatory complex; TNNT2 binds tropomyosin and TNNC1 binds calcium, making this an ideal AF3 multi-gene complex.
LMNA	EMD, SUN1, SUN2, TMEM43, MLIP/TMEM201	UniProt lists LMNA interactions including SUN2/TMEM43 and other nuclear-envelope partners; EMD/SUN proteins are also plausible complex-mapping partners.
6) Codex-ready implementation plan
Goal

Build an AlphaFold 3 interaction-aware structural annotation layer for the CatBoost cardiogenetic variant classifier.

The new feature family should not replace FoldX ΔΔG. It should add interface/partner/ligand context:

FoldX ΔΔG = predicted stability effect
AF3 complex mapping = predicted interface/location context

The output should be a feature table that can merge into the existing variant-level matrix by:

variant_id
gene
UniProt accession
protein position
reference amino acid
alternate amino acid
Phase 0 — Recover the exact gene panel from the local data

Codex should not rely on pasted prose for the final gene list.

Use the current best local matrix:

datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv

and, if available, the final analysis-ready matrix:

datasets/modeling/processed/...analysis_ready...

Codex should identify the gene column using this priority order:

gene
symbol
gene_symbol
vep_symbol
hgnc_symbol
clinvar_gene

Then generate:

results/af3/gene_panel_audit.tsv

with columns:

gene
n_total_variants
n_missense
n_plp
n_blb
n_vus
n_hgvsp_nonmissing
n_protein_position_nonmissing
n_uniprot_mapped
panel_bucket
ttn_separate_flag
excluded_or_quarantine_flag

Expected from your docs:

analysis-ready matrix: 85,677 rows
formal QC table genes explicitly printed: 27
local panel mentioned: 53 genes
FoldX final MS coverage: 313 rows / 0.37%

The code should also check:

Is KNCH2 present?
Is KCNH2 present?
Are both present?

and write a warning if KNCH2 appears but KCNH2 does not.

Phase 1 — Select AF3 target genes

Do not run every gene blindly.

Create a priority score:

priority_score =
  3 * n_plp_missense
+ 2 * n_vus_missense
+ 2 * n_catboost_alphaMissense_discordant
+ 2 * n_foldx_missing_but_uniprot_mapped
+ 1 * known_complex_available
- 5 * excluded_or_quarantine_flag

Start with three groups:

Group A: first pilot
KCNQ1-KCNE1-CALM-Ca2+
troponin/thin filament complex
SCN5A-CALM
Group B: second pilot
RYR2-FKBP1B or RYR2 domain-FKBP1B
MYBPC3-MYH7/actin fragment
LMNA-EMD/SUN2/TMEM43
Group C: automated expansion

Top 10 to 20 genes from the matrix after filtering for:

n_missense >= 50
n_uniprot_mapped >= 30
known_partner_count >= 1
not excluded/quarantined
Phase 2 — Retrieve canonical protein sequences

For each target gene and each partner gene, Codex should retrieve the canonical reviewed human UniProt entry.

Preferred query:

UniProt REST API
query = gene_exact:GENE AND organism_id:9606 AND reviewed:true
format = json

UniProt provides programmatic APIs for querying and retrieving protein data, including UniProtKB entries and sequences.

Save:

datasets/feature_sources/alphafold3/interim/uniprot_gene_to_accession.tsv
datasets/feature_sources/alphafold3/interim/uniprot_sequences.fasta
datasets/feature_sources/alphafold3/interim/uniprot_features.jsonl

Columns:

gene
uniprot_accession
reviewed
canonical
protein_name
sequence
length
isoform_note
source_url
retrieval_date

QC rules:

Use reviewed Swiss-Prot entries first.
Prefer canonical isoform.
Do not use unreviewed TrEMBL unless no reviewed human entry exists.
If multiple reviewed entries exist, require manual review.
Phase 3 — Retrieve candidate protein partners

Codex should collect partners from four sources, then rank them.

Source priority
Manual curated whitelist for obvious cardiac complexes.
Reactome for pathway/complex participants.
UniProt/IntAct for curated physical interactions.
STRING only for high-confidence backup, because STRING includes functional associations, not only direct binding.

Reactome has a content service API for programmatic access to curated pathway and complex data, and STRING has an API for accessing interaction networks. UniProt binary interactions are derived from IntAct, and IntAct curates experimentally derived interaction data from the literature.

Create:

datasets/feature_sources/alphafold3/interim/candidate_partners.tsv

Columns:

target_gene
target_uniprot
partner_gene
partner_uniprot
partner_type
source_database
evidence_type
evidence_score
pubmed_ids
same_pathway_flag
known_physical_interaction_flag
cardiac_relevance_flag
protein_length
token_feasible_full_length
recommended_sequence_mode
manual_review_status

Partner type values:

protein
protein_domain
peptide
ion
ligand
DNA
RNA
modification

Ranking rule:

physical curated interaction > Reactome complex participant > high-confidence STRING association
human evidence > non-human/by-similarity evidence
cardiac/sarcomere/channel relevance > generic cellular interaction
short/domain-feasible partner > huge full-length partner
Phase 4 — Retrieve ions, ligands, and non-protein entities

For ions and ligands, Codex should not guess.

Use:

UniProt binding-site/cofactor annotations
Reactome reaction participants
RCSB/PDB co-crystal ligands
ChEBI / CCD identifiers

Examples:

Ca2+ for TNNC1/troponin and RYR2-related jobs
Mg2+/ATP for myosin/ATPase-related jobs if modeling MYH7 motor-domain biology
K+ or Na+ only if the AF3 input design is explicitly channel/pore focused

For AlphaFold Server, ions/ligands should be entered using supported chemical component identifiers where available, or manually through the server interface if supported. AF3 supports ligands and ions as input molecular entities.

Create:

datasets/feature_sources/alphafold3/interim/candidate_nonprotein_entities.tsv

Columns:

target_gene
complex_name
entity_name
entity_type
ccd_or_chebi_id
evidence_source
evidence_note
include_in_first_pass
Phase 5 — Design AF3 jobs

Each AF3 job should represent one biological hypothesis:

target protein/domain + partner protein/domain + optional ion/ligand

Do not submit 10 unrelated proteins in one job. Submit molecules that are expected to form one complex.

Create:

datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv

Columns:

job_id
complex_name
target_gene
target_uniprot
target_sequence_mode
target_range_start
target_range_end
partner_genes
partner_uniprots
partner_ranges
ions
ligands
total_tokens_estimated
reason_for_job
evidence_sources
priority
manual_ready_for_submission
af3_status
output_path

Sequence mode values:

full_length
domain
peptide
fragment

Rules:

If total tokens > 5000, trim to known interacting domains.
For RYR2, do not submit full RYR2 + partners. Use domain-level jobs.
For MYBPC3/MYH7/thin filament jobs, avoid huge full filament assemblies. Use domains/fragments.
For troponin, one multi-chain complex can annotate several genes at once.
Phase 6 — Run AF3 and store outputs

For each completed server job, download:

model.cif
summary JSON / confidence JSON if available
PAE file / confidence metrics file
input sequences
server metadata
screenshots or exported job details if needed for audit

Store:

datasets/feature_sources/alphafold3/raw/{job_id}/

with:

input.fasta
input_entities.json
model_1.cif
model_2.cif
model_3.cif
model_4.cif
model_5.cif
confidence.json
pae.json
server_metadata.json
README.md

Also create:

datasets/feature_sources/alphafold3/interim/af3_completed_jobs.tsv
Phase 7 — Parse AF3 structures and extract features

Use Python libraries:

gemmi
biopython
pandas
numpy
scipy
freesasa

For each complex, identify:

target chain
partner chain(s)
ligand/ion coordinates
residue numbering map
pLDDT/B-factor fields
PAE matrix if available

Extract residue-level features:

af3_complex_available
af3_complex_name
af3_target_chain
af3_residue_mapped
af3_local_plddt
af3_min_distance_to_any_partner
af3_min_distance_to_each_partner
af3_interface_residue_4A
af3_interface_residue_5A
af3_interface_residue_8A
af3_partner_contact_count_5A
af3_min_distance_to_ligand
af3_min_distance_to_ion
af3_near_ligand_5A
af3_near_ion_5A
af3_mean_partner_pae_local
af3_max_partner_confidence
af3_interface_surface_area_if_available

Definitions:

interface_residue_5A = any heavy atom in target residue within 5 Å of any heavy atom in partner chain

near_ligand_5A = any target residue heavy atom within 5 Å of ligand atom

local_plddt = mean pLDDT over target residue atoms

local_partner_pae = mean PAE between target residue and partner-chain residues within contact window

Save:

datasets/feature_sources/alphafold3/processed/af3_residue_features.tsv
Phase 8 — Map AF3 residue features to your variant table

Merge by:

gene
uniprot_accession
protein_position
reference_amino_acid

Before merging, require reference amino-acid QC:

variant_ref_aa == uniprot_sequence[protein_position]

If mismatch:

af3_mapping_status = ref_aa_mismatch

Do not impute.

Create:

datasets/feature_sources/alphafold3/processed/af3_variant_features.tsv

Columns:

variant_id
gene
hgvsp
protein_position
ref_aa
alt_aa
af3_mapping_status
af3_complex_available
af3_any_interface_residue_5A
af3_min_distance_to_partner
af3_min_distance_to_ligand
af3_min_distance_to_ion
af3_local_plddt
af3_mean_interface_pae
af3_best_complex_name
af3_n_complexes_covering_residue
af3_feature_missing_reason

Missing reasons:

no_uniprot_mapping
no_protein_position
not_missense
gene_not_selected_for_af3
residue_not_in_modeled_domain
ref_aa_mismatch
af3_job_failed
low_confidence_structure
Phase 9 — Test whether AF3 helps the CatBoost model

Do not claim broad performance improvement unless tested properly.

Run three comparisons:

A. Full-table model
Baseline features
vs
Baseline + AF3 features

Expected issue: most rows will have missing AF3 features, so global AUROC may barely move.

B. AF3-covered subset

Only variants where:

af3_complex_available == true
af3_mapping_status == ok

This tests whether AF3 improves performance where it actually applies.

C. Mechanistic case study

Most publishable:

P/LP variants with low FoldX ΔΔG but high AF3 interface signal
VUS predicted pathogenic-like and near AF3 interface
CatBoost/AlphaMissense discordant variants explained by AF3 interface position
Benign variants far from interfaces despite structural mapping

This is likely stronger than chasing a tiny AUROC improvement.

10) Best first 30 AF3 runs

Do not spend the first 30 on random genes.

I would use:

1. KCNQ1 + KCNE1 + CALM + Ca2+
2. KCNQ1 + KCNE1 only
3. KCNQ1 cytoplasmic C-terminal domain + CALM
4. SCN5A IQ/C-terminal domain + CALM
5. SCN5A ankyrin-binding region + ANK3 domain/peptide
6. SCN5A + PKP2 domain/fragment if feasible
7. TNNI3 + TNNT2 + TNNC1 + Ca2+
8. TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+
9. TNNT2 + TPM1 + ACTC1
10. TNNC1 + Ca2+
11. MYBPC3 C0-C2 region + MYH7 fragment
12. MYBPC3 C-terminal region + myosin fragment
13. MYBPC3 domain + ACTC1/F-actin fragment
14. MYH7 motor domain + ATP/Mg2+ if AF3 ligand setup is feasible
15. RYR2 selected domain + FKBP1B
16. RYR2 selected domain + FKBP1A
17. RYR2 selected domain + CALM
18. CASQ2 + TRDN
19. RYR2 luminal/domain fragment + CASQ2/TRDN region if justified
20. LMNA + EMD
21. LMNA + SUN1 domain
22. LMNA + SUN2 domain
23. LMNA + TMEM43 domain
24. PLN + ATP2A2/SERCA2a domain if ATP2A2 is in your expanded biological target set
25. CALM1 + RYR2 peptide/domain
26. CALM2 + SCN5A peptide/domain
27. CALM3 + KCNQ1 peptide/domain
28. ACTC1 + TPM1
29. ACTC1 + TNNT2/TNNI3/TNNC1
30. Repeat the highest-priority complex with alternate seed if first output confidence is weak

Some of these require manual domain trimming. That is fine. The manifest should document exact residue ranges.

Bottom line

Yes, this is the right direction, but not as “run AF3 on all variants.”

The clean plan is:

Use AF3 to model selected wild-type cardiac protein complexes.
Map many variants onto those complexes.
Extract interface/ligand/ion-distance features.
Merge those features into the CatBoost matrix.
Evaluate whether interface-aware features explain pathogenic/VUS behavior beyond ΔΔG.

This fits your current manuscript perfectly because your existing structure layer is sparse and stability-focused, while AF3 adds a distinct interaction-aware mechanism.

I want this (Better strategy

Use AF3 to predict wild-type complexes, then map many variants onto each complex.

That means:

1 AF3 job = 1 biological complex

not

1 AF3 job = 1 variant

Example:

KCNQ1 + KCNE1 + CALM + Ca2+

One job can then annotate many KCNQ1 and KCNE1 variants by asking whether each variant lies at the interface or ligand/ion-binding region) , my question is how can I generate this (nteraction-aware AF3 input:

KCNQ1 + KCNE1 + CALM + Ca2+) and also would the model give me this? (This gives you a complex. Then you can map variants to:

KCNQ1-KCNE1 interface
KCNQ1-CALM interface
distance to Ca2+
local confidence
interface PAE) or the coordinate of this, and then I will need to build a engine to do the calcuations locally right? which is fine, but the alphafold will output some stuff about this (KCNQ1-KCNE1 interface
KCNQ1-CALM interface
distance to Ca2+
local confidence
interface PAE) 
Thought for 33s

Yes. You understand it correctly.

For your study, AlphaFold 3 gives you the predicted complex structure and confidence outputs, but you still need to build the local feature-extraction engine that converts those outputs into CatBoost-ready columns.

Your current structure layer is sparse: your matrix has 85,677 rows, but FoldX ΔΔG is available for only 313 / 85,677 rows, 0.37%, so this AF3-complex approach is attractive because one wild-type complex can annotate many variants by position rather than requiring one run per mutation.

1) How do you generate KCNQ1 + KCNE1 + CALM + Ca2+?

You need to give AlphaFold Server the molecular entities in the complex.

For this example, the entities would be:

Protein 1: KCNQ1 sequence
Protein 2: KCNE1 sequence
Protein 3: CALM1 sequence, or CALM2/CALM3 depending on your biological choice
Ion: Ca2+

In the AlphaFold Server UI, you can do this manually by adding multiple protein entities. The server documentation says that for proteins you can paste a single-letter amino-acid sequence or FASTA, and if you have multiple proteins you can put all sequences into one FASTA file. It also says ligands, ions, and modifications are selected from the interface using CCD codes.

So the practical version is:

>KCNQ1_HUMAN
[canonical KCNQ1 amino-acid sequence]

>KCNE1_HUMAN
[canonical KCNE1 amino-acid sequence]

>CALM1_HUMAN
[canonical calmodulin-1 amino-acid sequence]

Then in the AlphaFold Server interface, add Ca2+ as an ion/ligand entity.

If you were using local AlphaFold 3 JSON rather than the web UI, ions are represented as ligands. The AlphaFold 3 input documentation says ions are treated as ligands, for example magnesium is entered as a ligand with ccdCodes: ["MG"]; for calcium, the CCD code is generally CA.

A conceptual local JSON-style version would look like this:

{
  "name": "KCNQ1_KCNE1_CALM1_CA",
  "modelSeeds": [1],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "KCNQ1_SEQUENCE_HERE",
        "description": "KCNQ1"
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "KCNE1_SEQUENCE_HERE",
        "description": "KCNE1"
      }
    },
    {
      "protein": {
        "id": "C",
        "sequence": "CALM1_SEQUENCE_HERE",
        "description": "CALM1"
      }
    },
    {
      "ligand": {
        "id": "D",
        "ccdCodes": ["CA"],
        "description": "calcium ion"
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 4
}

For the web server, you may not need to write this JSON yourself. You can use the UI, paste the FASTA sequences, then add the ion through the entity bar. The JSON matters later because the downloaded job request can preserve exactly what you submitted.

2) Does AF3 directly give you this?

Partly yes, partly no.

AF3 gives you directly

AlphaFold Server gives you:

Predicted 3D coordinates of the complex
mmCIF structure files
pLDDT confidence
PAE data
pTM
ipTM
chain-level and chain-pair confidence summaries
job input JSON

The server documentation says each job gives five predicted structures as .cif files, summary confidence JSON files, full-data JSON files containing detailed confidence metrics such as full PAE data, the job request JSON, and terms of use.

The results page also displays backbone visualization colored by pLDDT, a PAE plot, pTM, and ipTM for complexes.

The AlphaFold 3 output documentation confirms that pLDDT is a per-atom confidence score, PAE estimates positional/orientation error between tokens, and ipTM estimates confidence in the relative placement of subunits within a complex.

AF3 does not directly give you a finished table like this

It will not automatically output:

KCNQ1-KCNE1 interface residue: yes/no
KCNQ1-CALM interface residue: yes/no
distance to Ca2+
variant is 4.2 Å from KCNE1
variant is interface-disrupting

You calculate those locally from the .cif and .json files.

So the correct mental model is:

AlphaFold 3 output = predicted coordinates + confidence data

Your local engine = converts coordinates/confidence into variant-level features
3) What exactly can your local engine calculate?

From the AF3 outputs, your local code can calculate:

A. KCNQ1-KCNE1 interface

Use the .cif coordinates.

For each KCNQ1 residue:

minimum distance from any KCNQ1 residue atom
to any KCNE1 atom

Then define:

KCNQ1-KCNE1 interface residue = yes
if min distance <= 5 Å

You can also use 4 Å or 8 Å as sensitivity thresholds.

Output columns:

af3_kcnq1_kcne1_interface_5A
af3_kcnq1_kcne1_min_dist_A
af3_kcnq1_kcne1_contact_count_5A
B. KCNQ1-CALM interface

Same idea, but partner chain is CALM.

minimum distance from KCNQ1 residue
to CALM chain

Output columns:

af3_kcnq1_calm_interface_5A
af3_kcnq1_calm_min_dist_A
af3_kcnq1_calm_contact_count_5A
C. Distance to Ca2+

The .cif file contains atom coordinates for the ion if AF3 placed it.

For each KCNQ1 residue:

minimum distance from residue atoms
to Ca2+ atom

Output columns:

af3_min_dist_to_ca_A
af3_near_ca_5A
af3_near_ca_8A
D. Local confidence

pLDDT is in the AF3 output. For your variant residue, compute:

mean pLDDT of all atoms in that residue

Output columns:

af3_local_plddt
af3_local_plddt_window_3aa
af3_local_plddt_window_5aa
E. Interface PAE

PAE is in the JSON full-data file. For a KCNQ1 residue and partner chain, compute:

mean PAE between the KCNQ1 residue token
and the partner residues it contacts

Output columns:

af3_kcnq1_kcne1_local_pae
af3_kcnq1_calm_local_pae
af3_best_partner_pae

This is important because a residue can look close to a partner in the structure, but if the inter-chain PAE is poor, the predicted interface may not be reliable.

4) So yes, you need a local calculation engine

The local engine should do this:

Input:
1. Your variant table
2. UniProt residue mapping
3. AF3 .cif structure
4. AF3 confidence JSON / PAE JSON
5. Chain mapping file: A = KCNQ1, B = KCNE1, C = CALM1, D = Ca2+

Output:
Variant-level AF3 features merged into CatBoost table

For each variant:

variant gene = KCNQ1
protein position = 325
ref AA = R
alt AA = Q

Your engine finds residue 325 in chain A and calculates:

distance to KCNE1
distance to CALM1
distance to Ca2+
local pLDDT
local interface PAE
is interface residue yes/no

Then you merge those columns into your CatBoost matrix.

5) What AlphaFold itself gives vs what you calculate
Feature	Comes directly from AF3?	You calculate locally?
3D coordinates of KCNQ1/KCNE1/CALM/Ca2+	Yes	No
Chain IDs	Yes	You map them to gene names
pLDDT	Yes	You aggregate per residue/window
PAE matrix	Yes	You extract relevant residue-partner values
pTM/ipTM	Yes	You store as complex-level QC
KCNQ1-KCNE1 interface yes/no	No	Yes
KCNQ1-CALM interface yes/no	No	Yes
Distance to Ca2+	No	Yes
Contact count	No	Yes
Variant-level CatBoost feature table	No	Yes
6) The exact pipeline for one complex

For KCNQ1 + KCNE1 + CALM1 + Ca2+:

Step 1: Get sequences

From UniProt:

KCNQ1 canonical protein sequence
KCNE1 canonical protein sequence
CALM1 canonical protein sequence
Step 2: Submit AF3 job

In AlphaFold Server:

Protein entity: KCNQ1
Protein entity: KCNE1
Protein entity: CALM1
Ion entity: Ca2+

Job name:

KCNQ1_KCNE1_CALM1_CA_WT
Step 3: Download outputs

Save:

model_0.cif
summary_confidences_0.json
full_data_0.json
job_request.json

The top-ranked model is usually rank 0, but keep all five models for QC.

Step 4: Create chain map

Something like:

A = KCNQ1
B = KCNE1
C = CALM1
D = CA

Do not assume this blindly. Check the output .cif and job_request.json.

Step 5: Parse structure locally

Use Python:

gemmi or BioPython for mmCIF parsing
numpy/scipy for distance calculations
pandas for feature table
Step 6: Create residue-level table

Example:

complex_name | target_gene | chain | residue_number | aa | local_plddt | min_dist_KCNE1 | min_dist_CALM1 | min_dist_CA | interface_KCNE1_5A | interface_CALM1_5A
Step 7: Merge to variants

Merge by:

gene + protein_position + ref_aa

Then variant-level output:

variant_id | gene | hgvsp | protein_position | af3_complex | af3_interface_KCNE1_5A | af3_min_dist_KCNE1 | af3_interface_CALM1_5A | af3_min_dist_CALM1 | af3_min_dist_CA | af3_local_plddt | af3_local_pae
7) Very important caution

Do not interpret AF3 as saying:

This variant definitely disrupts binding.

At this stage, it says:

This variant lies at or near a predicted interface / ligand / ion-binding region,
and the predicted complex has certain confidence metrics.

That is still valuable. For your model, that becomes an interaction-aware annotation, not an ACMG-level functional assay.

8) Best wording for your Methods

You could describe it like this:

For selected cardiogenetic proteins with well-established molecular partners, wild-type AlphaFold 3 complex predictions were generated using canonical UniProt protein sequences and relevant non-protein entities, such as calcium ions where biologically justified. Predicted complexes were downloaded as mmCIF and JSON confidence files. Variants were mapped to protein residues using the existing UniProt/protein-position mapping pipeline. For each mapped residue, we calculated distance to partner chains, ligand or ion proximity, local pLDDT, and interface-level PAE-derived confidence metrics. These features were treated as exploratory interaction-aware structural annotations and merged into the CatBoost feature matrix.

That is the cleanest version. It makes clear that AF3 gives the structure and confidence outputs, while your code creates the variant-level ML features.


ok final planning below: 

Yes, I agree with you. **Delete the prioritization-score idea.** For your goal, the better plan is:

> Build a **full-panel AF3 interaction map**: every gene in your cardiogenetic panel gets at least one biologically justified AF3 complex attempt, then variants are mapped onto those complex structures.

That is different from “prioritize 5–10 genes.” Your current project already has a large, mature matrix, with **85,677 analysis-ready variant rows**, **42,990 binary training rows**, and very sparse FoldX ΔΔG coverage, **313 / 85,677 rows, 0.37%**. That makes a broad AF3 interface-feature layer reasonable as a serious extension, not just a cute pilot. 

## The key correction

Do **not** do this:

```text
Pick only top 10–20 genes
```

Do this instead:

```text
Recover the full ~52–53 gene panel
For every gene, find at least one best-supported partner/complex/ligand context
Run the non-redundant set of AF3 complex jobs
Map all eligible variants onto those predicted complexes
```

The word **non-redundant** matters. For example, one troponin/thin-filament AF3 job can cover **TNNI3, TNNT2, TNNC1, TPM1, and ACTC1** together. So you may have ~52 genes, but not necessarily 52 separate jobs. Some genes share the same biological complex.

AlphaFold Server can model structures made of multiple proteins, nucleic acids, ligands, ions, and modifications, and it tries to model interactions between all molecules supplied in the job input. Each job has a **5,000-token limit**, and each user has **30 jobs/day**. ([EMBL-EBI][1])

So your timeline is realistic:

| Phase                      | Goal                                                                         | Output                                      |
| -------------------------- | ---------------------------------------------------------------------------- | ------------------------------------------- |
| 2 days figuring things out | Recover full gene list, find partners, get sequences, build AF3 job manifest | `af3_full_panel_job_manifest.tsv`           |
| 2 days production          | Submit ~31–60 curated AF3 jobs, download outputs                             | `/datasets/feature_sources/alphafold3/raw/` |
| After production           | Parse structures locally and create variant-level features                   | `af3_variant_features.tsv`                  |

I would not frame this as “prioritization.” I would frame it as **full-panel complex coverage with feasibility/QC flags**.

---

# Revised plan for Codex

## Goal

Build a full-panel AlphaFold 3 complex-mapping layer for the cardiogenetic CatBoost variant classifier.

This layer should cover the entire cardiogenetic gene panel, not only selected top-priority genes. Each gene should receive one of these outcomes:

```text
1. Covered by a direct AF3 protein-protein complex
2. Covered by a shared multi-gene complex
3. Covered by a protein-ligand/ion complex
4. Covered by a domain-level complex because full-length input exceeds 5,000 tokens
5. No defensible curated complex found; documented as no_curated_partner
```

The aim is not to run AF3 per variant. The aim is to run AF3 per **wild-type biological complex**, then map many variants onto that complex.

---

# Phase 0 — Recover the exact full panel

Codex should use the local matrix, not pasted prose.

Use:

```text
datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv
```

and, if available, the final processed analysis-ready table.

The uploaded manuscript confirms the current analysis-ready table is **85,677 rows and 420 columns**, while the intermediate table was **86,889 rows**. Manuscript-facing results should use the 85,677-row table. 

Codex should identify the gene column using this priority order:

```text
gene
symbol
gene_symbol
vep_symbol
hgnc_symbol
clinvar_gene
```

Generate:

```text
results/af3/gene_panel_audit.tsv
```

Columns:

```text
gene
n_total_variants
n_missense
n_plp
n_blb
n_vus
n_hgvsp_nonmissing
n_protein_position_nonmissing
n_uniprot_mapped
panel_bucket
ttn_separate_flag
excluded_or_quarantine_flag
af3_coverage_plan
af3_complex_group
af3_manual_review_needed
```

Also check:

```text
Is KNCH2 present?
Is KCNH2 present?
Are both present?
```

This is important because your SI/prose says **KNCH2**, but that may be a typo for **KCNH2**. The uploaded SI lists **AKAP9, DMD, FPGT, FPGT-TNNI3K, KCNE1B, and KNCH2** as excluded, and lists **ANK2, SOS1, TNNI3K, TRPM4, and TTN** as sensitivity/separate-handling genes. 

---

# Phase 1 — Remove prioritization and replace with full-panel coverage

Do **not** use this anymore:

```text
priority_score =
  3 * n_plp_missense
+ 2 * n_vus_missense
...
```

Replace it with:

```text
coverage_status =
  full_length_complex
  domain_level_complex
  ligand_or_ion_context
  shared_complex
  no_curated_partner
  manual_review_needed
```

Codex should not ask:

```text
Which genes are most important?
```

Codex should ask:

```text
For every gene, what is the best defensible AF3 biological context?
```

That means every gene gets reviewed.

---

# Phase 2 — Build a full-panel interaction manifest

Create:

```text
results/af3/af3_full_panel_partner_manifest.tsv
```

Columns:

```text
target_gene
target_uniprot
target_protein_name
target_length
target_panel_status
complex_group_id
complex_name
partner_gene_or_entity
partner_type
partner_uniprot_or_identifier
partner_length_or_token_count
interaction_source
interaction_evidence_type
pubmed_ids
reactome_pathway_or_complex
uniprot_interaction_note
intact_or_string_id
ligand_or_ion_id
recommended_af3_input_mode
domain_start
domain_end
reason_for_domain_trimming
token_estimate
manual_review_status
final_include
exclusion_reason
```

Partner types:

```text
protein
protein_domain
peptide
ion
ligand
DNA
RNA
PTM
none_found
```

For each gene, Codex should search in this order:

```text
1. Existing manual whitelist from cardiogenetics knowledge
2. UniProt reviewed human entry
3. Reactome complex/pathway partners
4. IntAct / curated physical interactions
5. PDB/RCSB co-complex evidence
6. STRING only as weak backup, not enough alone unless manually accepted
```

This is not prioritization. It is **full-panel partner assignment**.

---

# Phase 3 — Avoid redundant jobs by using complex groups

This is the key efficiency trick.

Instead of doing:

```text
TNNI3 + partner
TNNT2 + partner
TNNC1 + partner
TPM1 + partner
ACTC1 + partner
```

do:

```text
TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+
```

One job can cover multiple genes.

So Codex should build:

```text
results/af3/af3_complex_groups.tsv
```

Columns:

```text
complex_group_id
complex_name
genes_covered
protein_entities
ligand_entities
ion_entities
dna_entities
rna_entities
total_tokens_estimated
full_length_feasible
domain_trim_required
manual_submission_name
expected_variant_coverage
notes
```

Example rows:

```text
complex_group_id: AF3_COMPLEX_KCNQ1_KCNE1_CALM_CA
complex_name: KCNQ1-KCNE1-calmodulin-calcium complex
genes_covered: KCNQ1;KCNE1;CALM1;CALM2;CALM3
entities: KCNQ1, KCNE1, CALM1, Ca2+
```

```text
complex_group_id: AF3_COMPLEX_TROPONIN_THIN_FILAMENT_CA
complex_name: Cardiac troponin/thin-filament calcium complex
genes_covered: TNNI3;TNNT2;TNNC1;TPM1;ACTC1
entities: TNNI3, TNNT2, TNNC1, TPM1, ACTC1, Ca2+
```

```text
complex_group_id: AF3_COMPLEX_LMNA_EMD
complex_name: Lamin A/C-emerin complex
genes_covered: LMNA;EMD
entities: LMNA domain, EMD domain
```

So your “52-ish genes” may become maybe **30–55 AF3 jobs**, depending on shared complexes and domain trimming. With 30 jobs/day, that is still roughly **2 production days**.

---

# Phase 4 — Retrieve sequences

Codex only has the gene symbol from the matrix, so it needs to retrieve canonical protein sequences.

For each target gene and partner gene:

1. Query UniProt for **reviewed human Swiss-Prot** entry.
2. Prefer canonical isoform.
3. Save sequence, accession, protein length, and isoform notes.
4. Warn if multiple reviewed entries exist.
5. Warn if the target variant protein position does not match the UniProt canonical sequence.

Create:

```text
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
```

Columns:

```text
gene
uniprot_accession
entry_name
protein_name
organism
reviewed
canonical_sequence
sequence_length
isoform_note
retrieval_date
retrieval_source
manual_review_needed
```

Also create FASTA:

```text
datasets/feature_sources/alphafold3/interim/af3_all_entities.fasta
```

---

# Phase 5 — Build AF3 job manifest

Create:

```text
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv
```

Columns:

```text
job_id
complex_group_id
job_name
genes_covered
entity_1_type
entity_1_name
entity_1_sequence_or_id
entity_1_range
entity_2_type
entity_2_name
entity_2_sequence_or_id
entity_2_range
entity_3_type
entity_3_name
entity_3_sequence_or_id
entity_3_range
ligands
ions
total_tokens_estimated
full_length_or_domain
reason_for_job
evidence_summary
manual_ready_for_submission
submission_day
submission_order
af3_status
af3_output_zip_path
```

Rules:

```text
No job should exceed 5,000 tokens.
Full-length is preferred when feasible.
Domain-level jobs are allowed when full-length is impossible or biologically inappropriate.
Do not submit unrelated proteins together just to save quota.
Every job must represent one plausible biological complex.
Every gene must have a coverage_status, even if no curated partner is found.
```

The AF3 input format supports proteins, RNA, DNA, and ligands; ions are treated as ligands in the standalone AlphaFold 3 JSON format, for example magnesium can be specified with `ccdCodes: ["MG"]`. ([GitHub][2])

For the server UI, you can still submit through the interface and download the `job_request.json` afterward. AlphaFold Server output includes the job request JSON, five mmCIF structures, summary confidence JSON files, and full-data JSON files containing detailed confidence metrics such as full PAE data. ([EMBL-EBI][3])

---

# Phase 6 — Production schedule

Assuming ~52 jobs:

## Preparation days 1–2

Codex should generate these files:

```text
results/af3/gene_panel_audit.tsv
results/af3/af3_full_panel_partner_manifest.tsv
results/af3/af3_complex_groups.tsv
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv
datasets/feature_sources/alphafold3/jobs/manual_submission_packets/
```

Each manual submission packet should contain:

```text
job_name.txt
entities.fasta
ligands_and_ions.txt
notes_for_server_ui.txt
expected_chain_map.tsv
```

## Production day 1

Submit first 30 AF3 jobs.

Download each output zip immediately.

Save to:

```text
datasets/feature_sources/alphafold3/raw/{job_id}/
```

## Production day 2

Submit remaining jobs.

Download each output zip.

Update:

```text
datasets/feature_sources/alphafold3/interim/af3_completed_jobs.tsv
```

---

# Phase 7 — Local feature extraction engine

AlphaFold gives you the structure and confidence outputs. Your local code turns those outputs into CatBoost-ready features.

For each AF3 job, parse:

```text
fold_<job_name>_model_0.cif
fold_<job_name>_summary_confidences_0.json
fold_<job_name>_full_data_0.json
fold_<job_name>_job_request.json
```

The server ranks structures using confidence metrics including pTM and ipTM, and its JSON summaries include chain-level and chain-pair confidence metrics such as chain-pair ipTM and chain-pair PAE summaries. ([EMBL-EBI][3])

Create residue-level features:

```text
af3_complex_group_id
af3_job_id
gene
uniprot_accession
chain_id
residue_number
residue_aa
local_plddt
min_distance_to_any_partner_A
min_distance_to_partner_gene_A
interface_any_partner_4A
interface_any_partner_5A
interface_any_partner_8A
contact_count_any_partner_5A
min_distance_to_ion_A
near_ion_5A
min_distance_to_ligand_A
near_ligand_5A
chain_pair_iptm
chain_pair_pae_min
local_interface_pae_mean
af3_model_rank
af3_complex_confidence_flag
```

Then merge onto variants by:

```text
gene
uniprot_accession
protein_position
reference_amino_acid
```

Create:

```text
datasets/feature_sources/alphafold3/processed/af3_variant_features.tsv
```

Columns:

```text
variant_id
gene
hgvsp
protein_position
ref_aa
alt_aa
af3_mapping_status
af3_complex_available
af3_complex_group_id
af3_best_job_id
af3_residue_in_any_interface_5A
af3_min_distance_to_any_partner_A
af3_min_distance_to_best_partner_A
af3_best_partner_gene
af3_near_ion_5A
af3_near_ligand_5A
af3_local_plddt
af3_local_interface_pae_mean
af3_chain_pair_iptm
af3_feature_missing_reason
```

Missing reasons:

```text
not_missense
no_protein_position
no_uniprot_mapping
ref_aa_mismatch
gene_not_in_af3_panel
no_curated_partner
domain_not_covering_residue
af3_job_failed
low_confidence_complex
```

This preserves the same logic as your existing feature-source design: every structural feature needs explicit status and missing-reason fields. Your project already treats protein mapping errors and FoldX mismatches this way, which is important because FoldX values were not silently imputed. 

---

# Phase 8 — Model comparison

Do not make the primary paper claim depend on global AUROC moving by a lot. Since AF3 features will still be missing for many rows, the strongest analysis is:

```text
Baseline CatBoost
vs
Baseline CatBoost + AF3 full-panel complex features
```

Report this in three layers:

| Analysis                        | Why                                                   |
| ------------------------------- | ----------------------------------------------------- |
| Full matrix                     | Shows whether AF3 features help overall               |
| AF3-covered subset              | Tests whether AF3 helps where it actually applies     |
| Case-based mechanistic analysis | Shows variants with low ΔΔG but high interface signal |

The most publishable group is:

```text
Pathogenic variants with low FoldX ΔΔG but high AF3 interface signal
VUS predicted pathogenic-like and located at AF3-predicted interfaces
CatBoost/AlphaMissense discordant variants explained by interface proximity
Benign variants structurally mapped but far from interfaces
```

This is exactly the gap in your current model: FoldX ΔΔG asks whether the mutation destabilizes the fold, while AF3 complex mapping asks whether the residue sits at a biologically meaningful interaction surface.

---

# Updated Codex instruction block

Give Codex something like this:

```text
We are changing the AlphaFold 3 plan.

Do not prioritize only 5–10 genes. Build a full-panel AF3 complex-mapping workflow for the entire cardiogenetic gene panel.

Objective:
For every gene in the local variant matrix, assign at least one AF3 coverage status:
1. direct protein-protein complex,
2. shared multi-gene complex,
3. protein-ligand/ion context,
4. domain-level complex because full-length is not feasible,
5. no curated partner found.

Use the current best local matrix:
datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_vep.tsv

First, recover the exact gene list from the matrix and generate:
results/af3/gene_panel_audit.tsv

Then build:
results/af3/af3_full_panel_partner_manifest.tsv
results/af3/af3_complex_groups.tsv
datasets/feature_sources/alphafold3/interim/uniprot_gene_to_sequence.tsv
datasets/feature_sources/alphafold3/jobs/af3_job_manifest.tsv

Do not use a priority score to select genes. Every gene must be reviewed.

For each gene, retrieve canonical reviewed human UniProt sequence. Then identify known partners using, in order:
1. manual cardiac complex whitelist,
2. UniProt curated interactions,
3. Reactome complexes/pathways,
4. IntAct or other curated physical-interaction evidence,
5. PDB/RCSB co-complex evidence,
6. STRING only as weak backup requiring manual review.

Group genes into shared biological complexes when appropriate, so one AF3 job can cover multiple genes. Example:
TNNI3 + TNNT2 + TNNC1 + TPM1 + ACTC1 + Ca2+ can cover multiple thin-filament genes.

Each AF3 job must represent one plausible biological complex. Do not combine unrelated proteins only to save quota.

Each job must stay under 5,000 tokens. Full-length proteins are preferred when feasible, but domain-level jobs are allowed for very large proteins or when biology supports a specific interaction domain.

For each job, create a manual submission packet:
job_name.txt
entities.fasta
ligands_and_ions.txt
notes_for_server_ui.txt
expected_chain_map.tsv

After outputs are downloaded, parse mmCIF and JSON files locally and generate residue-level and variant-level features:
af3_residue_features.tsv
af3_variant_features.tsv

Feature extraction should calculate:
- local pLDDT
- minimum distance to partner chain
- interface residue flags at 4, 5, and 8 Å
- contact counts
- minimum distance to ligand/ion
- chain-pair ipTM
- chain-pair PAE summaries
- local interface PAE where possible

Merge AF3 features back to the CatBoost matrix by gene, UniProt accession, protein position, and reference amino acid. Require reference-AA QC. If the reference amino acid does not match the UniProt sequence, mark af3_mapping_status = ref_aa_mismatch and do not impute.

Output all missingness explicitly with af3_feature_missing_reason.
```

## Bottom line

Your revised plan should be:

```text
Full-panel AF3 complex coverage
not
prioritized AF3 pilot
```

The important safeguard is that **full-panel** does not mean “throw random proteins into AF3.” It means every gene is assigned a curated biological context, and the final output is a structured interaction-aware feature layer for your 85,677-row CatBoost matrix.

[1]: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/ "AlphaFold Server: Your gateway to AlphaFold 3 | AlphaFold"
[2]: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md "alphafold3/docs/input.md at main · google-deepmind/alphafold3 · GitHub"
[3]: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/interpreting-results-from-alphafold-server/ "Interpreting results from AlphaFold Server | AlphaFold"


Yes. The **normal protein-protein complex structures** are helpful for mutated variants because they let you ask:

> “Is this mutated residue located in a biologically important place?”

You are **not** directly asking AF3:

> “Does this exact mutation destroy binding?”

Unless you run the mutant complex too.

Instead, for the full 52-ish gene panel, the broad scalable version is:

> “Map each variant onto a normal/wild-type AF3 complex and extract whether that residue sits at an interface, ligand site, ion site, buried core, low-confidence region, flexible region, or high-confidence functional surface.”

That gives you **variant-level structural context features** for CatBoost.

Your current FoldX/structure layer is sparse: your manuscript table says the final matrix has **85,677 rows**, but FoldX ΔΔG covers only **313 / 85,677 rows, 0.37%**. Your current model is mostly powered by VEP, SpliceAI, gnomAD, dbNSFP, AlphaMissense, conservation, and ClinGen-type features, while protein-structure features are useful but sparse. 

## Why wild-type complexes help mutated variants

Think of it like this:

```text
Variant A: KCNQ1 p.ArgXXXGln
```

Without AF3, CatBoost may know:

```text
gene = KCNQ1
consequence = missense
AlphaMissense = high
REVEL = high
gnomAD AF = rare
domain = maybe known
```

With AF3 complex mapping, CatBoost can also know:

```text
This residue is 3.8 Å from KCNE1
This residue is at the KCNQ1-KCNE1 interface
This residue is 6.2 Å from calmodulin
This region has high pLDDT
The KCNQ1-KCNE1 interface has good inter-chain confidence
The residue is exposed at a partner-binding surface
```

That is valuable because many pathogenic variants do **not** necessarily destabilize the whole protein. They may instead disrupt:

```text
protein-protein binding
ion binding
ligand binding
channel gating regions
sarcomere complex assembly
nuclear-envelope interactions
calmodulin regulation
thin-filament regulation
```

So AF3 complex mapping gives you a different biological mechanism than ΔΔG.

## ΔΔG vs AF3 complex features

| Question                                                         | Feature family              | Meaning                  |
| ---------------------------------------------------------------- | --------------------------- | ------------------------ |
| Does the mutation destabilize the folded protein?                | FoldX ΔΔG                   | Stability effect         |
| Is the residue buried or exposed?                                | FreeSASA / DSSP / structure | Structural environment   |
| Is the residue in a functional domain?                           | UniProt / ClinGen / VEP     | Known biological context |
| Is the residue at a partner interface?                           | AF3 complex mapping         | Interaction context      |
| Is the residue near Ca²⁺, ATP, Mg²⁺, DNA/RNA, or another entity? | AF3 ligand/ion proximity    | Functional-site context  |
| Is the predicted interface reliable?                             | pLDDT / PAE / ipTM          | Confidence/QC            |

AF3 outputs include predicted structures plus confidence files. AlphaFold Server output includes predicted `.cif` structures, summary confidence JSON files, full-data JSON files with detailed metrics such as full PAE, and the original job request JSON. ([EMBL-EBI][1]) AlphaFold 3 confidence outputs include pLDDT, PAE, pTM, and ipTM-type metrics; pLDDT is a local per-atom confidence score, while PAE helps assess relative positional confidence between parts of the model. ([GitHub][2])

## What features you can extract locally for CatBoost

### 1) Basic AF3 coverage/status features

These are important because CatBoost handles missingness, but you need explicit missingness flags.

```text
af3_complex_available
af3_complex_group_id
af3_best_complex_name
af3_n_complexes_covering_gene
af3_n_complexes_covering_residue
af3_mapping_status
af3_feature_missing_reason
af3_domain_level_model_flag
af3_full_length_model_flag
```

Useful missing reasons:

```text
not_missense
no_hgvsp
no_protein_position
no_uniprot_mapping
ref_aa_mismatch
gene_not_modelled
residue_outside_modelled_domain
af3_job_failed
low_confidence_complex
no_curated_partner
```

These prevent silent garbage features.

---

### 2) Local residue confidence features

From AF3 pLDDT:

```text
af3_local_plddt
af3_local_plddt_window_3aa
af3_local_plddt_window_5aa
af3_local_plddt_window_11aa
af3_low_confidence_region_flag
af3_high_confidence_region_flag
```

Why useful:

A variant in a high-confidence structured interface is easier to interpret than a variant in a low-confidence disordered tail. Low pLDDT does **not** mean pathogenic. It often means disorder/flexibility. So CatBoost should learn its meaning rather than you manually overinterpret it.

---

### 3) Protein-protein interface features

These are probably the highest-value AF3 features.

For each variant residue:

```text
af3_in_any_partner_interface_4A
af3_in_any_partner_interface_5A
af3_in_any_partner_interface_8A
af3_min_distance_to_any_partner_A
af3_contact_count_any_partner_5A
af3_contact_count_any_partner_8A
af3_partner_surface_contact_area
af3_interface_partner_count
```

Partner-specific examples:

```text
af3_kcnq1_kcne1_interface_5A
af3_kcnq1_kcne1_min_distance_A
af3_kcnq1_calm_interface_5A
af3_scn5a_calm_interface_5A
af3_lmna_emd_interface_5A
af3_tnni3_tnnt2_interface_5A
```

Why useful:

A missense variant at a physical interaction surface may be pathogenic even if the protein fold is stable. This is exactly the mechanism that ΔΔG can miss.

---

### 4) Interface confidence features

A residue may appear close to a partner, but the interface may be low-confidence. You need confidence-aware features.

```text
af3_complex_iptm
af3_complex_ptm
af3_chain_pair_iptm
af3_chain_pair_pae_mean
af3_chain_pair_pae_min
af3_local_interface_pae_mean
af3_local_interface_pae_median
af3_interface_confident_flag
af3_interface_low_confidence_flag
```

A good derived feature:

```text
af3_confident_interface_residue =
    af3_in_any_partner_interface_5A == true
    AND af3_local_plddt >= 70
    AND af3_local_interface_pae_mean <= threshold
```

This is safer than using distance alone.

---

### 5) Ion and ligand proximity features

For jobs with calcium, magnesium, ATP, heme, etc.

```text
af3_near_any_ion_4A
af3_near_any_ion_5A
af3_min_distance_to_any_ion_A
af3_near_ca_5A
af3_min_distance_to_ca_A
af3_near_mg_5A
af3_min_distance_to_mg_A
af3_near_atp_5A
af3_min_distance_to_atp_A
af3_near_ligand_5A
af3_min_distance_to_ligand_A
```

Why useful:

For genes like **TNNC1/CALM/RYR2/KCNQ1/MYH7**, ion or ligand proximity may capture calcium regulation, ATPase biology, or channel-related functional mechanisms.

---

### 6) Buried/exposed/interface environment features

From AF3 coordinates, DSSP, FreeSASA, or your own surface calculations:

```text
af3_residue_sasa
af3_relative_sasa
af3_buried_flag
af3_exposed_flag
af3_interface_and_buried_flag
af3_interface_and_exposed_flag
af3_secondary_structure
af3_helix_flag
af3_sheet_flag
af3_loop_flag
af3_transmembrane_region_flag
```

Why useful:

A mutation at a buried core can destabilize structure. A mutation on an exposed interface can disrupt binding. A mutation in a flexible/disordered region may need different interpretation.

---

### 7) Mutation-aware structural features

The wild-type complex tells you where the residue is. The variant tells you what changed.

So combine AF3 location with amino-acid substitution properties:

```text
aa_ref
aa_alt
aa_charge_change
aa_polarity_change
aa_aromatic_change
aa_size_change
aa_hydrophobicity_change
aa_grantham_distance
aa_blosum62_score
```

Then interaction-specific hybrid features:

```text
af3_interface_charge_change
af3_interface_hydrophobicity_change
af3_interface_size_change
af3_interface_to_proline
af3_interface_to_glycine
af3_buried_large_to_small
af3_buried_small_to_large
af3_near_ca_charge_change
af3_near_ligand_polarity_change
```

These are important because **wild-type AF3 alone does not know the mutant effect**. The mutation-aware feature is what tells CatBoost:

> “This interface residue changed from a positively charged arginine to a neutral glutamine.”

That is much more meaningful than just:

> “Residue is at interface.”

---

### 8) Domain + complex overlap features

Your model already uses domains and protein-position context, but AF3 lets you combine domain and interaction context.

```text
af3_interface_in_uniprot_domain
af3_interface_in_transmembrane_domain
af3_interface_in_channel_pore_region
af3_interface_in_calmodulin_binding_region
af3_interface_in_sarcomere_domain
af3_interface_in_lamin_tail_domain
af3_domain_distance_to_partner_A
af3_variant_domain_complex_name
```

Example:

```text
KCNQ1 variant in transmembrane/channel region
AND near KCNE1 interface
AND rare in gnomAD
```

This is much more biologically informative than gene name alone.

---

### 9) Multi-complex aggregation features

Some proteins will appear in more than one AF3 complex. For example, KCNQ1 might be modeled with KCNE1 and calmodulin.

For each variant, aggregate across all complexes that cover that residue:

```text
af3_n_interfaces_hit
af3_min_distance_to_any_partner_across_complexes
af3_best_interface_confidence
af3_any_confident_interface_hit
af3_best_partner_gene
af3_best_complex_group
af3_mean_local_plddt_across_complexes
af3_max_contact_count_across_complexes
```

This makes the feature table cleaner because CatBoost does not need 500 partner-specific columns for every gene.

You can still keep partner-specific columns for key genes, but the model-facing table should include general features too.

---

### 10) Comparison features against FoldX / AlphaMissense

These may be your most publishable mechanistic features.

```text
af3_interface_high_and_foldx_ddg_low
af3_interface_high_and_alphamissense_high
af3_interface_high_and_alphamissense_low
af3_not_interface_but_foldx_ddg_high
af3_pathogenic_like_discordance_flag
```

Examples:

```text
Variant has low FoldX ΔΔG but lies at a confident protein-protein interface.
Variant has high AlphaMissense but no obvious stability/interface signal.
Variant is VUS, CatBoost pathogenic-like, and lies near Ca2+ or partner interface.
```

Your pasted AF3 plan already framed this as the most publishable angle: VUS predicted pathogenic-like near AF3 interfaces, CatBoost/AlphaMissense discordant variants, and pathogenic variants with low ΔΔG but high interface signal. 

## How helpful will this be?

I would be realistic:

### It may not massively improve full-matrix AUROC

Because many of your 85,677 variants are not clean missense/protein-position-mapped variants. Your current model already performs very well externally, and your strongest features are already broad-coverage features like VEP, SpliceAI, gnomAD, dbNSFP, AlphaMissense, and ClinGen. 

So I would **not** promise:

```text
AF3 will increase AUROC from 0.977 to 0.99.
```

### It can be very helpful in the AF3-covered missense subset

The best test is:

```text
Baseline CatBoost
vs
Baseline + AF3 features
```

in three datasets:

```text
1. Full matrix
2. Missense/protein-position-mapped subset
3. AF3-covered interface/complex subset
```

The likely strongest improvement is in subset 2 and 3, especially for:

```text
missense variants
VUS rows
rare variants
dbNSFP/AlphaMissense discordant rows
FoldX-missing rows
variants in channel/sarcomere/calmodulin/nuclear-envelope complexes
```

### It may help interpretability more than raw discrimination

This is important. Even if AUROC barely changes, AF3 could still strengthen the paper by explaining **why** certain variants are high-priority.

For example:

```text
CatBoost calls a VUS pathogenic-like.
AF3 shows it lies at a confident KCNQ1-KCNE1 interface.
FoldX ΔΔG is low.
Interpretation: the likely mechanism may be interaction disruption, not global protein destabilization.
```

That is clinically and biologically useful.

## The best CatBoost feature set to start with

I would start with this compact model-facing set:

```text
af3_complex_available
af3_mapping_status
af3_n_complexes_covering_residue
af3_any_interface_5A
af3_any_interface_8A
af3_min_distance_to_any_partner_A
af3_contact_count_any_partner_5A
af3_best_partner_gene
af3_best_complex_group
af3_local_plddt
af3_local_plddt_window_5aa
af3_chain_pair_iptm_best
af3_local_interface_pae_mean_best
af3_confident_interface_flag
af3_near_any_ion_5A
af3_min_distance_to_any_ion_A
af3_near_any_ligand_5A
af3_min_distance_to_any_ligand_A
af3_relative_sasa
af3_buried_flag
af3_exposed_flag
af3_secondary_structure
af3_interface_charge_change
af3_interface_size_change
af3_interface_hydrophobicity_change
af3_buried_size_change
af3_near_ion_charge_change
af3_interface_high_and_foldx_ddg_low
af3_interface_high_and_alphamissense_discordant
af3_feature_missing_reason
```

Then later, add partner-specific columns for high-value complexes:

```text
af3_kcnq1_kcne1_interface_5A
af3_kcnq1_calm_interface_5A
af3_scn5a_calm_interface_5A
af3_ryr2_fkbp1b_interface_5A
af3_troponin_complex_interface_5A
af3_lmna_emd_interface_5A
```

## The exact local logic

For each AF3 complex:

```text
1. Download AF3 mmCIF + JSON outputs.
2. Parse chains and map chain IDs to genes.
3. For every residue in every target chain:
   - calculate min distance to every partner chain
   - calculate contact count
   - calculate ion/ligand distance
   - extract pLDDT
   - extract PAE/ipTM confidence where possible
4. Create residue-level AF3 table.
5. Merge residue-level table to variant table by:
   gene + UniProt accession + protein position + reference amino acid
6. Add mutation-aware features using ref AA and alt AA.
7. Train CatBoost:
   baseline vs baseline + AF3 features.
```

## Bottom line

The normal AF3 protein-protein complexes are helpful because they create a **map of biological interaction surfaces**. Your variants are then mapped onto that map.

The core feature is not:

```text
AF3 says this mutant breaks binding.
```

The correct feature is:

```text
This variant affects a residue located at a confident protein-protein / ion / ligand interface.
```

Then CatBoost learns whether those interface-aware features add predictive value beyond AlphaMissense, REVEL, VEP consequence, gnomAD rarity, ClinGen context, and FoldX ΔΔG.

[1]: https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/interpreting-results-from-alphafold-server/?utm_source=chatgpt.com "Interpreting results from AlphaFold Server"
[2]: https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md?utm_source=chatgpt.com "alphafold3/docs/output.md at main"
