# ClinVar and pathogenicity

**What it is:** ClinVar is NCBI's public database where labs and researchers submit variants with a claimed clinical significance. "Pathogenic" means the evidence supports the variant causing disease.

**Biology / intuition:** Most variants are harmless. A pathogenic one disrupts a gene's function (e.g. a BRCA1 frameshift raises breast/ovarian cancer risk). Classification follows the ACMG/AMP guidelines (Richards et al. 2015) into five tiers: Pathogenic, Likely pathogenic, Uncertain significance (VUS), Likely benign, Benign. Submitters can **disagree**, and ClinVar records a review status (stars) showing how much evidence and consensus there is.

**When it's wrong to use:**
- A hit in ClinVar is evidence, not a diagnosis. Context matters (zygosity, family history, inheritance).
- Many variants are VUS. "Not in ClinVar" ≠ benign.
- Entries get reclassified over time, so record the ClinVar release date.
- Different genome builds (GRCh37 vs GRCh38) give different positions. Matching by rsID alone can mislead for multi-allelic sites.

**Where in our code:** `KNOWN_VARIANTS` in [variant_annotation.py](../../src/bioagent/pipelines/variant_annotation.py): a hardcoded dictionary of 8 rsIDs, **not** a live ClinVar query. Every entry says `missense_variant`, which I have not verified against ClinVar. The validation rule says to check them before claiming "ClinVar annotation".

**Likely interview questions**
1. *How do you annotate variants with clinical significance?* — Join on the variant's rsID/position to ClinVar, or run Ensembl VEP/ANNOVAR. My pipeline currently uses a small curated table and I'm replacing it with a real ClinVar extract. *(Earn in W4.)*
2. *What does "likely pathogenic" mean?* — Under ACMG/AMP it means about ≥90% certainty the variant is disease-causing, a tier below fully pathogenic. VUS means the evidence is insufficient either way.
