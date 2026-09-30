# Analysis-input provenance and reproducibility scope

## Files in this directory

hpa_expression_class_membership.tsv
    Exact frozen per-gene expression record used by the analysis, covering 20,151 genes with the
    post-mitotic and proliferative expression scores, their log ratio and the broad expression class
    assigned from them. Verified to reproduce the post_mitotic_enriched covariate of the
    gene-context model for every one of the 4,247 model genes with perfect agreement. The source
    release and retrieval date were not recorded at original retrieval and are not reconstructed
    here.

mondo_gene_condition_flags.tsv
    Exact frozen per-gene disease-annotation record used by the analysis, covering 4,401 genes with
    the cancer, hereditary cancer syndrome and neurodegenerative condition flags derived from MONDO
    terms mapped through ClinVar phenotype annotations, together with the broad expression class.
    This is the input to the neurological_disease_gene covariate of the gene-context model.

gene_context_model_input.tsv
    Per-gene accessible and eligible state counts with the model covariates, for all 4,247 genes
    entering the gene-context model. This is analysis input rather than a reader-facing result
    table.

## What these files reproduce, and what they do not

The broad three-class expression screen, the post-mitotic and neurological covariates and the
gene-context interaction model are fully reproducible from the files above.

Two per-gene assignments used in screens reported in this study are not recoverable and are
therefore not distributed. The eight single-cell expression classes shown in Figure 3a were derived
from the Human Protein Atlas single-cell type consensus resource in a working directory that was
subsequently lost, and the frozen record above carries the broad class rather than the per-cell-type
assignments. The ten top-level MONDO disease categories used in the category screen were derived in
the same working directory, and the frozen record above carries three condition flags rather than
the ten category assignments. Both screens are therefore distributed at the level of the computed
per-class and per-category result, in Supplementary Table S3a, and Figure 3a regenerates from that
frozen record rather than from per-gene memberships. Re-deriving either assignment from current
releases of those resources would change the annotation vintage and would not reproduce the
analysis as run.
