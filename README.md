# 8-Oxoguanine transcriptional miscoding: reproducibility repository

Companion repository for the manuscript "Disease-associated protein states accessible to 8-oxoguanine transcriptional miscoding without genomic mutation".

This repository rebuilds the publication-facing supplementary tables, regenerates all six figures
and recomputes every numerical result the manuscript reports, starting from the derived data
distributed here. Primary source files are not redistributed.

## Contents

`data/supplementary_data/` holds the 7 publication-facing datasets, which correspond to the
six supplementary table items S1 to S6 because S3 is split into two worksheets:

| dataset | rows | columns | content |
|---|---|---|---|
| S1 | 48 | 6 | Genetic-code consequences of coding-sense C>A changes |
| S2 | 3,910 | 23 | Accessible disease-associated protein states |
| S3a | 21 | 13 | Expression-class and disease-category screens |
| S3b | 15 | 11 | Gene-context model estimates and specification sensitivity |
| S4 | 945 | 16 | Pathway accessibility analysis |
| S5 | 78 | 23 | Directly observed proteomic states |
| S6 | 23 | 17 | Oxidized-guanine map analysis |

The same directory holds `Supplementary_Table_Index.csv` and `data_dictionary.csv`, which describes
every field of every dataset.

`data/analysis_inputs/` holds the frozen per-gene inputs the analysis consumed, namely the Human
Protein Atlas expression record, the MONDO disease-annotation record and the gene-context model
input, together with `PROVENANCE_NOTE.md` describing their scope.

`data/figure_inputs/` holds one frozen JSON record per figure plus the legend records, so figures
regenerate without re-querying any external resource.

`data/provenance/` holds `SOURCE_MANIFEST.tsv`, which identifies every primary resource by
accession, release and retrieval date, the strand interpretation of each oxidized-guanine map
resource, and the ClinVar completeness audit inputs.

`code/analysis/` holds the loading, claim-recomputation and quality-control modules. `code/figures/`
holds one builder per figure plus the shared style module.

## Figures reproduced

Figure 1 genetic code and chemistry, Figure 2 disease convergence, Figure 3 gene context,
Figure 4 proteogenomic observation, Supplementary Figure S1 specification sensitivity and
Supplementary Figure S2 oxidized-guanine map evidence.

## What this repository can and cannot recompute

It recomputes 50 reported quantities directly from the distributed datasets, rebuilds the
workbook from the canonical delimited tables, regenerates every figure from its frozen input and
runs structural, dictionary, semantic, retired-reference, figure-geometry and checksum checks.

It cannot recompute six quantities, which are reported in the claim audit as not machine checkable.
These are counts published by the source proteogenomic study, intermediates of the extraction step
that reads primary source files, and the route-matched background comparison, which is enumerated
over all retained transcripts rather than over the distributed datasets. The eight single-cell
expression classes of Figure 3a and the ten top-level disease categories of the category screen
regenerate from their frozen per-class records rather than from per-gene memberships, because those
memberships were not retained. `data/analysis_inputs/PROVENANCE_NOTE.md` states this precisely.

## Creating the environment

    conda env create -f environment.yml
    conda activate oxog_tm

## Regenerating everything

    python reproduce_all.py

Individual stages are available as `python reproduce_figures.py` and `python run_qc.py`.

Figures use Arial when it is installed and DejaVu Sans otherwise. Set `OXOG_FIGURE_FONT` to force a
family, which is how the fallback path is exercised.

## Expected output

    1. supplementary datasets
       7 datasets, 5,040 rows
    2. data dictionary
       109 fields described
    3. consolidated workbook
       107,961 cells written, 0 mismatches
    4. figures
       6 figures regenerated, 0 geometry failures
    5. quality control
       PASS
    PASS

A non-zero exit status means at least one check failed, and `audit/QC_REPORT.txt` names every
failure.
