"""Recompute every numerical quantity the manuscript reports."""
from __future__ import annotations

from decimal import ROUND_HALF_UP, Decimal

import pandas as pd

from . import tables
from .tables import as_bool, load

NEURONAL = r"synap|neurotransmitter|NMDA|neuron|axon|dendrit|glutamatergic|GABAergic"
NOT_NEURONAL = r"\bDNA\b|meiot|meiosis|recombinat|strand exchange|homologous"
HEMOGLOBIN = ["HBA2", "HBB", "HBD"]


def _round(value, places="0.01") -> str:
    """Conventional round-half-up to the reported precision."""
    return str(Decimal(str(value)).quantize(Decimal(places), rounding=ROUND_HALF_UP))


def neuronal_mask(frame):
    """The prespecified synaptic or neuronal gene-set classification."""
    return (frame.pathway_name.str.contains(NEURONAL, case=False, na=False)
            & ~frame.pathway_name.str.contains(NOT_NEURONAL, case=False, na=False))


def recompute() -> pd.DataFrame:
    """Compare every reported quantity against the distributed datasets."""
    s1, s2, s3a, s3b = load("S1"), load("S2"), load("S3a"), load("S3b")
    s4, s5, s6 = load("S4"), load("S5"), load("S6")
    model = tables.analysis_input("gene_context_model_input.tsv")

    consequence = s1.consequence_class.str.lower()
    missense = s1[consequence == "missense"]
    pairs = set(zip(missense.reference_amino_acid, missense.product_amino_acid))
    residues = set(missense.reference_amino_acid) | set(missense.product_amino_acid)
    one_way = sum(1 for r in residues
                  if (r in missense.reference_amino_acid.values)
                  != (r in missense.product_amino_acid.values))

    higher = as_bool(s2.higher_confidence_subset)
    detected = as_bool(s2.any_map_detection)
    route = s2.route_class

    neurons = s3a[s3a.screen.str.contains("single-cell", case=False, na=False)
                  & s3a.term.str.contains("neuron", case=False, na=False)].iloc[0]
    nervous = s3a[s3a.screen.str.contains("disease ontology", case=False, na=False)
                  & s3a.term.str.contains("nervous", case=False, na=False)].iloc[0]

    def interaction(stratum):
        rows = s6[s6.analysis.str.contains("interaction", na=False)
                  & s6.stratum.str.contains(stratum, case=False, na=False)]
        return rows.iloc[0]

    primary, supported, confident = (interaction("primary"),
                                     interaction("oxidized-guanine-supported subset"),
                                     interaction("higher-confidence supported"))
    neuronal = neuronal_mask(s4)
    significant_neuronal = s4[neuronal & (s4.q_value < 0.05)]

    checks = [
        ("C01", "C>A codon-position events", 48, len(s1)),
        ("C02", "missense outcomes", 33, int((consequence == "missense").sum())),
        ("C03", "synonymous outcomes", 11, int((consequence == "synonymous").sum())),
        ("C04", "stop-codon outcomes", 4, int(consequence.str.contains("stop").sum())),
        ("C05", "directional amino-acid substitutions", 18, len(pairs)),
        ("C06", "replaceable residues", 11, missense.reference_amino_acid.nunique()),
        ("C07", "introducible residues", 13, missense.product_amino_acid.nunique()),
        ("C08", "residues accessible in one direction only", 8, one_way),
        ("C09", "accessible disease-associated states", 3910, len(s2)),
        ("C10", "genes carrying an accessible state", 1333, s2.gene.nunique()),
        ("C11", "percentage of eligible missense states", "5.7",
         _round(100 * model.accessible_states.sum() / model.eligible_states.sum(), "0.1")),
        ("C12", "higher-confidence states", 977, int(higher.sum())),
        ("C13", "higher-confidence genes", 403, s2[higher].gene.nunique()),
        ("C14", "exact nucleotide route", 3115, int(route.str.contains("xact").sum())),
        ("C15", "same-position alternate-base route", 584,
         int(route.str.contains("ame-position").sum())),
        ("C16", "different-position route", 211, int(route.str.contains("ifferent").sum())),
        ("C17", "different-position genes", 172,
         s2[route.str.contains("ifferent")].gene.nunique()),
        ("C18", "different-position higher-confidence", 48,
         int((route.str.contains("ifferent") & higher).sum())),
        ("C19", "phenylalanine to leucine states", 173,
         int(((s2.reference_amino_acid == "F") & (s2.alternate_amino_acid == "L")
              & route.str.contains("ifferent")).sum())),
        ("C20", "serine to arginine states", 38,
         int(((s2.reference_amino_acid == "S") & (s2.alternate_amino_acid == "R")
              & route.str.contains("ifferent")).sum())),
        ("C21", "neurons accessible fraction", "7.67", _round(neurons.observed_pct)),
        ("C22", "neurons permutation expectation", "6.03", _round(neurons.null_pct)),
        ("C23", "neurons adjusted P value", "0.004", _round(neurons.q_value, "0.001")),
        ("C24", "nervous-system category adjusted P value", "0.21", _round(nervous.q_value)),
        ("C25", "joint interaction odds ratio", "1.49", _round(primary.odds_ratio)),
        ("C26", "joint interaction interval", "1.08 to 2.06",
         f"{_round(primary.ci_lower)} to {_round(primary.ci_upper)}"),
        ("C27", "joint interaction P value", "0.016", _round(primary.p_value, "0.001")),
        ("C28", "estimable gene sets", 945, int(s4.q_value.notna().sum())),
        ("C29", "gene sets passing correction", 102, int((s4.q_value < 0.05).sum())),
        ("C30", "enriched gene sets", 64,
         int(((s4.q_value < 0.05) & (s4.odds_ratio > 1)).sum())),
        ("C31", "depleted gene sets", 38,
         int(((s4.q_value < 0.05) & (s4.odds_ratio < 1)).sum())),
        ("C32", "synaptic or neuronal gene sets tested", 30, int(neuronal.sum())),
        ("C33", "synaptic or neuronal gene sets passing correction", 8,
         len(significant_neuronal)),
        ("C34", "synaptic or neuronal significant sets that are enriched", 8,
         int((significant_neuronal.odds_ratio > 1).sum())),
        ("C35", "map-detected accessible coordinates", 2146, int(detected.sum())),
        ("C36", "map-detected higher-confidence coordinates", 559,
         int((detected & higher).sum())),
        ("C37", "map-supported interaction odds ratio", "1.77", _round(supported.odds_ratio)),
        ("C38", "map-supported interaction interval", "1.22 to 2.57",
         f"{_round(supported.ci_lower)} to {_round(supported.ci_upper)}"),
        ("C39", "map-supported interaction P value", "0.0027",
         _round(supported.p_value, "0.0001")),
        ("C40", "higher-confidence supported odds ratio", "2.87", _round(confident.odds_ratio)),
        ("C41", "higher-confidence supported interval", "1.42 to 5.80",
         f"{_round(confident.ci_lower)} to {_round(confident.ci_upper)}"),
        ("C42", "observed proteomic states", 78, len(s5)),
        ("C43", "genes carrying observed states", 44, s5.gene.nunique()),
        ("C44", "observed states in hemoglobin genes", 25, int(s5.gene.isin(HEMOGLOBIN).sum())),
        ("C45", "endogenous lesion-support states", 29,
         int((s5.oxidized_guanine_support_category == "endogenous").sum())),
        ("C46", "induced-only lesion-support states", 9,
         int((s5.oxidized_guanine_support_category == "induced-only").sum())),
        ("C47", "states with no map detection", 39,
         int((s5.oxidized_guanine_support_category == "not detected").sum())),
        ("C48", "states with an unresolved lesion coordinate", 1,
         int((s5.oxidized_guanine_support_category == "unresolved").sum())),
        ("C49", "genes in the gene-context model", 4247, len(model)),
        ("C50", "accessible states in model genes", 3836, int(model.accessible_states.sum())),
    ]
    rows = [dict(claim_id=i, claim=text, expected_value=str(expected),
                 recomputed_value=str(observed),
                 status="PASS" if str(expected) == str(observed) else "FAIL")
            for i, text, expected, observed in checks]
    for claim_id, text, expected, note in [
        ("C51", "confidently localized sites in the proteogenomic source", 1955,
         "published count from the source study, which this archive does not redistribute"),
        ("C52", "source states reconciled to MANE Select", 1658,
         "intermediate of the source extraction, retained in provenance"),
        ("C53", "substitution-class compatible source states", 93,
         "intermediate of the source extraction, retained in provenance"),
        ("C54", "route-matched background positions", 638965,
         "enumerated over all retained transcripts by the extraction step"),
        ("C55", "background positions with map detection", 346415,
         "enumerated by the extraction step, retained in provenance"),
        ("C56", "background comparison odds ratio, interval and P value",
         "1.03 (0.96 to 1.09), P = 0.40",
         "computed from the two counts above, which are not recomputable from the distributed "
         "datasets alone")]:
        rows.append(dict(claim_id=claim_id, claim=text, expected_value=str(expected),
                         recomputed_value="not recomputable here",
                         status="NOT_MACHINE_CHECKABLE", notes=note))
    frame = pd.DataFrame(rows)
    if "notes" not in frame:
        frame["notes"] = ""
    return frame.fillna({"notes": ""}).sort_values("claim_id").reset_index(drop=True)
