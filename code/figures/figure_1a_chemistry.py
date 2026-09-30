"""Chemically verified base-pairing schematic for Figure 1a.

Structures come from PubChem canonical SMILES (guanine CID 135398634,
8-oxoguanine CID 135420630, adenine CID 190, cytosine CID 597) and are
verified by graph analysis rather than drawn from memory. Geometry follows
the Pol II elongation-complex structures of Damsma and Cramer
(J Biol Chem 2009, PDB 3I4M Watson-Crick and 3I4N Hoogsteen).
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

SMILES = {
    "guanine": "C1=NC2=C(N1)C(=O)NC(=N2)N",
    "8-oxoguanine": "C12=C(NC(=O)N1)N=C(NC2=O)N",
    "cytosine": "C1=C(NC(=O)N=C1)N",
    "adenine": "C1=NC2=NC=NC(=C2N1)N",
}


def purine_atoms(mol):
    """Assign purine numbering by ring topology, not by atom order."""
    rings = mol.GetRingInfo().AtomRings()
    five = [r for r in rings if len(r) == 5][0]
    six = [r for r in rings if len(r) == 6][0]
    fused = set(five) & set(six)
    out = {}
    c8 = [i for i in five if mol.GetAtomWithIdx(i).GetSymbol() == "C" and i not in fused]
    out["C8"] = c8[0]
    for i in six:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() != "C":
            continue
        oxy = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "O"]
        amine = [n.GetIdx() for n in atom.GetNeighbors()
                 if n.GetSymbol() == "N" and n.GetTotalNumHs() == 2]
        if oxy:
            out["C6"], out["O6"] = i, oxy[0]
        if amine:
            out["C2"], out["N2"] = i, amine[0]
    # C5 is the fused carbon bonded to C6; C4 is the other fused carbon
    c5 = [i for i in fused
          if mol.GetBondBetweenAtoms(i, out["C6"]) is not None]
    out["C5"] = c5[0]
    out["C4"] = [i for i in fused if i != out["C5"]][0]
    for name, anchor in (("N7", "C5"), ("N9", "C4")):
        cand = [n.GetIdx() for n in mol.GetAtomWithIdx(out[anchor]).GetNeighbors()
                if n.GetSymbol() == "N" and n.GetIdx() in five]
        out[name] = cand[0]
    # N1 sits between C2 and C6 in the six-membered ring
    n1 = [i for i in six if mol.GetAtomWithIdx(i).GetSymbol() == "N"
          and mol.GetBondBetweenAtoms(i, out["C6"]) is not None
          and mol.GetBondBetweenAtoms(i, out["C2"]) is not None]
    out["N1"] = n1[0]
    out["N3"] = [i for i in six if mol.GetAtomWithIdx(i).GetSymbol() == "N"
                 and i not in (out["N1"],)][0]
    oxo = [n.GetIdx() for n in mol.GetAtomWithIdx(out["C8"]).GetNeighbors()
           if n.GetSymbol() == "O"]
    out["O8"] = oxo[0] if oxo else None
    return out


def adenine_atoms(mol):
    """Assign adenine numbering; adenine carries an exocyclic amine at C6 and no O6."""
    rings = mol.GetRingInfo().AtomRings()
    five = [r for r in rings if len(r) == 5][0]
    six = [r for r in rings if len(r) == 6][0]
    fused = set(five) & set(six)
    out = {}
    for i in six:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() != "C" or i in fused:
            continue
        amine = [n.GetIdx() for n in atom.GetNeighbors()
                 if n.GetSymbol() == "N" and n.GetTotalNumHs() == 2]
        if amine:
            out["C6"], out["N6"] = i, amine[0]
        else:
            out["C2"] = i
    out["C5"] = [i for i in fused
                 if mol.GetBondBetweenAtoms(i, out["C6"]) is not None][0]
    out["C4"] = [i for i in fused if i != out["C5"]][0]
    out["N1"] = [i for i in six if mol.GetAtomWithIdx(i).GetSymbol() == "N"
                 and mol.GetBondBetweenAtoms(i, out["C6"]) is not None][0]
    out["N3"] = [i for i in six if mol.GetAtomWithIdx(i).GetSymbol() == "N"
                 and i != out["N1"]][0]
    for name, anchor in (("N7", "C5"), ("N9", "C4")):
        out[name] = [n.GetIdx() for n in mol.GetAtomWithIdx(out[anchor]).GetNeighbors()
                     if n.GetSymbol() == "N" and n.GetIdx() in five][0]
    return out


def pyrimidine_atoms(mol, kind):
    """Assign the pairing-edge atoms of cytosine."""
    six = mol.GetRingInfo().AtomRings()[0]
    out = {}
    for i in six:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol() == "C":
            oxy = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "O"]
            amine = [n.GetIdx() for n in atom.GetNeighbors()
                     if n.GetSymbol() == "N" and n.GetTotalNumHs() == 2]
            if oxy:
                out["C2"], out["O2"] = i, oxy[0]
            if amine:
                out["C4"], out["N4"] = i, amine[0]
    ring_n = [i for i in six if mol.GetAtomWithIdx(i).GetSymbol() == "N"]
    out["N3"] = [i for i in ring_n
                 if mol.GetBondBetweenAtoms(i, out["C2"]) is not None
                 and mol.GetBondBetweenAtoms(i, out["C4"]) is not None][0]
    out["N1"] = [i for i in ring_n if i != out["N3"]][0]
    return out


def verify():
    """Confirm the two purines differ only by the C8 carbonyl, with N7-H present."""
    mols = {k: Chem.MolFromSmiles(s) for k, s in SMILES.items()}
    g, o = mols["guanine"], mols["8-oxoguanine"]
    pg, po = purine_atoms(g), purine_atoms(o)
    checks = {
        "guanine formula": Chem.rdMolDescriptors.CalcMolFormula(g) == "C5H5N5O",
        "8-oxoguanine formula": Chem.rdMolDescriptors.CalcMolFormula(o) == "C5H5N5O2",
        "one extra oxygen": (sum(a.GetSymbol() == "O" for a in o.GetAtoms())
                             - sum(a.GetSymbol() == "O" for a in g.GetAtoms()) == 1),
        "C8 carbonyl in 8-oxoguanine": (po["O8"] is not None and
            o.GetBondBetweenAtoms(po["C8"], po["O8"]).GetBondTypeAsDouble() == 2.0),
        "guanine C8 unsubstituted": pg["O8"] is None,
        "8-oxoguanine N7 protonated": o.GetAtomWithIdx(po["N7"]).GetTotalNumHs() == 1,
        "8-oxoguanine N9 available for glycosidic bond":
            o.GetAtomWithIdx(po["N9"]).GetTotalNumHs() == 1,
        "O6 retained in both": pg["O6"] is not None and po["O6"] is not None,
        "exocyclic N2 retained in both": pg["N2"] is not None and po["N2"] is not None,
    }
    return mols, {"guanine": pg, "8-oxoguanine": po}, checks


def _coords(mol):
    AllChem.Compute2DCoords(mol)
    c = mol.GetConformer()
    return np.array([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y]
                     for i in range(mol.GetNumAtoms())])


def _place(xy, top, bottom, body_side, edge_x):
    """Stand the pairing edge vertical with `top` above `bottom`, body on `body_side`.

    The pairing edge is placed at x = edge_x so the two partners of a pair face
    each other across a gap. Mirroring in x is permitted because these are planar
    depictions with no stereocentres.
    """
    xy = xy - xy[[top, bottom]].mean(axis=0)
    v = xy[bottom] - xy[top]
    ang = -np.pi / 2 - np.arctan2(v[1], v[0])
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    xy = xy @ R.T
    centre_x = xy.mean(axis=0)[0]
    if np.sign(centre_x) != np.sign(body_side):
        xy[:, 0] = -xy[:, 0]
    edge = xy[[top, bottom]].mean(axis=0)[0]
    xy[:, 0] += edge_x - edge
    return xy


def build_pair(template, partner, t_edge, p_edge, gap=5.2):
    """Return combined molecule, coordinates and per-fragment atom offsets."""
    tm = Chem.MolFromSmiles(SMILES[template])
    pm = Chem.MolFromSmiles(SMILES[partner])
    txy, pxy = _coords(tm), _coords(pm)
    tmap = purine_atoms(tm)
    pmap = adenine_atoms(pm) if partner == "adenine" else pyrimidine_atoms(pm, partner)
    txy = _place(txy, tmap[t_edge[0]], tmap[t_edge[-1]], -1.0, -gap / 2)
    pxy = _place(pxy, pmap[p_edge[0]], pmap[p_edge[-1]], +1.0, +gap / 2)
    combo = Chem.CombineMols(tm, pm)
    off = tm.GetNumAtoms()
    conf = Chem.Conformer(combo.GetNumAtoms())
    for i, (x, y) in enumerate(np.vstack([txy, pxy])):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))
    combo.RemoveAllConformers()
    combo.AddConformer(conf, assignId=True)
    return combo, tmap, {k: v + off for k, v in pmap.items()}
