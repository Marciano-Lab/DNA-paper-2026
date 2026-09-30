"""Figure 1: the genetic code confines 8-oxoguanine transcriptional mutagenesis.

Panel a draws 2'-deoxyguanosine and 8-oxo-7,8-dihydro-2'-deoxyguanosine from SMILES carried in
the figure input, so no external drawing is reproduced and no network access is required.
"""
import matplotlib.lines as mlines
import io
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as mpatches
from figure_style import (apply_style, open_frame, panel_letter, nt_row, LGRAY,
                       FS_TINY, FS_MIN, FS_SMALL, FS_BODY, MONO, INK, LGRAY, GRAY)

CA='#6B4E9B'; OTH='#DDE1E4'; NON='#F8F8F9'; ACC='#B03A2E'; PUR='#6B4E9B'; MUT='#6B7178'


def _pair(template, partner, t_edge, p_edge):
    """Build a chemically verified template:incoming base pair (see figure_1a_chemistry)."""
    import figure_1a_chemistry as chem
    _, _, checks = chem.verify()
    assert all(checks.values()), [k for k, v in checks.items() if not v]
    return chem.build_pair(template, partner, t_edge, p_edge, gap=3.4)


def _pair_image(combo, px=(900, 520), bond=30.0):
    """Render a base pair and return the image plus an atom-to-pixel accessor."""
    import io
    from PIL import Image as _PIL
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Geometry import Point2D
    d = rdMolDraw2D.MolDraw2DCairo(*px)
    o = d.drawOptions()
    o.fixedBondLength = bond
    o.bondLineWidth = 2
    o.minFontSize = 16
    o.maxFontSize = 16
    o.clearBackground = False
    o.padding = 0.04
    rdMolDraw2D.PrepareAndDrawMolecule(d, combo)
    d.FinishDrawing()
    img = np.array(_PIL.open(io.BytesIO(d.GetDrawingText())).convert("RGBA"))
    conf = combo.GetConformer()

    def anchor(idx):
        p = conf.GetAtomPosition(idx)
        q = d.GetDrawCoords(Point2D(p.x, p.y))
        return np.array([q.x, q.y])

    return img, anchor


def _structures(smiles_by_key, px=1500):
    """Render the nucleosides from SMILES with a shared scaffold orientation.

    The oxidized nucleoside is depicted against the reference coordinates of the unmodified
    one, so the two panels differ only where the chemistry differs.
    """
    from PIL import Image
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    mols = {k: Chem.MolFromSmiles(v) for k, v in smiles_by_key.items()}
    ref = mols['dG']
    rdDepictor.Compute2DCoords(ref)
    core = Chem.MolFromSmarts(rdFMCS.FindMCS(list(mols.values()), timeout=20,
                                             ringMatchesRingOnly=True).smartsString)
    for key, mol in mols.items():
        if key == 'dG':
            continue
        rdDepictor.Compute2DCoords(mol)
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, ref, refPatt=core, acceptFailure=True)
        except Exception:
            pass
    out = {}
    for key, mol in mols.items():
        drawer = rdMolDraw2D.MolDraw2DCairo(px, int(px*0.80))
        opts = drawer.drawOptions()
        opts.bondLineWidth = 5
        opts.minFontSize = 46
        opts.maxFontSize = 46
        opts.padding = 0.06
        opts.clearBackground = True
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert('RGB')
        arr = np.asarray(img)
        ink = np.where((arr < 245).any(axis=2))
        if len(ink[0]):
            r0, r1 = ink[0].min(), ink[0].max(); c0, c1 = ink[1].min(), ink[1].max()
            img = img.crop((max(0, c0-8), max(0, r0-8),
                            min(arr.shape[1], c1+9), min(arr.shape[0], r1+9)))
        out[key] = img
    return out


def _databox(fig, ax, x0, x1, y0, y1):
    """Convert a data-coordinate box on ax into a figure-coordinate axes rectangle."""
    p0 = fig.transFigure.inverted().transform(ax.transData.transform((x0, y0)))
    p1 = fig.transFigure.inverted().transform(ax.transData.transform((x1, y1)))
    return [p0[0], p0[1], p1[0]-p0[0], p1[1]-p0[1]]


def build_F1(G):
    AA20,PAIRS,ANYP,REPL,INTRO,DIRN=G['AA20'],G['PAIRS'],G['ANYP'],G['REPL'],G['INTRO'],G['DIRN']
    CHEM=G['CHEM']
    apply_style(); MM=1/25.4
    fig=plt.figure(figsize=(174*MM,224*MM))
    Gs=gridspec.GridSpec(3,1,figure=fig,hspace=0.30,height_ratios=[1.16,1.52,0.44])

    # ---- a: template base pairing and the resulting transcript base
    aA=fig.add_subplot(Gs[0]); aA.set_axis_off()
    aA.set_xlim(0,100); aA.set_ylim(0,100)
    XP0,XP1=6,62; XR,XC=76,91
    aA.text((XP0+XP1)/2,97.5,'template base and incoming RNA base',ha='center',va='center',
            fontsize=FS_TINY,color=MUT)
    aA.text(XR,97.5,'transcript',ha='center',va='center',fontsize=FS_TINY,color=MUT)
    aA.text(XC,97.5,'coding sense',ha='center',va='center',fontsize=FS_TINY,color=MUT)
    aA.plot([4,95],[93.5,93.5],lw=0.7,color='#C9CDD2')
    PAIRSPEC=[(50,90,'guanine','cytosine',['O6','N1','N2'],['N4','N3','O2'],
               'Watson-Crick','G(anti):C','C','C',INK),
              (6,46,'8-oxoguanine','adenine',['O6','N7'],['N6','N1'],
               'Hoogsteen-like','8-oxoG(syn):A','A','A',PUR)]
    for ylo,yhi,tname,pname,tedge,pedge,geom,pairlab,rb,cb,col in PAIRSPEC:
        combo,tmap,pmap=_pair(tname,pname,tedge,pedge)
        img,anchor=_pair_image(combo)
        h=(yhi-ylo); w=h*img.shape[1]/img.shape[0]*(224/174)*(100/100)
        w=min(w,XP1-XP0)
        x0=XP0+((XP1-XP0)-w)/2
        aA.imshow(img,extent=(x0,x0+w,ylo,yhi),aspect='auto',zorder=2,interpolation='lanczos')
        to_data=lambda p:(x0+p[0]/img.shape[1]*w, yhi-p[1]/img.shape[0]*h)
        for a,b in zip([tmap[k] for k in tedge],[pmap[k] for k in pedge]):
            pa=np.array(to_data(anchor(a))); pb=np.array(to_data(anchor(b)))
            u=(pb-pa)/np.linalg.norm(pb-pa)
            aA.plot(*zip(pa+u*1.05,pb-u*1.05),ls=(0,(1,1.5)),lw=1.5,color=col,zorder=4)
        p9=np.array(to_data(anchor(tmap['N9'])))
        aA.annotate('deoxyribose',xy=(p9[0],p9[1]),xytext=(p9[0]-3.2,p9[1]-3.4),
                    fontsize=FS_TINY,color=MUT,ha='right',va='center',
                    arrowprops=dict(arrowstyle='-',lw=0.7,color='#B8BDC2'),zorder=5)
        yc=(ylo+yhi)/2
        aA.text(x0+w/2,yhi-1.0,f'{geom}   {pairlab}',ha='center',va='top',fontsize=FS_MIN,
                color=col,fontweight=('bold' if col!=INK else 'normal'),zorder=5)
        aA.annotate('',xy=(XR-4.5,yc),xytext=(XP1+1.5,yc),
                    arrowprops=dict(arrowstyle='-|>',lw=0.9,color='#B8BDC2',mutation_scale=8))
        aA.text(XR,yc,rb,ha='center',va='center',fontsize=13,family=MONO,color=col)
        aA.text(XC,yc,cb,ha='center',va='center',fontsize=14,family=MONO,color=col,
                fontweight=('bold' if col!=INK else 'normal'))
    ymid=(PAIRSPEC[0][0]+PAIRSPEC[1][1])/2
    aA.text((XP1+1.5+XR-4.5)/2,ymid,'RNA Pol II',ha='center',va='center',fontsize=FS_TINY,color=MUT)
    aA.plot([4,95],[3.5,3.5],lw=0.7,color='#C9CDD2')
    aA.text(XC,0.5,'C>A',ha='center',va='center',fontsize=FS_BODY,color=PUR,fontweight='bold')
    panel_letter(aA,'a')

    # ---- b: 20x20 amino-acid state matrix
    aB=fig.add_subplot(Gs[1]); open_frame(aB)
    n=len(AA20); grid=np.zeros((n,n))
    for i,a in enumerate(AA20):
        for j,b in enumerate(AA20):
            if a==b: grid[i,j]=np.nan
            elif (a,b) in PAIRS: grid[i,j]=2
            elif (a,b) in ANYP: grid[i,j]=1
    cm=mpl.colors.ListedColormap([NON,OTH,CA])
    aB.imshow(np.nan_to_num(grid),cmap=cm,vmin=0,vmax=2,interpolation='nearest')
    for i in range(n):
        aB.add_patch(mpatches.Rectangle((i-0.5,i-0.5),1,1,facecolor='white',edgecolor='none'))
    # outline the C>A cells so they dominate
    for i,a in enumerate(AA20):
        for j,b in enumerate(AA20):
            if (a,b) in PAIRS:
                aB.add_patch(mpatches.Rectangle((j-0.5,i-0.5),1,1,facecolor='none',
                                                edgecolor=INK,lw=0.9,zorder=4))
    aB.set_xticks(range(n)); aB.set_xticklabels(AA20,fontsize=FS_BODY+0.5,family=MONO)
    aB.set_yticks(range(n)); aB.set_yticklabels(AA20,fontsize=FS_BODY+0.5,family=MONO)
    aB.set_xlabel('Product residue',fontsize=FS_SMALL)
    aB.set_ylabel('Reference residue',fontsize=FS_SMALL)
    aB.tick_params(length=2)
    for s in aB.spines.values(): s.set_visible(False)
    hs=[mpatches.Patch(facecolor=CA,edgecolor=INK,lw=0.9,label=f'C>A ({len(PAIRS)})'),
        mpatches.Patch(facecolor=OTH,edgecolor='none',
                       label=f'other single base ({len(ANYP)-len(PAIRS)})'),
        mpatches.Patch(facecolor=NON,edgecolor='none',
                       label=f'no single-base route ({n*(n-1)-len(ANYP)})')]
    aB.legend(handles=hs,loc='lower left',bbox_to_anchor=(0.0,1.012),ncol=3,frameon=False,
              fontsize=FS_MIN,handlelength=1.0,columnspacing=1.4,handletextpad=0.4)
    panel_letter(aB,'b')

    # ---- c: directional accessibility over all twenty residues
    aC=fig.add_subplot(Gs[2]); aC.set_axis_off()
    xr=np.arange(len(AA20))
    ROWS=[('Residues that can be replaced',[a in REPL for a in AA20],CA,1.0),
          ('Residues that can be introduced',[a in INTRO for a in AA20],PUR,0.0)]
    only={a for a in AA20 if (a in REPL)!=(a in INTRO)}
    for lab,vals,col,yy in ROWS:
        aC.text(-0.9,yy,lab,fontsize=FS_MIN,color=INK,ha='right',va='center')
        for xi,(a,v) in enumerate(zip(AA20,vals)):
            if v:
                aC.scatter([xi],[yy],s=52,marker='s',color=col,linewidths=0,zorder=3)
            else:
                aC.scatter([xi],[yy],s=52,marker='s',facecolor='white',edgecolor=LGRAY,
                           linewidths=0.7,zorder=3)
    for xi,a in enumerate(AA20):
        aC.text(xi,-0.62,a,fontsize=FS_BODY,family=MONO,ha='center',va='center',
                color=INK if a in only else MUT)
        if a in only:
            aC.scatter([xi],[-0.62],s=150,marker='s',facecolor='none',edgecolor=LGRAY,
                       linewidths=0.6,zorder=1)
    aC.text(-0.9,-0.62,'One-way residues outlined',fontsize=FS_TINY,color=MUT,ha='right',va='center')
    aC.set_xlim(-9.2,len(AA20)-0.4); aC.set_ylim(-1.22,1.24)
    panel_letter(aC,'c',dx=-0.055,dy=0.88)

    return fig,(aA,aB,aC)
