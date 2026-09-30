"""Manuscript visual system. Journal-width scientific figures, no decorative fills.

Semantic colors (do not extend):
  BLUE   C>A-accessible states and general data
  PURPLE regulatory / phosphorylation outcome
  CORAL  altered nucleotide, stop, genomic pathogenic allele
  GRAY   reference and background
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtr

BLUE='#3C6E9F'; PURPLE='#6B4E8F'; CORAL='#B5462F'
INK='#1A1A1A'; MUTED='#5E6266'; GRAY='#9AA0A6'
LGRAY='#C9CDD1'; FAINT='#EFF1F2'; WHISPER='#F7F8F9'

# journal width 180 mm = 7.09 in.
# Type hierarchy at final printed size: panel letters 10 pt (the only bold text),
# axis labels 9 pt, tick labels 8 pt, annotations 7.5 pt floor, legends 8 pt.
FS_TINY_, FS_MIN, FS_SMALL, FS_BODY, FS_LABEL, FS_BIG = 8.0, 8.0, 9.0, 9.5, 10.0, 20.0
ARROW='\u2192'
# Arial throughout, including codons: every base is centred in its own tile or on a
# fixed grid, so a monospace face is not needed for column alignment.
import os

def resolve_font():
    """Return the figure font: Arial when installed, otherwise an explicit open fallback.

    Resolving once avoids re-requesting a missing family for every text object, which is
    what produces the font-warning flood and leaves the geometry environment dependent.
    Set OXOG_FIGURE_FONT to force a family, which is how the fallback path is exercised
    on machines where Arial is present.
    """
    forced = os.environ.get('OXOG_FIGURE_FONT')
    if forced:
        return forced
    from matplotlib import font_manager
    return 'Arial' if 'Arial' in {f.name for f in font_manager.fontManager.ttflist} else 'DejaVu Sans'


FONT = resolve_font()
FONT_FAMILY = FONT
MONO = FONT

def apply_style():
    mpl.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':[FONT],
        'axes.spines.top':False,'axes.spines.right':False,
        'axes.grid':False,'axes.edgecolor':INK,'axes.linewidth':0.6,
        'axes.labelcolor':INK,'text.color':INK,
        'xtick.color':INK,'ytick.color':INK,
        'xtick.labelsize':FS_MIN,'ytick.labelsize':FS_MIN,
        'xtick.major.width':0.6,'ytick.major.width':0.6,
        'xtick.major.size':2.5,'ytick.major.size':2.5,
        'axes.labelsize':FS_SMALL,'legend.fontsize':FS_MIN,
        'figure.facecolor':'white','savefig.facecolor':'white',
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        # Arial has no Unicode subscript glyphs; route mathtext through Arial so
        # chemical subscripts stay in the same face as the rest of the figure.
        'mathtext.fontset':'custom','mathtext.rm':FONT,
        'mathtext.it':FONT,'mathtext.bf':FONT,
        'mathtext.default':'regular'})

def bare(ax):
    """Strip an axes to a blank drawing surface."""
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

def open_frame(ax):
    """Left and bottom spines only."""
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(True); ax.spines['bottom'].set_visible(True)

def panel_letter(ax,letter,dx=-0.085,dy=1.03):
    """Panel letter, regular weight, one fixed offset for every panel."""
    ax.text(dx,dy,letter,transform=ax.transAxes,fontsize=FS_LABEL,
            fontweight='normal',va='bottom',ha='left',color=INK)

def nt_row(ax,x,y,seq,w,h,hi=None,hicol=None,size=FS_BODY,aa=None,
           arrow=None,aacol=None,lab=None,labcol=None,labdx=None):
    """One codon as outlined cells: black letters on white, thin colored outline on the
    altered base only. Optional subtle arrow to the encoded residue."""
    import matplotlib.patches as mp
    if arrow is None: arrow=w*1.5
    for i,ch in enumerate(seq):
        cx=x+i*w
        onhi=(hi is not None and i==hi)
        ax.add_patch(mp.Rectangle((cx-w*0.44,y-h*0.5),w*0.88,h,facecolor='white',
                     edgecolor=(hicol or CORAL) if onhi else LGRAY,
                     linewidth=1.0 if onhi else 0.5,zorder=2,clip_on=False,
                     joinstyle='miter',capstyle='projecting'))
        ax.text(cx,y,ch,fontsize=size,color=INK,ha='center',va='center',zorder=3,clip_on=False)
    xe=x+(len(seq)-1)*w+w*0.44
    if aa is not None:
        x0=xe+w*0.30; x1=x0+arrow
        ax.annotate('',xy=(x1,y),xytext=(x0,y),zorder=2,
                    arrowprops=dict(arrowstyle='-|>',lw=0.6,color=GRAY,
                                    shrinkA=0,shrinkB=0,mutation_scale=6))
        ax.text(x1+w*0.30,y,aa,fontsize=size,color=aacol or INK,ha='left',va='center',zorder=3)
    if lab is not None:
        ax.text(x-w*0.75 if labdx is None else x+labdx,y,lab,fontsize=FS_TINY_,
                color=labcol or MUTED,ha='right',va='center',zorder=3)
    return xe

def codon(ax,x,y,seq,hi=None,size=FS_BODY,ha='left',hicol=CORAL):
    """Monospace codon; color only the highlighted index."""
    t=ax.transData if ha=='left' else ax.transData
    tx=ax.text(x,y,'',transform=t)
    for i,ch in enumerate(seq):
        c=hicol if (hi is not None and i==hi) else INK
        w=ax.text(x,y,ch,transform=tx.get_transform(),family=MONO,fontsize=size,
                  color=c,fontweight='normal',va='center')
        tx=w
        tx._transform=mtr.offset_copy(w.get_transform(),x=size*0.66,y=0,
                                      units='points',fig=ax.figure)
    return tx

def seqrow(ax,x,y,seq,hi=None,size=FS_BODY,step=0.42,hicol=CORAL):
    """Draw a sequence on a fixed grid; returns x of each character."""
    xs=[]
    for i,ch in enumerate(seq):
        xx=x+i*step
        c=hicol if (hi is not None and i==hi) else INK
        ax.text(xx,y,ch,family=MONO,fontsize=size,color=c,ha='center',va='center',
                fontweight='normal')
        xs.append(xx)
    return xs

def check(fig,axes):
    """Return (overlaps, undersized) for a render-time QC gate."""
    r=fig.canvas.get_renderer()
    tk={id(t) for ax in axes for t in ax.get_xticklabels()+ax.get_yticklabels()}
    ts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
        if t.get_text().strip() and t.get_visible() and id(t) not in tk
        and len(t.get_text())>2]
    ov=[(a.get_text()[:20],b.get_text()[:20]) for i,(a,ba) in enumerate(ts)
        for b,bb in ts[i+1:] if ba.overlaps(bb)]
    small=[t.get_text()[:24] for t in fig.findobj(mpl.text.Text)
           if t.get_text().strip() and t.get_visible() and t.get_fontsize()<FS_MIN-0.01]
    return ov,small

FS_TINY = FS_TINY_  # annotation floor (7.5 pt) used by the figure builders

def no_overlap(fig,axes):
    return check(fig,axes)[0]

def type_floor(fig,floor=None):
    f=FS_MIN-0.01 if floor is None else floor
    return [t.get_text()[:24] for t in fig.findobj(mpl.text.Text)
            if t.get_text().strip() and t.get_visible() and t.get_fontsize()<f]


FONT_NOTE = ("Figures use Arial when it is installed and DejaVu Sans otherwise, resolved once "
             "at import. No font files are redistributed with this archive.")

def save_figure(fig, axes, stem, formats=('png', 'pdf', 'svg')):
    """Write a finished figure to the publication formats used in the manuscript."""
    for ext in formats:
        fig.savefig(f'{stem}.{ext}', dpi=600 if ext == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    return [f'{stem}.{ext}' for ext in formats]
