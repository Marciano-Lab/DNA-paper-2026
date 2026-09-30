"""Supplementary Figure S1: specification sensitivity of the gene-context interaction."""
import numpy as np, matplotlib.pyplot as plt
from figure_style import (apply_style, open_frame, BLUE, CORAL, MUTED, LGRAY, FS_TINY_, FS_MIN, FS_SMALL)
FS_TINY=FS_TINY_

def build_S1(G):
    apply_style(); MM=1/25.4
    R=G['ROWS']
    fig=plt.figure(figsize=(150*MM,58*MM))
    ax=fig.add_subplot(111); open_frame(ax)
    fig.subplots_adjust(left=0.315,right=0.70,top=0.90,bottom=0.235)
    y=np.arange(len(R))[::-1]
    ax.axvline(1.0,lw=0.7,color=LGRAY,zorder=1)
    for yy,d in zip(y,R):
        c=CORAL if d['primary'] else BLUE
        ax.plot([d['lo'],d['hi']],[yy,yy],lw=1.3 if d['primary'] else 1.0,color=c,
                solid_capstyle='round',zorder=2,alpha=1.0 if d['primary'] else 0.8)
        ax.scatter([d['or']],[yy],s=38 if d['primary'] else 24,color=c,linewidths=0,zorder=3)
        ax.text(1.03,yy,f"{d['or']:.2f} ({d['lo']:.2f}, {d['hi']:.2f})",
                transform=ax.get_yaxis_transform(),fontsize=FS_TINY,
                color=c if d['primary'] else MUTED,va='center',ha='left')
    ax.set_yticks(y); ax.set_yticklabels([d['label'] for d in R],fontsize=FS_MIN)
    ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.set_xlim(0.85,3.95); ax.set_ylim(-0.65,len(R)-0.15)
    ax.set_xticks([1,2,3])
    ax.set_xlabel('Post-mitotic by neurological interaction odds ratio',fontsize=FS_SMALL)
    ax.text(1.03,len(R)-0.14,'odds ratio (95% CI)',transform=ax.get_yaxis_transform(),
            fontsize=FS_TINY,color=MUTED,va='bottom',ha='left')
    return fig,(ax,)
