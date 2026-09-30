"""Figure 4: C>A-accessible amino-acid states observed in human proteomic data."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib import gridspec
from figure_style import (apply_style, open_frame, panel_letter, BLUE, CORAL, INK, MUTED,
                          GRAY, LGRAY, FAINT, FS_TINY_, FS_MIN, FS_SMALL)
FS_TINY=FS_TINY_

def build_F4(G):
    apply_style(); MM=1/25.4
    FUN,SUB,GEN=G['FUN'],G['SUB'],G['GEN']
    fig=plt.figure(figsize=(160*MM,124*MM))
    Gs=gridspec.GridSpec(2,2,figure=fig,hspace=0.34,wspace=0.90,height_ratios=[0.50,1.00],
                         left=0.295,right=0.965,top=0.955,bottom=0.085)

    # ---- a: stepped filtering cascade from source sites to the accessible subset
    aA=fig.add_subplot(Gs[0,:]); aA.set_axis_off()
    n=len(FUN); bw=0.50; bh=0.40; gap=0.15; step=0.085
    for i,d in enumerate(FUN):
        last=(i==n-1)
        x0=i*step; top=-(i*(bh+gap))
        aA.add_patch(plt.Rectangle((x0,top-bh),bw,bh,transform=aA.transData,
                     facecolor=(CORAL if last else FAINT),
                     edgecolor=(CORAL if last else LGRAY),linewidth=0.8,zorder=2))
        aA.text(x0+bw/2,top-bh/2,f"{d['n']:,}",fontsize=FS_MIN,
                color=('white' if last else INK),va='center',ha='center',zorder=3,
                fontweight='bold')
        aA.text(x0-0.022,top-bh/2,d['label'],fontsize=FS_MIN,color=INK,
                va='center',ha='right',zorder=3,fontweight=('bold' if last else 'normal'))
        if not last:
            xc=(i+1)*step+0.05
            aA.annotate('',xy=(xc,top-bh-gap+0.02),xytext=(xc,top-bh-0.02),
                        xycoords='data',textcoords='data',
                        arrowprops=dict(arrowstyle='-|>',color=GRAY,lw=0.8,
                                        shrinkA=0,shrinkB=0,mutation_scale=6))
    aA.set_xlim(-1.02,(n-1)*step+bw+0.02)
    aA.set_ylim(-(n*(bh+gap))+gap-0.04,0.04)
    panel_letter(aA,'a')

    # ---- b: substitution classes among the observed accessible states
    aB=fig.add_subplot(Gs[1,0]); open_frame(aB)
    y2=np.arange(len(SUB))[::-1]
    for yy,d in zip(y2,SUB):
        aB.hlines(yy,0,d['n'],lw=0.9,color=LGRAY,zorder=1)
        aB.scatter([d['n']],[yy],s=34,color=BLUE,linewidths=0,zorder=3)
    aB.set_yticks(y2); aB.set_yticklabels([d['label'] for d in SUB],fontsize=FS_MIN)
    aB.tick_params(axis='y',length=0); aB.spines['left'].set_visible(False)
    aB.set_xlim(0,48); aB.set_xticks([0,10,20,30,40])
    aB.set_ylim(-0.7,len(SUB)-0.3)
    aB.set_xlabel('Observed states',fontsize=FS_SMALL)
    panel_letter(aB,'b')

    # ---- c: genes carrying the observed states
    aC=fig.add_subplot(Gs[1,1]); open_frame(aC)
    y3=np.arange(len(GEN))[::-1]
    for yy,d in zip(y3,GEN):
        col=CORAL if d['hb'] else BLUE
        aC.hlines(yy,0,d['n'],lw=0.9,color=LGRAY,zorder=1)
        aC.scatter([d['n']],[yy],s=34,color=col,linewidths=0,zorder=3)
    aC.set_yticks(y3)
    aC.set_yticklabels([d['gene'] for d in GEN],fontsize=FS_MIN,fontstyle='italic')
    aC.tick_params(axis='y',length=0); aC.spines['left'].set_visible(False)
    aC.set_xlim(0,14.5); aC.set_xticks([0,5,10])
    aC.set_ylim(-0.7,len(GEN)-0.3)
    aC.set_xlabel('Observed states',fontsize=FS_SMALL)
    aC.scatter([],[],s=34,color=CORAL,label='hemoglobin cluster')
    aC.scatter([],[],s=34,color=BLUE,label='other genes')
    aC.legend(loc='lower right',bbox_to_anchor=(1.04,-0.02),frameon=False,fontsize=FS_TINY,
              handlelength=0.9,handletextpad=0.4,labelspacing=0.28)
    panel_letter(aC,'c')
    return fig,(aA,aB,aC)
