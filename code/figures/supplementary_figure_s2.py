"""Supplementary Figure S2: oxidative-lesion evidence at accessible coordinates."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib import gridspec
from figure_style import (apply_style, open_frame, panel_letter, BLUE, CORAL, INK, MUTED,
                          GRAY, LGRAY, FAINT, FS_TINY_, FS_MIN, FS_SMALL)
FS_TINY=FS_TINY_

SHORT={'endogenous conditions': 'endogenous',
       'induced conditions only': 'induced only',
       'no detection': 'not detected',
       'lesion coordinate unresolved': 'unresolved'}

def build_S2(G):
    apply_style(); MM=1/25.4
    BARS,CATS,ST=G['BARS'],G['CATS'],G['STAT']
    fig=plt.figure(figsize=(150*MM,66*MM))
    Gs=gridspec.GridSpec(1,2,figure=fig,wspace=0.62,width_ratios=[0.88,1.00],
                         left=0.215,right=0.975,top=0.86,bottom=0.235)

    aA=fig.add_subplot(Gs[0,0]); open_frame(aA)
    y=np.arange(len(BARS))[::-1]
    for yy,d in zip(y,BARS):
        col=CORAL if d['hot'] else BLUE
        aA.plot([d['lo'],d['hi']],[yy,yy],lw=1.3,color=col,solid_capstyle='round',zorder=2)
        aA.scatter([d['pct']],[yy],s=40,color=col,linewidths=0,zorder=3)
        aA.text(d['pct'],yy+0.30,f"{d['k']:,} of {d['n']:,}",fontsize=FS_TINY,color=MUTED,
                ha='center',va='bottom')
    aA.set_yticks(y); aA.set_yticklabels([d['label'] for d in BARS],fontsize=FS_MIN,
                                         linespacing=1.35)
    aA.tick_params(axis='y',length=0); aA.spines['left'].set_visible(False)
    aA.set_xlim(53.1,56.5); aA.set_ylim(-1.55,len(BARS)-0.16)
    aA.set_xlabel('Coordinates with oxidized-guanine detection (%)',fontsize=FS_SMALL)
    aA.text(0.99,0.015,f"OR {ST['orv']:.2f} ({ST['lo']:.2f} to {ST['hi']:.2f})\n"
            f"P = {ST['p']:.2f}",transform=aA.transAxes,fontsize=FS_TINY,color=MUTED,
            ha='right',va='bottom',linespacing=1.5)
    panel_letter(aA,'a')

    aB=fig.add_subplot(Gs[0,1]); aB.set_axis_off()
    aB.set_xlim(0,100); aB.set_ylim(0,100)
    tot=G['NTOT']; x=0.0; ybar=52.0; hbar=17.0
    cols=[CORAL,'#E8A598',FAINT,LGRAY]
    for d,col in zip(CATS,cols):
        w=100.0*d['n']/tot
        aB.add_patch(plt.Rectangle((x,ybar),w,hbar,facecolor=col,edgecolor='white',lw=0.8,zorder=2))
        if w>=6.0:
            aB.text(x+w/2,ybar+hbar/2,f"{d['n']}",fontsize=FS_MIN,
                    color=('white' if col==CORAL else INK),ha='center',va='center',zorder=3)
        else:
            aB.plot([x+w/2,x+w/2],[ybar+hbar+0.8,ybar+hbar+5.2],lw=0.6,color='#D6D9DC',zorder=1)
            aB.text(x+w/2,ybar+hbar+6.0,f"{d['n']}",fontsize=FS_MIN,color=INK,
                    ha='center',va='bottom',zorder=3)
        x+=w
    xl=0.0
    for i,(d,col) in enumerate(zip(CATS,cols)):
        w=100.0*d['n']/tot
        yy=ybar-7.0-(i%2)*13.0
        aB.plot([xl+w/2,xl+w/2],[ybar-1.0,yy+4.0],lw=0.6,color='#D6D9DC',zorder=1)
        aB.text(xl+w/2,yy,SHORT[d['label']],fontsize=FS_TINY,color=MUTED,
                ha='center',va='top',linespacing=1.3,zorder=3)
        xl+=w
    aB.text(0,ybar+hbar+12.5,f'n = {tot} observed accessible states',fontsize=FS_TINY,
            color=MUTED,ha='left',va='bottom')
    panel_letter(aB,'b',dx=-0.045,dy=0.99)
    return fig,(aA,aB)
