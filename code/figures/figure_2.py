import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import gridspec
from figure_style import (apply_style, open_frame, panel_letter, nt_row,
                       FS_TINY, FS_MIN, FS_SMALL, FS_BODY, MONO, INK, LGRAY, GRAY, MUTED)

BLU='#3E6D9C'; ACC='#B03A2E'; PUR='#6B4E9B'; GRY='#9AA0A6'; MUT='#6B7178'; PALE='#DCDEE1'

def build_F2(G):
    SUB, CVC, GJ, FAM = G['SUB'], G['CVC'], G['GJ'], G['FAM']
    apply_style(); MM=1/25.4
    fig=plt.figure(figsize=(150*MM,168*MM))
    Gs=gridspec.GridSpec(2,2,figure=fig,hspace=0.50,wspace=0.80,
                         height_ratios=[1.75,0.92],
                         width_ratios=[0.80,1.20])

    # ---- a: which substitutions carry the disease-associated load
    aA=fig.add_subplot(Gs[0,0]); open_frame(aA)
    ks=[k for k,_ in SUB][::-1]
    allv=[v[0] for _,v in SUB][::-1]; consv=[v[1] for _,v in SUB][::-1]
    y=np.arange(len(ks))
    aA.hlines(y,consv,allv,lw=1.0,color=PALE,zorder=1)
    aA.scatter(allv,y,s=26,facecolor='white',edgecolor=BLU,linewidths=1.0,zorder=3,
               label='all accessible')
    aA.scatter(consv,y,s=26,color=BLU,linewidths=0,zorder=3,label='higher-confidence subset')
    aA.set_yticks(y); aA.set_yticklabels(ks,fontsize=FS_TINY,family=MONO)
    aA.set_xlabel('Accessible disease-associated states',fontsize=FS_SMALL)
    aA.legend(handles=[
        mpl.lines.Line2D([],[],marker='o',color='none',markerfacecolor='white',
                         markeredgecolor=BLU,markeredgewidth=1.0,markersize=5.4,
                         label='all accessible states'),
        mpl.lines.Line2D([],[],marker='o',color='none',markerfacecolor=BLU,
                         markeredgecolor=BLU,markersize=5.4,label='higher-confidence subset')],
        loc='upper center',bbox_to_anchor=(0.5,-0.155),ncol=2,frameon=False,
        fontsize=FS_TINY,handlelength=0.9,handletextpad=0.4,columnspacing=1.6)
    aA.set_xlim(0,max(allv)*1.14); aA.set_ylim(-0.8,len(ks)-0.2)
    aA.tick_params(axis='y',length=0); aA.tick_params(axis='x',length=2)
    panel_letter(aA,'a')

    # ---- b: one worked different-position state, nucleotide by nucleotide
    aB=fig.add_subplot(Gs[0,1]); aB.set_axis_off()
    aB.set_xticks([]); aB.set_yticks([])
    aB.set_xlim(0,100); aB.set_ylim(0,100)
    aB.text(0,95,f"{GJ['gene']} {GJ['wt']}{GJ['pos']}{GJ['mut']}",fontsize=FS_SMALL,color=INK)
    rows=[('reference',GJ['ref'],GJ['ref_aa'],None,None),
          ('genomic allele',GJ['gen'],GJ['aa'],GJ['gen_pos'],ACC),
          ('C>A transcript',GJ['ca'],GJ['aa'],GJ['ca_pos'],PUR)]
    for (lab,cod,aa,hi,col),yy in zip(rows,[62,40,18]):
        nt_row(aB,44,yy,cod,w=10.5,h=11.5,hi=None if hi is None else hi-1,hicol=col,
               size=FS_SMALL,aa=aa,arrow=12.0,aacol=INK,lab=lab,labcol=MUT,labdx=-9.5)
    panel_letter(aB,'b')

    # ---- c: the three nucleotide routes, linear scale
    aC=fig.add_subplot(Gs[1,0]); open_frame(aC)
    labs=[r[0] for r in CVC][::-1]
    tot=[r[1] for r in CVC][::-1]; cns=[r[2] for r in CVC][::-1]
    yc=np.arange(len(labs))
    for yy,t,c,lab in zip(yc,tot,cns,labs):
        hot='Different codon position' in lab
        col=PUR if hot else PALE
        aC.barh(yy,t,height=0.48,color=col,edgecolor='none',zorder=2)
        aC.text(1.02,yy,f'{t:,} ({c:,})',transform=aC.get_yaxis_transform(),
                fontsize=(FS_MIN if hot else FS_TINY),color=INK,va='center',ha='left',
                fontweight=('bold' if hot else 'normal'))
    aC.text(1.02,len(labs)-0.30,'total\n(higher confidence)',transform=aC.get_yaxis_transform(),
            fontsize=FS_TINY,color=MUT,va='bottom',ha='left')
    aC.set_yticks(yc); aC.set_yticklabels(labs,fontsize=FS_MIN,linespacing=1.3)
    aC.set_xlabel('Accessible disease-associated states',fontsize=FS_SMALL)
    aC.set_xlim(0,3450); aC.set_ylim(-0.6,len(labs)+0.10)
    aC.set_xticks([0,1000,2000,3000])
    aC.tick_params(axis='y',length=0); aC.spines['left'].set_visible(False)
    panel_letter(aC,'c')

    # ---- d: the two codon configurations that permit different-position convergence
    aD=fig.add_subplot(Gs[1,1]); aD.set_axis_off()
    aD.set_xticks([]); aD.set_yticks([])
    aD.set_xlim(0,100); aD.set_ylim(0,100)
    AA3={'F':'Phe','L':'Leu','S':'Ser','R':'Arg'}
    CODE={'TTC':'F','TTA':'L','CTC':'L','AGC':'S','AGA':'R','CGC':'R'}
    for (n,sub,ref,ca,gen,cap,cgp),ytop in zip(FAM,[99,47]):
        short=sub.replace(' to ','\u2192')
        aD.text(0,ytop,f'{short} ({n})',fontsize=FS_MIN,color=INK,ha='left',va='center')
        for cod,yy,hi,col in [(ref,ytop-14,None,None),(ca,ytop-26,cap,PUR),
                              (gen,ytop-38,cgp,ACC)]:
            nt_row(aD,38,yy,cod,w=11.4,h=12.4,hi=None if hi is None else hi-1,hicol=col,
                   size=FS_SMALL,aa=AA3[CODE[cod]],arrow=9.0,aacol=INK,
                   lab={ref:'reference',ca:'C>A transcript',gen:'genomic allele'}[cod],
                   labcol=MUT,labdx=-8.0)
    panel_letter(aD,'d')

    return fig,(aA,aB,aC,aD)
