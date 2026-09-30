"""Figure 3: gene context of C>A-accessible disease-associated amino-acid states."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib import gridspec
from figure_style import (apply_style, open_frame, panel_letter, BLUE, CORAL, INK, MUTED,
                          GRAY, LGRAY, FAINT, FS_TINY_, FS_MIN, FS_SMALL)
FS_TINY=FS_TINY_
FS_BODY3C=11.0

def build_F3(G):
    apply_style(); MM=1/25.4
    CC,DC,EG,PW=G['CC'],G['DC'],G['EG'],G['PW']
    IT=G['INT']; META=G['PWMETA']
    fig=plt.figure(figsize=(168*MM,194*MM))
    Gs=gridspec.GridSpec(3,1,figure=fig,hspace=0.40,
                         height_ratios=[1.00,0.92,2.05],left=0.320,right=0.965,
                         top=0.958,bottom=0.052)

    # ---- a: expression classes against permutation nulls
    aA=fig.add_subplot(Gs[0,0]); open_frame(aA)
    _pa=aA.get_position(); aA.set_position([_pa.x0,_pa.y0,0.62*(_pa.x1-_pa.x0),_pa.height])
    y=np.arange(len(CC))[::-1]
    for yy,d in zip(y,CC):
        aA.plot([d['lo'],d['hi']],[yy,yy],lw=4.0,color=FAINT,solid_capstyle='butt',zorder=1)
        aA.plot([d['null'],d['null']],[yy-0.26,yy+0.26],lw=1.0,color=GRAY,zorder=2)
        aA.scatter([d['obs']],[yy],s=26,color=BLUE,linewidths=0,zorder=3)
        if d['hot']:
            aA.text(d['obs']+0.24,yy,'\u2020',fontsize=FS_MIN,color=INK,ha='left',va='center')
    aA.set_yticks(y); aA.set_yticklabels([d['label'] for d in CC],fontsize=FS_MIN)
    aA.tick_params(axis='y',length=0); aA.spines['left'].set_visible(False)
    aA.set_xlim(3.2,8.9); aA.set_ylim(-0.7,len(CC)+0.45); aA.set_xticks([4,5,6,7,8])
    aA.set_xlabel('Accessible fraction (%)',fontsize=FS_SMALL)
    aA.plot([],[],lw=4.0,color=FAINT,label='95% permutation interval')
    aA.plot([],[],lw=1.0,color=GRAY,label='permutation mean')
    aA.legend(loc='upper left',bbox_to_anchor=(-0.015,-0.235),ncol=2,frameon=False,
              fontsize=FS_TINY,handlelength=1.0,handletextpad=0.45,columnspacing=1.0,
              labelspacing=0.24)
    panel_letter(aA,'a')


    # ---- c: interaction between the two gene-context factors
    aC=fig.add_subplot(Gs[1,0]); open_frame(aC)
    byk={('yes','yes'):EG[0],('yes','no'):EG[1],('no','yes'):EG[2],('no','no'):EG[3]}
    xs={'no':0.0,'yes':1.0}
    for pm,col,mk in [('yes',CORAL,'o'),('no',BLUE,'s')]:
        xv=[xs[nu] for nu in ['no','yes']]
        yv=[byk[(pm,nu)]['pct'] for nu in ['no','yes']]
        aC.plot(xv,yv,lw=1.5,color=col,zorder=2,solid_capstyle='round')
        for nu in ['no','yes']:
            d=byk[(pm,nu)]
            aC.plot([xs[nu],xs[nu]],[d['lo'],d['hi']],lw=1.1,color=col,zorder=3,
                    solid_capstyle='round')
            aC.scatter([xs[nu]],[d['pct']],s=44,color=col,marker=mk,linewidths=0,zorder=4)
            aC.text(xs[nu]+(0.040 if nu=='yes' else -0.040),
                    d['pct']+(-0.28 if pm=='yes' else 0.28),f"{d['pct']:.2f}",
                    fontsize=FS_TINY,color=MUTED,va=('top' if pm=='yes' else 'bottom'),
                    ha=('left' if nu=='yes' else 'right'),zorder=5)
        aC.plot([],[],lw=1.5,color=col,marker=mk,markersize=5,
                label=('post-mitotic expression' if pm=='yes' else 'other expression classes'))
    aC.set_xticks([0,1]); aC.set_xticklabels(['no','yes'],fontsize=FS_MIN)
    aC.set_xlim(-0.22,1.98); aC.set_ylim(3.6,9.6)
    aC.set_xlabel('Neurological disease annotation',fontsize=FS_SMALL)
    aC.set_ylabel('Accessible fraction of\neligible states (%)',fontsize=FS_SMALL,linespacing=1.4)
    aC.text(0.995,0.985,f"post-mitotic \u00d7 neurological\nOR {IT['orv']:.2f}\n"
            f"95% CI {IT['lo']:.2f} to {IT['hi']:.2f}\nP = {IT['p']:.3f}",
            transform=aC.transAxes,fontsize=FS_TINY,color=INK,ha='right',va='top',
            linespacing=1.5)
    aC.legend(loc='upper left',bbox_to_anchor=(-0.015,0.99),frameon=False,
              fontsize=FS_TINY,handlelength=1.5,handletextpad=0.5,labelspacing=0.28)
    panel_letter(aC,'b',dx=-0.140,dy=1.10)

    # ---- d: pathway-level accessibility
    aD=fig.add_subplot(Gs[2,0]); open_frame(aD)
    _p=aD.get_position(); aD.set_position([0.605,_p.y0,0.895-0.605,_p.height])
    aD.set_xscale('log'); aD.axvline(1.0,lw=0.7,color=LGRAY,zorder=1)
    aD.text(1.0,1.012,'null = 1',transform=aD.get_xaxis_transform(),fontsize=FS_TINY,
            color=MUTED,ha='center',va='bottom')
    fams=[]
    for d in PW:
        if d['family'] not in fams: fams.append(d['family'])
    FAMCOL={'synaptic and neuronal signaling':CORAL,'other enriched processes':BLUE,
            'depleted processes':GRAY}
    FAMSHORT={'synaptic and neuronal signaling':'synaptic and neuronal',
              'other enriched processes':'other enriched','depleted processes':'depleted'}
    yy=[]; lab=[]; cur=0.0
    for f in fams:
        for d in [x for x in PW if x['family']==f]:
            cur-=1.58
            yy.append(cur); lab.append(d)
        cur-=0.9
    for y_,d in zip(yy,lab):
        col=FAMCOL.get(d['family'],BLUE)
        aD.plot([d['lo'],d['hi']],[y_,y_],lw=1.2,color=col,solid_capstyle='round',zorder=2)
        aD.scatter([d['or']],[y_],s=26,color=col,linewidths=0,zorder=3)
        aD.text(1.008,y_,f"{d['n']}",transform=aD.get_yaxis_transform(),fontsize=FS_TINY,
                color=MUTED,va='center',ha='left')
    for f in fams:
        aD.scatter([],[],s=26,color=FAMCOL.get(f,BLUE),label=FAMSHORT.get(f,f))
    aD.legend(loc='upper center',bbox_to_anchor=(0.45,-0.115),ncol=3,frameon=False,
              fontsize=FS_TINY,handlelength=0.8,handletextpad=0.35,columnspacing=1.1,
              scatterpoints=1)
    aD.set_yticks(yy)
    aD.set_yticklabels([d['label'] for d in lab],fontsize=FS_TINY)
    aD.tick_params(axis='y',length=0); aD.spines['left'].set_visible(False)
    aD.set_ylim(min(yy)-1.3,max(yy)+1.9)
    aD.set_xlim(0.020,4.6); aD.set_xticks([0.05,0.2,1.0,2.0])
    aD.set_xticklabels(['0.05','0.2','1.0','2.0'])
    aD.set_xlabel('Odds ratio of accessibility (log scale)',fontsize=FS_SMALL)
    aD.text(1.008,max(yy)+1.45,'genes',transform=aD.get_yaxis_transform(),
            fontsize=FS_TINY,color=MUTED,va='center',ha='left')
    panel_letter(aD,'c')
    return fig,(aA,aC,aD)
