#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#====change the candels catalog to Guo 2013
#Created on Tue Jan 23 08:16:14 2018
#@author: ashley

import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np 
from math import radians,cos,sin,asin,sqrt
import math
from astropy.io import fits


def lamda2v(x):
    return 3e18/x


def obs2emit(v,z):
    return (1+z)*v


def find_all_index(arr,item):
    return[i for i,a in enumerate(arr) if a==item]
    
    
def nummatch(f,v):
    f_da=[i for i in f if i>0]
    f_fu=[i for i in f if i<=0]
    index=[]
    for i in range(0,len(f_da),1):
        index_fu=f.index(f_da[i])
        index.append(index_fu)#dayuling de zuobiao 
    v_final=[]
    if len(f_fu)>0:
        for i in range(0,len(index),1):
            index1=index[i]
            v_final.append(v[index1])
        return v_final
    else:
        return v
    
    
def mag2flux(m):
    return 10**((m+48.6)/(-2.5))


def flux2vL(f,v,d):
    return 4*math.pi*d*d*f*v


def red2dis(z):
    return 3E8*z*3.08568E24/(6.773999999999998e4)


def dis(ra,dec,ra1,dec1):
    ra,dec,ra1,dec1=map(radians,[ra,dec,ra1,dec1])
    dra=ra-ra1
    ddec=dec-dec1
    a=sin(ddec/2)**2+cos(dec)*cos(dec1)*sin(dra/2)**2
    c=2*asin(sqrt(a))
    return c


def flux2flux(flux):
    f=10**((-73.6/2.5)+math.log10(flux))
    return f


f=open('/home/ashley/Link_sed/catalog/fig7_R06SED.dat','r')
p=f.readlines()
p0=[]
for n in range(1,len(p),1):
    p0.append(p[n])
lognu=[n.strip().split('    ')[0] for n in p0]
lognulnnu=[n.strip().split('    ')[1] for n in p0]

#delete the GEMS data===add the radio 
lamda1main=[658,900,1249,2163,3600,2.14e8]#nm    #v/ir/ir/nir/nir/'/'
fmain=[3e17/i for i in lamda1main]

"""
lamda_info get
"""
lamda_candels1=[]
lamda_ecdfs1=[]
lamda_tenis1=[]
lamda_galex1=[]
f1=open('/home/ashley/Link_sed/catalog/eff_lamda2014.2','r')
p1=f1.readlines()
for i in range(2,56,1):
    if i>=20 and i<=48:
        lamda_ecdfs1.append(p1[i])
    elif 49<=i<=50:
        lamda_tenis1.append(p1[i])
    elif 51<=i<=52:
        lamda_galex1.append(p1[i])
lamda_candels1=[]
f2=open('/home/ashley/Link_sed/catalog/eff_lamda2014.1.1','r')
p2=f2.readlines()
for i in range(0,17,1):
    lamda_candels1.append(p2[i])

lamda_candels=[n.strip().split('\t')[1] for n in lamda_candels1]
lamda_ecdfs=[n.strip().split('\t')[1] for n in lamda_ecdfs1]
lamda_tenis=[n.strip().split('\t')[1] for n in lamda_tenis1]
lamda_galex=[n.strip().split('\t')[1] for n in lamda_galex1]
f_candels=[lamda2v(n) for n in map(float,lamda_candels)]
f_ecdfs=[lamda2v(n) for n in map(float,lamda_ecdfs)]
f_tenis=[lamda2v(n) for n in map(float,lamda_tenis)]
f_galex=[lamda2v(n) for n in map(float,lamda_galex)]


#catalog read
main=fits.open('/home/ashley/Link_sed/catalog/maincat.fits')
CANDELS=fits.open('/home/ashley/Link_sed/catalog/CANDELS.GOODSS.F160W.v1.fits')
ECDFS=fits.open('/home/ashley/Link_sed/catalog/ecdfs_BVRdet_Subaru_v11.fits')
TENIS=fits.open('/home/ashley/Link_sed/catalog/TENIS.fits')
GALEX=fits.open('/home/ashley/Link_sed/catalog/CDFS_00-xd-mcat.fits')

"""
ra_dec info 
"""
ra_main=np.ndarray.tolist(main[1].data.field('RA'))
dec_main=np.ndarray.tolist(main[1].data.field('DEC'))
cpcat=np.ndarray.tolist(main[1].data.field('CP_CAT'))
ra_candels=np.ndarray.tolist(CANDELS[1].data.field('ra'))
dec_candels=np.ndarray.tolist(CANDELS[1].data.field('dec'))
ra_ecdfs=np.ndarray.tolist(ECDFS[1].data.field('RA'))
dec_ecdfs=np.ndarray.tolist(ECDFS[1].data.field('DEC'))
ra_tenis=np.ndarray.tolist(TENIS[1].data.field('ra'))
dec_tenis=np.ndarray.tolist(TENIS[1].data.field('dec'))
ra_galex=np.ndarray.tolist(GALEX[1].data.field('alpha_j2000'))
dec_galex=np.ndarray.tolist(GALEX[1].data.field('delta_j2000'))

"""
stdout
output = sys.stdout
f=open('/home/ashley/Link_sed/result/sed4.7/id','a+')
sys.stdout = f
"""

"""
color
"""
#===='tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive'
c_main=['tab:orange','tab:red','tab:pink','tab:pink','tab:brown','tab:green']#visible yellow   3.6-8.0 grey
c_candels=['tab:blue','tab:blue','tab:pink','tab:blue','tab:blue','tab:blue','tab:blue','tab:red','tab:red','tab:red','tab:pink','tab:pink','tab:pink','tab:brown','tab:brown','tab:brown','tab:brown']
c_ecdfs=['tab:purple','tab:purple','tab:orange','tab:orange','tab:orange','tab:red','tab:red','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:olive','tab:brown','tab:brown','tab:brown','tab:brown']
c_tenis=['tab:pink','tab:pink']
c_galex=['tab:purple','tab:purple']


#ECDFS APER Flux ---> total Flux
#F_TOT = F_APER * (BVR_FLUX_AUTO / F_BVR) * TOTCOR
#F_corrected=F_TOT * EXTCOR * COLCOR
TOTCOR=np.ndarray.tolist(ECDFS[1].data.field('TOTCOR'))
BVR_FLUX_AUTO=np.ndarray.tolist(ECDFS[1].data.field('BVR_FLUX_AUTO'))
F_BVR=np.ndarray.tolist(ECDFS[1].data.field('F_BVR'))

fluxcor=fits.open('/home/ashley/Link_sed/catalog/flux_correct.fits')
extcor1=np.ndarray.tolist(fluxcor[1].data.field('EXTCOR'))
colcor1=np.ndarray.tolist(fluxcor[1].data.field('COLCOR'))
EXTCOR2=[]
COLCOR2=[]
for h in range(1,33,1):
    EXTCOR2.append(extcor1[h])
    COLCOR2.append(colcor1[h])

EXTCOR=[float(n) for n in EXTCOR2]
COLCOR=[float(n) for n in COLCOR2]  
"""
judge
"""
min_dis=0.5*math.pi/(3600*180)#2.42406840554768e-06

m=1
i=39
"""
draw pic
"""
while i==39:
    
    plt.figure(m)
    min_y=[]
    """
    ra   dec need the info of catpat
    """
    if cpcat[i-1]=='CANDELS':
        ra_main_one=main[1].data[i-1][29]
        dec_main_one=main[1].data[i-1][30]
    if cpcat[i-1]=='GEMS   ':
        ra_main_one=main[1].data[i-1][26]
        dec_main_one=main[1].data[i-1][27]
    if cpcat[i-1]=='WFI    ':
        ra_main_one=main[1].data[i-1][20]
        dec_main_one=main[1].data[i-1][21]
    if cpcat[i-1]=='TENIS  ':
        ra_main_one=main[1].data[i-1][32]
        dec_main_one=main[1].data[i-1][33]
    if cpcat[i-1]=='SEDS   ':
        ra_main_one=main[1].data[i-1][35]
        dec_main_one=main[1].data[i-1][36]
    if cpcat[i-1]=='GOODS-S':
        ra_main_one=main[1].data[i-1][23]
        dec_main_one=main[1].data[i-1][24]
    if cpcat[i-1]=='VLA    ':
        ra_main_one=main[1].data[i-1][38]
        dec_main_one=main[1].data[i-1][39]
    if cpcat[i-1]=='...':
        ra_main_one=ra_main[i-1]
        dec_main_one=dec_main[i-1]

    index_main=[]
    index_candels=[]
    
    distance_main_candels=[]
    distance_main_ecdfs=[]
    distance_main_tenis=[]
    distance_main_galex=[]
    distance_main_irac=[]
    for j in range(0,len(ra_candels),1):
        ra_candels_one=ra_candels[j]
        dec_candels_one=dec_candels[j]
        distance_main_candels.append(dis(ra_main_one,dec_main_one,ra_candels_one,dec_candels_one))
    for j in range(0,len(ra_ecdfs),1):
        ra_ecdfs_one=ra_ecdfs[j]
        dec_ecdfs_one=dec_ecdfs[j]
        distance_main_ecdfs.append(dis(ra_main_one,dec_main_one,ra_ecdfs_one,dec_ecdfs_one))
    for j in range(0,len(ra_tenis),1):
        ra_tenis_one=ra_tenis[j]
        dec_tenis_one=dec_tenis[j]
        distance_main_tenis.append(dis(ra_main_one,dec_main_one,ra_tenis_one,dec_tenis_one))
    
    for j in range(0,len(ra_galex),1):
        ra_galex_one=ra_galex[j]
        dec_galex_one=dec_galex[j]
        distance_main_galex.append(dis(ra_main_one,dec_main_one,ra_galex_one,dec_galex_one))

    dis_candels=np.min(distance_main_candels)
    dis_ecdfs=np.min(distance_main_ecdfs)
    dis_tenis=np.min(distance_main_tenis)
    index_dis_tenis=distance_main_tenis.index(dis_tenis)
    dis_galex=np.min(distance_main_galex)
    
    id_mcetg=[i] # i is id

    if main[1].data[i-1][50]>0:#redshift final
        
        z_final=main[1].data[i-1][50]
        D=red2dis(z_final)
        
        mag=[]
        for j in range(22,26,3):
            mag.append(main[1].data[i-1][j])
        for j in range(31,42,3):
            mag.append(main[1].data[i-1][j])
        mag_exist=[n for n in mag if n!= -1]#choose 
        flux_main=[mag2flux(n) for n in mag_exist]
        v_main1=nummatch(mag,fmain)
        v_main=[obs2emit(n,z_final) for n in v_main1]
        c_main1=nummatch(mag,c_main)
        vLmai=[]
        for j in range(0,len(flux_main),1):
            v1mai=v_main[j]
            fluxmai=flux_main[j]
            vLmai.append(flux2vL(fluxmai,v1mai,D))
        logvmai=[math.log10(n) for n in v_main]
        logvLmai=[math.log10(n) for n in vLmai]
        if len(logvLmai)==0:
            i=i+1
            continue
        else:
            mplt.markers.MarkerStyle(fillstyle='none')
            min_y0=np.min(logvLmai)
            min_y.append(min_y0)
            star =plt.scatter(logvmai,logvLmai,marker='o',s=5,facecolors='none',edgecolor=c_main1,label='maincat')
            plt.legend(loc=1)
            

        if dis_candels <= min_dis:
            index_candels=distance_main_candels.index(dis_candels)
            id_mcetg.append(index_candels)
            flux_candels1=[]
            flux_candels0=[]
            flux_infocandels=[]
            for j in range(7,40,2):
                flux_candels1.append(CANDELS[1].data[index_candels][j])
                flux_infocandels.append(CANDELS[1].data[index_candels][j])
                flux_infocandels.append(CANDELS[1].data[index_candels][j+1])
            print('flux_candels1',flux_candels1)
            print('flux_infocandels',flux_infocandels)
            flux_candels0=[1e-29*n for n in flux_candels1 if n>0]
            v_candels=nummatch(flux_candels1,f_candels)
            c_candels1=nummatch(flux_candels1,c_candels)
            if len(flux_candels0)>0:#所有flux>0
                v_canrest=[]
                for j in range(0,len(v_candels),1):
                    v1can=obs2emit(v_candels[j],z_final)
                    v_canrest.append(v1can)#rest v>0
                vLcan=[]
                for n in range(0,len(flux_candels0),1):
                    vcan=v_canrest[n]
                    fluxcan=flux_candels0[n]
                    vLcan.append(flux2vL(fluxcan,vcan,D))
                logvcan=[math.log10(n) for n in v_canrest]
                logvLcan=[math.log10(n) for n in vLcan]
                min_y1=np.min(logvLcan)
                min_y.append(min_y1)
                mplt.markers.MarkerStyle(fillstyle='none')
                hline =plt.scatter(logvcan,logvLcan,marker='.',s=4,c=c_candels1,label='CANDELS')
                plt.legend(loc=1)
        else:
            id_mcetg.append(-1)
            min_y.append(1000)
            
        if dis_tenis <= min_dis:
            index_tenis=distance_main_tenis.index(dis_tenis)
            id_mcetg.append(index_tenis)
            flux_tenis=[]
            flux_tenis0=[]
            flux_infotenis=[]
            for j in range(3,5,1):
                flux_tenis.append(TENIS[1].data[index_tenis][j])
                flux_infotenis.append(TENIS[1].data[index_tenis][j])
            print('flux_infotenis:',flux_infotenis)
            flux_tenis0=[1e-29*n for n in flux_tenis if n>0]
            if len(flux_tenis0)>0:
                v_tenis=nummatch(flux_tenis,f_tenis)#v_tenis>0
                c_tenis1=nummatch(flux_tenis,c_tenis)
                v_tenrest=[]
                for j in range(0,len(flux_tenis0),1):
                    v1ten=obs2emit(v_tenis[j],z_final)
                    v_tenrest.append(v1ten)
                vLten=[]
                for n in range(0,len(flux_tenis0),1):
                    vten=v_tenrest[n]
                    fluxten=flux_tenis0[n]
                    vLten.append(flux2vL(fluxten,vten,D))
                logvten=[math.log10(n) for n in v_tenrest]
                logvLten=[math.log10(n) for n in vLten]
                min_y3=np.min(logvLten)
                min_y.append(min_y3)
                mplt.markers.MarkerStyle(fillstyle='none')
                plus =plt.scatter(logvten,logvLten,marker='+',s=4,color=c_tenis1,label='TENIS')
                plt.legend(loc=0)
        else:
            id_mcetg.append(-1)
            min_y.append(1000)

        if dis_ecdfs <= min_dis:
            index_ecdfs=distance_main_ecdfs.index(dis_ecdfs)
            id_mcetg.append(index_ecdfs)
            flux_ecdfs=[]
            flux_ecdfs0=[]
            totcor=TOTCOR[index_ecdfs]
            bvr_flux_auto=BVR_FLUX_AUTO[index_ecdfs]
            
            f_bvr=F_BVR[index_ecdfs]
            flux_infoecdfs=[]
            for j in range(14,28,2):
                
                flux_ecdfs.append(ECDFS[1].data[index_ecdfs][j])
                flux_infoecdfs.append(ECDFS[1].data[index_ecdfs][j])
                flux_infoecdfs.append(ECDFS[1].data[index_ecdfs][j+1])
            for j in range(34,78,2):
                
                flux_ecdfs.append(ECDFS[1].data[index_ecdfs][j])
                flux_infoecdfs.append(ECDFS[1].data[index_ecdfs][j])
                flux_infoecdfs.append(ECDFS[1].data[index_ecdfs][j+1])
            print('flux_infoecdfs',flux_infoecdfs)
            flux_ecdfs_raw=[1e-29*n for n in flux_ecdfs if n > 0]
            extcor=nummatch(flux_ecdfs,EXTCOR)
            colcor=nummatch(flux_ecdfs,COLCOR)
            v_ecdfs=nummatch(flux_ecdfs,f_ecdfs)
            c_ecdfs1=nummatch(flux_ecdfs,c_ecdfs)
            flux_ecdfs0=[]
            
            for k in range(0,len(flux_ecdfs_raw),1):
                f_tot=flux_ecdfs_raw[k]*float(bvr_flux_auto/f_bvr) * totcor
                f_corrected=f_tot * extcor[k] * colcor[k]
                flux_ecdfs0.append(f_corrected)
            '''
            flux_ecdfs0 is the flux_correct
            '''
            if len(flux_ecdfs0)>0:
                v_ecdrest=[]
                for j in range(0,len(v_ecdfs),1):
                    v1ecd=obs2emit(v_ecdfs[j],z_final)
                    v_ecdrest.append(v1ecd)
                vLecd=[]
                for n in range(0,len(flux_ecdfs0),1):
                    vecd=v_ecdrest[n]
                    fluxecd=flux_ecdfs0[n]
                    vLecd.append(flux2vL(fluxecd,vecd,D))
                logvecd=[math.log10(n) for n in v_ecdrest]
                logvLecd=[math.log10(n) for n in vLecd]
                min_y2=np.min(logvLecd)
                min_y.append(min_y2)
                mplt.markers.MarkerStyle(fillstyle='none')
                vline =plt.scatter(logvecd,logvLecd,marker='|',s=6,color=c_ecdfs1,label='ECDFS')
                plt.legend(loc=1)
        else:
            id_mcetg.append(-1)
            min_y.append(1000)            
   
        if dis_galex <= min_dis:
            index_galex=distance_main_galex.index(dis_galex)
            id_mcetg.append(index_galex)
            flux_galex=[]
            flux_galex0=[]
            flux_infogalex=[]
            for j in range(34,39,4):
                flux_galex.append(GALEX[1].data[index_galex][j])
                flux_infogalex.append(GALEX[1].data[index_galex][j])
                flux_infogalex.append(GALEX[1].data[index_galex][j+1])
            print('flux_infogalex',flux_infogalex)
            flux_galex0=[1e-29*n for n in flux_galex if n>0]
            if len(flux_galex0)>0:
                v_galex=nummatch(flux_galex,f_galex)#v_tenis>0
                c_galex1=nummatch(flux_galex,c_galex)
                v_galrest1=[]

                for j in range(0,len(flux_galex0),1):
                    v1gal=obs2emit(v_galex[j],z_final)
                    v_galrest1.append(v1gal)
                v_galrest=list(reversed(v_galrest1))
                vLgal=[]
                for n in range(0,len(flux_galex0),1):
                    vgal=v_galrest[n]
                    fluxgal=flux_galex0[n]
                    vLgal.append(flux2vL(fluxgal,vgal,D))
                logvgal=[math.log10(n) for n in v_galrest]
                logvLgal=[math.log10(n) for n in vLgal]
                min_y4=np.min(logvLgal)
                min_y.append(min_y4)
                mplt.markers.MarkerStyle(fillstyle='none')
                plus =plt.scatter(logvgal,logvLgal,marker='.',s=5,color=c_galex1,label='GALEX')
                plt.legend(loc=0)
        else:
            id_mcetg.append(-1)
            min_y.append(1000)
     
        min_yy=np.min(min_y)-1
        
        plt.plot(lognu,lognulnnu,linestyle='dashed')
        plt.xlim(13,17)
        plt.ylim(min_yy,47)
        plt.text(13.5,46,'id:'+' '+str(main[1].data[i-1][0]))
        plt.text(13.5,45.5,'z_final:'+' '+str(main[1].data[i-1][50]))

        plt.xlabel('log'+' '+'$v_{rest}$'+'(Hz)')
        plt.ylabel('log'+' '+'$v$ $L_v$'+'($erg$ $s^{-1}$)')
        picname='/home/ashley/Link_sed/result/sed4.9/'+str(main[1].data[i-1][0])+'.png'
        plt.savefig(picname,format='png',dpi=200)
        i=i+1
        m=m+1
        plt.close()
    else:
        i=i+1
main.close()
CANDELS.close()
ECDFS.close()
TENIS.close()
GALEX.close()
f.close()







