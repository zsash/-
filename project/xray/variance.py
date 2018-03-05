#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:37:01 2017

@author: ashley
"""

#==================this code is used to calculate the variance of one galaxy and then draw the diatance_variance pic.
#==================the asic_num>=2 galaxy can calculate the variance
#======nan
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import pandas as pa
from matplotlib.patches import Rectangle
import matplotlib.patches as pt
import math
"""
file read
"""
#file
file=fits.open('/mnt/wwn-0x5001b444a4ec7c48/xray/data/resultge2.fits')#match all the data from SDSS and CXO
name=np.ndarray.tolist(file[1].data.field('name_1'))
redshift=np.ndarray.tolist(file[1].data.field('z'))
acis_num=np.ndarray.tolist(file[1].data.field('acis_num_1'))
#file1
file1=fits.open('/mnt/wwn-0x5001b444a4ec7c48/xray/data/dataforyang2.fits')#3770=acis_num>=1   830 acis_num>=2
flux=np.ndarray.tolist(file1[1].data.field('flux_aper_b'))
flux_hilim=np.ndarray.tolist(file1[1].data.field('flux_aper_hilim_b'))
flux_lolim=np.ndarray.tolist(file1[1].data.field('flux_aper_lolim_b'))
name1=np.ndarray.tolist(file1[1].data.field('name'))
acis_num1=np.ndarray.tolist(file1[1].data.field('acis_num'))


#=============pick up the name of the acis_num >=2
#=============svar >3 xianzhu guangbian 
name_2=[]
redshift_2=[]
acis_2=[]
for i in range(0,len(name),1):
    if acis_num[i] >= 2:#830
        #print(i)
        name_2.append(name[i])
        redshift_2.append(redshift[i])
        acis_2.append(acis_num[i])
    else:
        i+=1

D=[]
fluxratio=[]  
variance=[]
mean_L=[]
max_L=[]
min_L=[]
excess_var=[]
id_fluxeq0=[]
Sd2=[]
svar=[]
#print('len(name_3)',len(name_3))

for i in range(0,len(name_2),1):
    indexname=name1.index(name_2[i])
    acisnum=acis_2[i]
    redshift_tmp=redshift[i]
    
    dist=redshift_tmp*3E8*308568E16/70#cm
    dist_mpc=dist*3.24077929e-25
    D.append(dist_mpc)
    tmp_flux=[]
    #mp_flux_h=[]
    uncertainty_flux=[]
    tmpstd_flux=[]
    std=[]
    if math.isnan(flux[indexname])==1 and math.isnan(flux_lolim[indexname])==1 and math.isnan(flux_hilim[indexname])==1:
        while math.isnan(flux[indexname])==1:
            
            indexname+=1
    
    for j in range(indexname,indexname + acisnum,1):
        if flux[j]==0.0:
            id_fluxeq0.append(j)
            tmp_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(0.0)
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
            uncertainty=flux_hilim[j]/2
            uncertainty_flux.append(uncertainty)
         #   j=j+1
        elif flux[j]!=0.0 and math.isnan(flux_lolim[j])==1:
            tmp_flux.append(flux[j])
            uncertainty=flux_hilim[j]-flux[j]
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux[j])
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
            uncertainty_flux.append(uncertainty)
        else:
            tmp_flux.append(flux[j])    
            uncertainty=(flux_hilim[j]-flux_lolim[j])/2
            uncertainty_flux.append(uncertainty)
            tmpstd_flux.append(flux_hilim[j])
            tmpstd_flux.append(flux[j])
            tmpstd_flux.append(flux_lolim[j])
            std_tmp=np.std(tmpstd_flux)
            std.append(std_tmp)
    #print('i',i)
    #print('j',j)
    max_flux=np.max(tmp_flux)
    index_max=tmp_flux.index(max_flux)
    min_flux=np.min(tmp_flux)
    index_min=tmp_flux.index(min_flux)
    err_min1=uncertainty_flux[index_min]
    err_max1=uncertainty_flux[index_max]
    err_min=std[index_min]
    err_max=std[index_max] #Fmax and Fmin are the maximum and minimum EPIC fluxes of a unique source, with the corresponding (statistica
    svar_tmp=(max_flux-min_flux)/np.sqrt(err_min1**2+err_max1**2)
    svar.append(svar_tmp)
    
    ratio=max_flux/min_flux
    fluxratio.append(ratio)
    mean_flux=np.mean(tmp_flux)

    
    max_l=max_flux*dist**2*4*np.pi
    mean_l=mean_flux*dist**2*4*np.pi
    min_l=min_flux*dist**2*4*np.pi
    
    mean_L.append(mean_l)
    log_meanl=[np.log10(n) for n in mean_L]
    max_L.append(max_l)
    min_L.append(min_l)
    tmp_excessvar=0
    for n in range(0,len(tmp_flux),1):
        tmp_excessvar+=((tmp_flux[n]-mean_flux)**2-(uncertainty_flux[n])**2)/((len(tmp_flux))*(mean_flux)**2)
    excess_var.append(tmp_excessvar)# 
    if tmp_excessvar>=1:
        print('indexname',indexname)
        print(tmp_excessvar)
        print(tmp_flux)
        print(uncertainty_flux)
        print(mean_flux)
        print(i)
    tmp_sd1=0
    for n in range(0,len(tmp_flux),1):
        tmp_sd1+=((tmp_flux[n]-mean_flux)**2-(std[n])**2-tmp_excessvar)**2*mean_flux**2/(len(tmp_flux)-1)
    tmp_sd2=tmp_sd1/(mean_flux**2*np.sqrt(len(tmp_flux)))
    Sd2.append(tmp_sd2)
svargt3=[n for n in svar if n >=3]
print('num',len(svargt3))
index_info=[]
for i in range(0,len(svargt3),1):
    index_info.append(svar.index(svargt3[i]))
max_max=np.max(max_L)
index=max_L.index(max_max)
#print('index_maxl',index)

#print('svar',svargt3)   
Sd=[np.sqrt(n) for n in Sd2]    
meanl_div1040=[n/(10**40) for n in mean_L]
#sqrt_rms=[np.sqrt(n) for n in excess_var]
    #print(excess_var)
logmeanl=[np.log10(n) for n in mean_L]
logmeanlgt3=[]
excessgt3=[]
for n in range(0,len(index_info),1):
    logmeanlgt3.append(logmeanl[index_info[n]])
    excessgt3.append(excess_var[index_info[n]])

#=====change y to square
sigma1=[];sigma2=[];sigma3=[];sigma4=[];sigma5=[];sigma6=[];sigma7=[];sigma8=[];sigma9=[]
sigma10=[];sigma11=[]
meanl1=[];meanl2=[];meanl3=[];meanl4=[];meanl5=[];meanl6=[];meanl7=[];meanl8=[];meanl9=[];
meanl10=[];meanl11=[]
evar1=[];evar2=[];evar3=[];evar4=[];evar5=[];evar6=[]
for n in range(0,len(logmeanl),1):
    if 40>=logmeanl[n]>=39 and svar[n]>=3:
        sigma1.append(Sd2[n])
        meanl1.append(logmeanl[n])
        evar1.append(excess_var[n])
 
    if 41>=logmeanl[n]>=40 and svar[n]>=3:
        sigma2.append(Sd2[n])
        meanl2.append(logmeanl[n])
        evar2.append(excess_var[n])
 
    if 42>=logmeanl[n]>=41 and svar[n]>=3:
        sigma3.append(Sd2[n])
        meanl3.append(logmeanl[n])
        evar3.append(excess_var[n])
  
    if 43>=logmeanl[n]>=42 and svar[n]>=3:
        sigma4.append(Sd2[n])
        meanl4.append(logmeanl[n])
        evar4.append(excess_var[n])
   
    if 44>=logmeanl[n]>=43 and svar[n]>=3:
        sigma5.append(Sd2[n])
        meanl5.append(logmeanl[n])
        evar5.append(excess_var[n])
   
    if 45>=logmeanl[n]>=44 and svar[n]>=3:
        sigma6.append(Sd2[n])
        meanl6.append(logmeanl[n])
        evar6.append(excess_var[n])
mean_sigma1=np.sum(sigma1)/len(sigma1)
mean_evar1=np.sum(evar1)/len(evar1)
tmp_err1=0;tmp_err2=0;tmp_err3=0;tmp_err4=0;tmp_err5=0;tmp_err6=0;
print(len(evar1))
for n in range(0,len(evar1),1):
    tmp_err1+=(evar1[n]-mean_evar1)**2/(len(evar1)*(len(evar1)-1))
#print('tmp_err1',tmp_err1)
sqr_err1=np.sqrt(tmp_err1)
#print('sqr_err1',sqr_err1)
mean_ml1=np.mean(meanl1)

mean_sigma2=np.sum(sigma2)/len(sigma2)
mean_evar2=np.sum(evar2)/len(evar2)
print(len(evar2))
for n in range(0,len(evar2),1):
    tmp_err2+=(evar2[n]-mean_evar2)**2/(len(evar2)*(len(evar2)-1))
sqr_err2=np.sqrt(tmp_err2)
mean_ml2=np.mean(meanl2)

mean_sigma3=np.sum(sigma3)/len(sigma3)
mean_evar3=np.sum(evar3)/len(evar3)
for n in range(0,len(evar3),1):
    tmp_err3+=(evar3[n]-mean_evar3)**2/(len(evar3)*(len(evar3)-1))
sqr_err3=np.sqrt( tmp_err3)
mean_ml3=np.mean(meanl3)

mean_sigma4=np.sum(sigma4)/len(sigma4)
mean_evar4=np.sum(evar4)/len(evar4)
for n in range(0,len(evar4),1):
    tmp_err4+=(evar4[n]-mean_evar4)**2/(len(evar4)*(len(evar4)-1))
sqr_err4=np.sqrt(tmp_err4) 
mean_ml4=np.mean(meanl4)

mean_sigma5=np.sum(sigma5)/len(sigma5)
mean_evar5=np.sum(evar5)/len(evar5)
for n in range(0,len(evar5),1):
    tmp_err5+=(evar5[n]-mean_evar5)**2/(len(evar5)*(len(evar5)-1))
sqr_err5=np.sqrt(tmp_err5) 
mean_ml5=np.mean(meanl5)

mean_sigma6=np.sum(sigma6)/len(sigma6)
mean_evar6=np.sum(evar6)/len(evar6)
for n in range(0,len(evar6),1):
    tmp_err6+=(evar6[n]-mean_evar6)**2/(len(evar6)*(len(evar6)-1))
sqr_err6=np.sqrt(tmp_err6)
mean_ml6=np.mean(meanl6)

xlolim=[40-mean_ml1,41-mean_ml2,42-mean_ml3,43-mean_ml4,44-mean_ml5,45-mean_ml6]#need to change 
#print(len(xuplim))
xuplim=[mean_ml1-39,mean_ml2-40,mean_ml3-41,mean_ml4-42,mean_ml5-43,mean_ml6-44]
yuplim=[sqr_err1,sqr_err2,sqr_err3,sqr_err4,sqr_err5,sqr_err6]

ylolim=[sqr_err1,sqr_err2,sqr_err3,sqr_err4,sqr_err5,sqr_err6]

x=[mean_ml1,mean_ml2,mean_ml3,mean_ml4,mean_ml5,mean_ml6]
y=[mean_evar1,mean_evar2,mean_evar3,mean_evar4,mean_evar5,mean_evar6]
  

B=pt.ArrowStyle("Fancy", head_length=.4, head_width=.4, tail_width=.4)
a=plt.scatter(logmeanlgt3,excessgt3,s=4,facecolors='none',edgecolor='grey',label='acis_num>=2')

plt.legend(loc=0)
#
plt.scatter(x,y,s=4,facecolor='none',edgecolor='r')
plt.errorbar(x,y,ms=4,yerr=ylolim,lolims=True,xlolims=True,xerr=xlolim,markerfacecolor='none',markeredgecolor='r',fmt='o',c='r')
plt.errorbar(x,y,ms=4,uplims=True,yerr=yuplim,xuplims=True,xerr=xuplim,markerfacecolor='none',markeredgecolor='r',fmt='o',c='r')

plt.xlabel('Mean '+'$logL_{x}$'+' '+'  '+'(0.5-7keV)'+' '+'(erg $s^{-1}$)')
plt.ylabel('$\sigma_{rms}^2$')

#plt.ylim(-0.05,0.3)
plt.savefig('/mnt/wwn-0x5001b444a4ec7c48/xray/result/excess_var&svargt31.png',dpi=400,format='png')  

ratiogt10=[]
ratiogt5=[]
ratiolt5=[]
minlgt10=[]
minlgt5=[]
minllt5=[]
for i in range(0,len(fluxratio),1):
    if fluxratio[i]>=10:
        ratiogt10.append(fluxratio[i])
        minlgt10.append(min_L[i])
    if fluxratio[i]>=5:
        ratiogt5.append(fluxratio[i])
        minlgt5.append(min_L[i])
    if fluxratio[i]<5:
        ratiolt5.append(fluxratio[i])
        minllt5.append(min_L[i])
max_fluxratio=np.max(fluxratio)
index1=fluxratio.index(max_fluxratio)
print(index1)
logminl=[np.log10(n) for n in min_L]
logminlgt5=[np.log10(n) for n in minlgt5]
logminllt5=[np.log10(n) for n in minllt5]
sminl=pa.Series(logminl)
sminlgt5=pa.Series(logminlgt5)
sminllt5=pa.Series(logminllt5)

log_minlgt10=[np.log10(n) for n in minlgt10]
plt.figure()
Sminl=sminl.dropna()

Sminlgt5=sminlgt5.dropna()
Sminllt5=sminllt5.dropna()
bin1=[39,39.2,39.4,39.6,39.8,40,40.2,40.4,40.6,40.8,41, 41.2, 41.4, 41.6, 41.8,42, 42.2, 42.4, 42.6, 42.8,43, 43.2, 43.4, 43.6, 43.8,44, 44.2, 44.4, 44.6, 44.8]

plt.hist(Sminl,bins=bin1,facecolor='none',edgecolor='g')
plt.hist(log_minlgt10,bins=bin1,facecolor='none',edgecolor='b') 
handels=[Rectangle((0,0),1,1,color=c) for c in ['g','b']]
labels=["all","ratio>10"]
plt.legend(handels,labels)   
plt.xlabel('All log_minL')
plt.title('fluxratio')
plt.yscale('log')
plt.savefig('/mnt/wwn-0x5001b444a4ec7c48/xray/result/histall.png',dpi=200,format='png')

plt.figure()
plt.hist(Sminlgt5,bins=bin1,facecolor='none',edgecolor='r')
plt.hist(Sminllt5,bins=bin1,facecolor='none',edgecolor='y')
handels1=[Rectangle((0,0),1,1,color=c) for c in ['r','y']]
labels1=["ratio>=5","ratio<5"]
plt.legend(handels1,labels1)
plt.xlabel('log_minL')
plt.title('flux_ratio')
plt.yscale('log')
plt.savefig('/mnt/wwn-0x5001b444a4ec7c48/xray/result/hist_div5.png',dpi=200,format='png')

file.close()
file1.close()
