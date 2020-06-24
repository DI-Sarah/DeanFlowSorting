#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:51:33 2020

@author: duclosivetich
"""


from scipy import *
from scipy . integrate import odeint 
import matplotlib.pyplot as plt 
import numpy as np
import math

a=2.3e-6
a2=7.32e-6
rho=1000
H=30e-6
U=260e-3
Uf=260e-3
Q=260e-3/H
W=100e-6
pphi=np.linspace(0,0.05,5000) 
rc=np.linspace(0.003,0.00125,5000)                         #Vitesse constante longitudinale
dt=0.00001
m=1050*4/3*np.pi*(a/2)**3 
m2=1050*4/3*np.pi*(a2/2)**3
mu = 0.001
pr0=-5e-6
pr02=-5e-6
vr=0
vr2=0

'''def Fl(p,rho,U,r,W):
    z2=p/W
    Cl=np.abs(np.sin(10*np.pi/3/W*p-5*np.pi/3))
    if (z2<1/2 and z2>1/5) or z2>4/5 :
        
        return -rho*U*U*Cl*16*r*r*r*r/W/W
        #return -rho*Uf*Uf*Cl*2*a*2*a*2*a*2*a/H/H
    else :
        return rho*U*U*Cl*16*r*r*r*r/W/W'''
        

def Fl(p,rho,U,W,r):
    z2=p/W
    Cl=2.18*z2-21.23*z2*z2*z2
    return 27*np.pi*r*r*r*r*rho*U*U*Cl/(2*W*W)
    
def Fd(p,mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    return 5.4e-4*np.pi*mu*De**(1.63)*r


def methodepas(pphi,pr0,pr02,vr,vr2):
    tabr0=np.array([pr0])
    tabr02=np.array([pr02])
    tab_vr=np.array([vr])
    tab_vr2=np.array([vr2])
    for i in range(4999): 
        rcb=rc[i]
        vr= ( Fl(pr0,rho,Uf,W,a)-Fd(pr0,mu,a,rcb,rho,Uf,H))/m*dt
        vr2= (Fl(pr02,rho,Uf,W,a2)-Fd(pr02,mu,a2,rcb,rho,Uf,H))/m2*dt              #Utilisation du PFD sur la particule 1
        pr0=pr0+vr*dt
        pr02=pr02+vr2*dt
        tabr0=np.append(tabr0,pr0)
        tabr02=np.append(tabr02,pr02)
        tab_vr=np.append(tab_vr,vr)
        tab_vr2=np.append(tab_vr2,vr2)
        
    return tabr0,tabr02  
        
        
pos,pos2=methodepas(pphi,pr0,pr02,vr,vr2)
plt.plot(pphi,pos)
plt.plot(pphi,pos2,'c')        
        
        
        
        
        
        
        