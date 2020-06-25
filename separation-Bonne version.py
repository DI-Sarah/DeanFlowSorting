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

a=1.9e-6
a2=7.32e-6
rho=1000
H=30e-6
U=260e-3
Uf=260e-3
Q=260e-3/H
W=100e-6
pphi=np.linspace(0,0.5,5000) 
rc=np.linspace(0.00305,0.0013,5000)                         #Vitesse constante longitudinale
dt=0.00001
m=1050*4/3*np.pi*(a/2)**3 
m2=1050*4/3*np.pi*(a2/2)**3
mu = 0.001
pr0b=0.00004
pr02b=0.000002
vrb=0
vr2b=0
pr=rc[0]
pr2=rc[0]

def Fl(p,rho,U,W,r):
    z2=p/W
    Cl=(np.sin(10*z2*np.pi))
    return 27*np.pi*r*r*r*r*rho*U*U*Cl/(2*W*W)
        
def Flb(p,rho,U,W,r):
    z2=p/W
    Cl=2.18*z2-21.23*z2*z2*z2
    return 27*np.pi*r*r*r*r*rho*U*U*Cl/(2*W*W)


    
def Fd(p,mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    return 5.4e-4*np.pi*mu*De**(1.63)*r


def methodepas(pphi,vrb,vr2b,pr0b,pr02b,pr,pr2):
    tabr0b=np.array([pr0b])
    tabr02b=np.array([pr02b])
    tabr=np.array([pr])
    tabr2=np.array([pr2])
    for i in range(4999): 
        rcb=rc[i]
        vrb=  (Flb(pr0b,rho,Uf,W,a)+Fd(pr0b,mu,a,rcb,rho,Uf,H))/m*dt
        vr2b= (Flb(pr02b,rho,Uf,W,a2)+Fd(pr02b,mu,a2,rcb,rho,Uf,H))/m2*dt          #Utilisation du PFD sur la particule 1
        pr0b=pr0b+vrb*dt
        pr02b=pr02b+vr2b*dt
        pr=pr0b+rcb
        pr2=pr02b+rcb
        tabr0b=np.append(tabr0b,pr0b)
        tabr02b=np.append(tabr02b,pr02b)
        tabr=np.append(tabr,pr)
        tabr2=np.append(tabr2,pr2)
        
    return tabr0b,tabr02b 


def spirale():
    tabt=np.linspace(10*np.pi,0,10001)
    tabrint=np.linspace(0.003,0.00125,10001)
    tabrext=tabrint+100e-6
    tabxI=tabrint*np.cos(tabt)
    tabyI=tabrint*np.sin(tabt)
    tabxE=tabrext*np.cos(tabt)
    tabyE=tabrext*np.sin(tabt)
    plt.plot(tabxI,tabyI,'k')
    plt.plot(tabxE,tabyE,'k')
#spirale()       
pos3,pos4=methodepas(pphi,vrb,vr2b,pr0b,pr02b,pr,pr2)
plt.plot(pphi,pos3,'r')
plt.plot(pphi,pos4,'y')
plt.show()        
       
        
        
        
        
        