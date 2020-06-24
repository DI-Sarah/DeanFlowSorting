#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 12:41:44 2020

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
pphi=np.linspace(0,0.05,100) 
rc=np.linspace(0.003,0.00125,100)                         #Vitesse constante longitudinale
dt=0.0000005
m=1050*4/3*np.pi*(a/2)**3 
m2=1050*4/3*np.pi*(a2/2)**3
mu = 0.001
pr0=0-0.000002
pr02=0
vr=0
vr2=0


def Fl(p,rho,U,W,r):
    z2=p/W
    Cl=2.18*z2-21.23*z2*z2*z2
    if 0-0.0000001>p and p>-4*W/10+0.0000001 or p>4*W/10+0.0000001:
        return -27*np.pi*2*r*2*2*r*2*r*r*rho*U*U*Cl/(2*H*H)
        #return -rho*Uf*Uf*Cl*2*a*2*a*2*a*2*a/H/H
    else :
        return 27*np.pi*2*2*2*2*r*r*r*r*rho*U*U*Cl/(2*H*H)
        #return rho*Uf*Uf*Cl*2*a*2*a*2*a*2*a/H/H

def Fd(p,mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    return 5.4e-4*np.pi*mu*De**(1.63)*r


def methodepas(pphi,pr0,pr02,vr,vr2,rc):
    tabr0=np.array([pr0])
    tabr02=np.array([pr02])
    tab_vr=np.array([vr])
    tab_vr2=np.array([vr2])
    for i in range(99): 
        rcb=rc[i]
        vr= (Fl(pr0,rho,Uf,H,a,W)+Fd(pr0,mu,a,rcb,rho,Uf,H))/m*dt 
        vr2= (Fl(pr02,rho,Uf,H,a2,W)+Fd(pr02,mu,a2,rcb,rho,Uf,H))/m2*dt              #Utilisation du PFD sur la particule 1
        pr0=pr0+vr*dt
        pr02=pr02+vr2*dt
        tabr0=np.append(tabr0,pr0)
        tabr02=np.append(tabr02,pr02)
        tab_vr=np.append(tab_vr,vr)
        tab_vr2=np.append(tab_vr2,vr2)
        
    return tabr0,tabr02

#axes=plt.gca()
#axes.set_ylim(-50e-6,50e-6)
'''pos,pos2=methodepas(pphi,pr0,pr02,vr,vr2,rc)
plt.plot(pphi,pos)
plt.plot(pphi,pos2,'c')'''




'''pr''(t)+Fl-Fd=0
pr'(t)=p(t)
p'(t)=Fd-Fl'''

def integrale(y,pphi,rho,U,H,a,mu,Uf,rc):
    tabdydt=np.array([])
    for i in range(100):
        rcb=rc[i]                                             #indice indiquant l'avancement de la spirale
        pr,p=y
        dydt=[p,Fl(pr0,rho,Uf,H,a,W)+Fd(pr0,mu,a,rcb,rho,Uf,H)]
        tabdydt=np.append(tabdydt,dydt)
    print(tabdydt)
    return dydt

'''axes=plt.gca()
axes.set_ylim(-10e-6,10e-6)
y0=[pr0,0.0]
sol1=odeint(integrale,y0,pphi,args=(rho,U,H,a,mu,Uf,rc))
x1=(sol1[:,0])
plt.plot(pphi,x1)
plt.show()'''



x=np.linspace(-W/2,W/2,1000)
y=Fl(x,rho,U,W,a)
z=Fd(x,mu,a,rc[0],rho,U,H)
z1=np.ones(1000)*z
plt.plot(x,z1)
plt.plot(x,y)
plt.show()()





