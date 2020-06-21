#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 11:41:08 2020

@author: duclosivetich
"""
import matplotlib.pyplot as plt 
import numpy as np
import math
from scipy import *
from scipy . integrate import odeint 

a=2e-6
a2=8e-6
rho=1000
H=30e-6
U=260e-3
Q=260e-3/H
W=100e-6
pr0=[H/3,W/3]
pphi=np.linspace(0,10*np.pi,1000) 
rc=np.linspace(0.003,0.00125,100)
vphi=0.27                                         #Vitesse constante longitudinale
dt=0.0005
m=1050*4/3*np.pi*(a/2)**3 
m2=1050*4/3*np.pi*(a2/2)**3
rc=0.003
mu = 0.001
pr0=15e-6
pr02=15e-6


def Ul(p,pphi,rho,U,a,H):
    z2=p/H
    Cl=2.18*z2-21.23*z2*z2*z2
    return rho*U*U*a*a*a*Cl/(3*np.pi*mu*H*H)

def metpas (pr0,a):
    tabr=np.array([pr0])
    for i in range(999):
        pr0=pr0+Ul(pr0,pphi,rho,U,a,H)*dt
        tabr=np.append(tabr,pr0)
    return tabr


   
solution=odeint(Ul,pr0,pphi,args=(rho,U,a,H))
solution2=odeint(Ul,pr02,pphi,args=(rho,U,a2,H))
pos=solution[:,0]
pos2=solution2[:,0]
plt.plot(pphi,pos)
plt.plot(pphi,pos2,'r')




ym=metpas(pr0,a)
y2=metpas(pr02,a2)
x=pphi
plt.plot(x,ym,'c')
plt.plot(x,y2,'g')
plt.show()


