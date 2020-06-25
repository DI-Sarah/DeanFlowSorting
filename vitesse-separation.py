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

a=1.9e-6
a2=5.32e-6
rho=1000
U=260e-3
W=100e-6
pphi=np.linspace(0,0.5,10000) 
rc=np.linspace(0.003,0.00125,100)
vphi=0.27                                         #Vitesse constante longitudinale
dt=0.001
m=1050*4/3*np.pi*(a/2)**3 
m2=1050*4/3*np.pi*(a2/2)**3
mu = 0.001
pr02=-0.000045


def Ul(p,pphi,rho,U,a,W):
    z2=p/W
    Cl=2.18*z2-21.23*z2*z2*z2
    return rho*U*U*a*a*a*Cl/(3*np.pi*mu*W*W)

def metpas (pr0,a):
    tabr=np.array([pr0])
    for i in range(9999):
        pr0=pr0+Ul(pr0,pphi,rho,U,a,W)*dt
        tabr=np.append(tabr,pr0)
    return tabr




y2=metpas(pr02,a2)
x=pphi
plt.plot(x,y2,'g')
plt.show()


