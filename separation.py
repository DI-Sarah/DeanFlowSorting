import matplotlib.pyplot as plt 
import numpy as np
import math
from scipy.integrate import odeint




H=20e-6
W=100e-6                                          #Hauteur du tube
mu = 0.001 #
r=1e-6                                           #Taille particule 1
r2=2e-6                                          #Taille particule 2
a=350e-6                                 #Espace entre chaque tour de la spirale
m=1050*4/3*np.pi*(r/2)**3                         #Masse particule
rho=1000
Uf=130e-3
Umax=2*Uf
m2=1050*4/3*np.pi*(r2/2)**3

pphi=np.linspace(10*np.pi,0,100)                                   #Position initale de la particule selon phi
#rc=(a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)   
    #Rayon de courbure inital de la spirale
rc=np.linspace(0.003,0.00125,100)
vphi=0.27                                         #Vitesse constante longitudinale
dt=0.00001                                        #Pas de temps
vr=0                                              #Vitesse latérale initiale de la particule 1
vr2=0                                             #Vitesse latérale initiale de la particule 2
pr0=H/2                                         #Position initale de la particule 1 au sein du tube
pr02=H/2                                          #Position initale de la particule 2 au sein du tube
pr=0.00305                                       #Position initale de la particule 1 
pr2=rc+pr02                                       #Position initale de la particule 2







#Net Lift Force
def Fl(rho,Umax,H,r,pr):
    if pr>H/5+0.00000001:
        Cl=np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    else :
        Cl =-np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    #print(rho*Umax*Umax*Cl*r*r*r*r/H/H)
    return rho*Umax*Umax*Cl*r*r*r*r/H/H



#Dean Force
def Fd(mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    #De=0.94
    #print(5.4e-4*np.pi*mu*De**(1.63)*r)
    return 5.4e-4*np.pi*mu*De**(1.63)*r
    




#Méthode en prenant un pas de temps arbitraire et en passant de l'accélération à la vitesse en utilisant ce pas

def methodepas(pr,pphi,pr2,pr0,pr02,vr,vr2,rc):
    tabr=np.array([pr0])
    tabphi=np.array([])
    tabrc=np.array([rc])
    tabr2=np.array([pr2])
    tab_vr=np.array([0])
    tab_vr2=np.array([0])
    tab_Fl=np.array([0])
    tab_Fd=np.array([0])
    tabr0=np.array([pr0])
    for i in range(99): 
        rcb=rc[i]
        pphib=pphi[i]
        vr= vr-Fl(rho,Umax,H,r,pr0)/m*dt+Fd(mu,r,rcb,rho,Uf,H)/m*dt              #Utilisation du PFD sur la particule 1
        vr2=vr2-Fl(rho,Umax,H,r2,pr02)/m2*dt+Fd(mu,r2,rcb,rho,Uf,H)/m2*dt        #Utilisation du PFD sur la particule 2
        pr0=pr0+vr*dt
        pr02=pr02+vr2*dt
        #rc= (a/(2*np.pi))*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)
        pr=pr0+rcb
        pr2=pr02+rcb
        #Mise à jour des tableaux à chaque passage
        tabrc=np.append(tabrc,rcb)
        tabr=np.append(tabr,pr0)
        tabphi=np.append(tabphi,pphib)
        #tabr2=np.append(tabr2,pr2)
        tab_vr=np.append(tab_vr,vr)
        tab_vr2=np.append(tab_vr2,vr2)
        tab_Fl=np.append(tab_Fl,Fl(rho,Umax,H,r,pr0))
        tab_Fd=np.append(tab_Fd,Fd(mu,r,rc,rho,Uf,H))

    return(tabr,tabr2,tabphi)



#Tracé la spirale d'archimède de 5 tours à bonne échelle

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
   
    #plt.plot(tabxM,tabyM,'k')



#Méthode en utilisant odeint qui permet d'intégrer notre équation différentielle : pr''(t)+Fl-Fd=0 où pr est la position de la particule
'''pr''(t)+Fl-Fd=0
pr'(t)=p(t)
p'(t)=Fd-Fl'''
def integrale(y,t,rho,Umax,H,r,mu,Uf,rc):
    for i in range(100):
        rcb=rc[i]                                             #indice indiquant l'avancement de la spirale
        pr,p=y
        dydt=[p,Fd(mu,r,rc[i],rho,Uf,H)-Fl(rho,Umax,H,r,pr0)]
    return dydt
   




#Fonction principale permettant d'appeler les fonctions annexes - execution du programme
def main(n):
    #spirale()
    if n==1 :                                                       #Utilisation de la méthode pas
        tab1,tab2,tab3=methodepas(pr,pphi,pr2,pr0,pr02,vr,vr2,rc)
        #axes=plt.gca()
        #axes.set_xlim(-0.004,0.004)
        #axes.set_ylim(-0.004,0.004)
        x=tab1
        y=pphi
        #x2=tab2*np.cos(tab3)
        plt.plot(y,x,"g")
    if n==2 :                                                       #Utilisation de la méthode intégrale
        y0=[pr0,0.0]
        sol1=odeint(integrale,y0,pphi,args=(rho,Umax,H,r,mu,Uf,rc))
        sol2=odeint(integrale,y0,pphi,args=(rho,Umax,H,r2,mu,Uf,rc))
        x1=(sol1[:,0]+rc)*np.cos(pphi)
        y1=(sol1[:,0]+rc)*np.sin(pphi)
        #x=phi
        #y=sol1[:,0]
        #y2=sol2[:,0]
        #x2=(sol2[:,0]+rc2)*np.cos(t)
        #y2=(sol2[:,0]+rc2)*np.sin(t)
        plt.plot(x1,y1,"r")
        #plt.plot(x,y2,"b")


main(1)
plt.show()
