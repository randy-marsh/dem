# -*- coding: utf-8 -*-
"""
Created on Wed May 09 18:26:56 2018
Diseño de ejes
@author: Juan
"""

import numpy as np

def edsodemberg(Ta, Tm, Ma, Mm, Se, Syt, n, kfsm=2.7, kfst=2.7, kfm=2.2, kft=2.2):
    """
    Ed-Sodemberg
    :return: diameter[m]
    """
    alternate = (1.0/Se)*np.sqrt(4*np.square(kfsm*Ma) + 3*np.square(kfst*Ta))
    medium = (1.0/Syt)*np.sqrt(4*np.square(kfm*Mm) + 3*np.square(kft*Tm))
    d = np.power((16.0*n/np.pi)*(alternate + medium), 1.0/3)
    return d

def Mmax(F,L):
    """
    Bending moment suposed symetric
    """
    return F*L/4.0

def _kb(d):
    """
    Kb param for Marin ec.
    """
    return np.power(d/7.62, -0.1133)


if __name__ == '__main__':
    Ta1=178
    Ta3=25.43
    kc=0.577
    d1 = 0.051707630003
    d3 = 0.0275182151866
    Ma1 = Mmax(1766,0.11)
    Ma3 = Mmax(567.72,0.1)
    Se=0.643718*_kb(27.5420078636)*kc*0.814*264e+6
#    Se=0.643718*0.8062*0.814*264e+6
    Syt=372e+6
    n=2
    kfsm= 1 + (1/(1 + 0.08/np.sqrt(43.03)))*(2.7 - 1)
    kft= 1 + (1/(1 + 0.08/np.sqrt(43.03)))*(2.2 - 1)
    d = edsodemberg(Ta=Ta3,Tm=0.0,Ma=Ma3,Mm=0.0,Se=Se,Syt=Syt,n=n, kfsm=kfsm,
                    kft=kft)
    print d
    