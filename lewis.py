# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 17:13:36 2018
Ecuacion de Lewis para el calculo del modulo
@author: Juan
"""
import sympy as sp
import numpy as np
import pickle
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


def selectkv(n, m, z):
    """
    dinamic factor
    
    :param n: revolutions per minute [rpm]
    :param m: module [m]
    :param z: theeth
    :return kv: kv facor
    """
    v = n*(2*np.pi/60)*(z*m)/2
    # hierro fundido, perfil moldeado
#    kv = (3.05 + v)/ 3.05
    # perfil cortado o fresado
    kv = (6.1 + v)/6.1
    # perfil generado con fresa madre o cepillado
#    kv = (3.56 + sp.sqrt(v))/3.56
    # perfil cepillado o esmerilado
#    kv = sp.sqrt((5.56 + sp.sqrt(v))/ 5.56)
    return kv


def selectF(m, n=5):
    """
    teeth width aproximation
    
    :param m: module [m]
    :param n: integer betwee 3 and 5
    :return F: teeth's width
    """
    # minimun value
    F = n*m*np.pi
    return F


def selectwt(P, n, m, z):
    """
    Gear tangential force
    
    :param P: input power [w]
    :param n: revolutions per minute [rpm]
    :param m: module [m]
    :param z: teeth  number
    """
    wt = (P/((2*np.pi*n)/60))/(z*m/2)
    return wt


def module(P, n, z, st, coef):
    """
    Calculate a module given
    
    :param P: input power [w]
    :param n: revolutions per minute [rpm]
    :param z: teeth  number
    :param st: bending strees [N/m^2]
    :param coef: safety factor
    :return sol: module [m]
    """
    Z = pickle.load(open('Z', 'rb'))
    y = pickle.load(open('y', 'rb'))
    model = make_pipeline(PolynomialFeatures(7), Ridge())
    model.fit(Z, y)
    Y = model.predict(z)
    m = sp.Symbol('m')
    kv = selectkv(n, m, z)
    F = selectF(m)
    wt = selectwt(P, n, m, z)
    sigma = kv*wt/(F*m*Y)
    print sigma
    sol = sp.solve(sigma - st/coef, m)
    return sol[0]




def selectY(z):
    """
    Lewis Shape factor
    
    :param z: teeth's number
    :returns Y: Lewis shape factor
    """
    def aproxY(z):
        """
        Polynomial aprox of Y
        
        :param z: teeth's number
        :returns Y: Lewis shape factor polynomial aproximation
        """
        model = make_pipeline(PolynomialFeatures(7), Ridge())
        model.fit(Z, y)
        return model.predict(z)

    Z = pickle.load(open('Z', 'rb'))
    y = pickle.load(open('y', 'rb'))
    if not y[np.where(Z == z)[0]].any():
        return aproxY(z)
    else:
        return y[np.where(Z == z)[0]]


def lewis(P, n, m, z, sigma, coef):
    """
    Calulate if given modules passes Lewis ec
    
    :param P: input power [w]
    :param n: revolutions per minute [rpm]
    :param m: module [m]
    :param z: teeth  number
    :param sigma: strees [N/m^2]
    :param coef: safety factor
    :return: True if passes False otherwise
    """
    Y = selectY(z)
    kv = selectkv(n, m, z)
    F = selectF(m)
    wt = selectwt(P, n, m, z)
    calsigma = kv*wt/(F*m*Y)
    if sigma/coef >= calsigma:
        return True
    else:
        return False
    
if __name__ == '__main__':
    P = 3728.5
    n = 1400
    zg = 16
    st = 274.69e+6
    coef= 4
    m = 0.004
    sol = module(P, n, zg, st, coef)
    print 'aprox module: ' + str(sol)
    if lewis(P, n, m, zg, st , coef):
        print 'module: ' + str(m) + ' mm passes'