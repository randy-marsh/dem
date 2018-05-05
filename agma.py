# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 17:45:45 2018
Procedimientos AGMA
@author: Juan
"""
import numpy as np
import sympy as sp
import warnings

from lewis import selectF, selectY, selectwt, selectkv

def _ks(F, Y, P):
    """
    size factor (ks)
    
    :param F: teeth width [m]
    :param Y: Lewis shape factor
    :param P: paso diamteral 1/m [m]^-1
    :return ks: ks factor
    """
    ks = 1.192*np.power(F*np.sqrt(Y)/P, 0.0535)
    if ks < 1.0:
        return 1.0
    else:
        return ks
    

def _kv(n, m, z, Qv):
    """
    dinamic factor
    
    :param n: revolutions per minute [rpm]
    :param m: module [m]
    :param z: theeth
    :param Qv: Qv number
    :return kv: kv factor
    """
    v = n*(2*np.pi/60)*(z*m)/2
    B = 0.25*np.power(12 - Qv, 2.0/3.0)
    A = 50 + 56*(1-B)
    kv = np.power(((A + np.sqrt(200*v))/A), B)
    return kv

### Kh
def selectcpf(F,d):
    """
    
    :param F: teeth width
    :param d: primitive diameter m*z
    """
    # conversion to inches
    F = F*(5.0/127)
    d = d*(5.0/127)
    # from page 17
    if F <= 1.0:
        return F/(10*d) - 0.025
    elif 1.0 < F <= 17:
        return F/(10*d) - 0.0375 + 0.0125*F
    else:
        return F/(10*d) - 0.0207*F + 0.000228*np.square(F)


def selectcma(F):
    # unidades comerciales cerradas
    A = 0.127
    B = 0.0158
    C = -0.930e-4
    return A +B*F + C*np.square(F)


def _kh( cpf, cma, cmc=0.8, cpm=1.0,  ce=1):
    """
    kh factor de distribucion de carga
    
    """
    return 1 + cmc*(cpf*cpm + cma*ce)

def sigmaagma(wt,ko,kv,F,m,kh,ks,kb,yj):
    """
    AGMA bending stress
    
    :param wt: tangential load [N]
    :param ko: overload factor
    :param kv: dinamic factor
    :param ks: size factor
    :param F: teeth width [m]
    :param m: gear module [m]
    :param kh: factor de distribucion de carga
    :param kb: factor de espesor del aro
    :param yj: yjgeometric factor
    :return sigma: AGMA bending stress [N/m^2]
    """
    return (wt*ko*kv*ks*kh*kb)/(F*m*yj)

def _zn(N):
    """
    life factor suposed worst case and N > 10e+7
    
    :param N: load circles
    :return zn: life factor
    """
    if N < 10e+7:
        warnings.warn("life circles below 10e+7 calculation may be accurate!!")
    return 1.4488*np.power(N, -0.023)


def _yn(N):
    """
    life factor suposed worst case and N > 10e+7
     :param N: load circles
    :return y: life factor
    """
    if N < 10e+6:
        warnings.warn("life circles below 10e+6 calculation may be accurate!!")
    return 1.6831*np.power(N, -0.0323)


def _yz(R):
    """
    """
    if 0.5 < R < 0.99:
        return 0.658 - 0.0759*np.log(1- R)
    elif 0.99 <= R <= 0.9999:
        return 0.5 - 0.109*np.log(1 - R)


def _sf(st, yn, sigma, yo=1.0, yz=1.0):
    """
    AGMA stress safety factor
    
    :param st: allowable stress [N/m^2]
    :param yn:
    :param yo: temperature factor
    :param yz: reliabitily factor
    :parma sigma: AGMA stress [N/m^2]
    returns safety factor
    """
    return ((st*yn)/(yo*yz))/(sigma)


def _sfc(sc, yo, yz, zn, sigmac, zw=1.0):
    """
    AGMA contact safety factor
    
    :param sc: allowable stress [N/m^2]
    :param yo: temperature factor
    :param yz: reliabitily factor
    :param zw: hardness factor
    :param sigmac: AGMA contact shear stress [N/m^2]
    :returns sfc: AGMA contact safety factor
    """
    return (sc*zn*zw)/(yo*yz*sigmac)
    

def _ze(E=207e+9, poisson=0.292):
    """
    suposed same material
    [E/(2π*(1-v^2))]^(1/2)
    
    :param E: Young module Pa
    :param poisson: poisson
    """
    return np.sqrt(E/(2*np.pi*(1-np.square(poisson))))
    
def _zi(phi_t,mg,mn=1.0):
    """
    geometric factor
    Z_i=(cos Φ_t senΦ_t)/(2m_N )  m_G/(m_G+1)
    
    :param phi_t: pressure angle [deg]
    :param mg: gear ratio Ng/Nd
    :param mn: load balance factor
    """
    return (np.cos(np.radians(phi_t))*np.sin(np.radians(phi_t))*mg)/(2*mn*(mg-1))

def _sigmac(ze,wt,ko,kv,ks,kh,dp,F,zi,zr=1.0):
    """
    AGMA contact stress
    
    :param ze:
    :param wt: tangential load [N]
    :param ko: overload factor
    :param kv: dinamic factor
    :param ks: size factor
    :param kh: factor de distribucion de carga
    :param dp: primitive diameter [m]
    :param F: teeth width [m]
    :param zi: geometric factor
    :param zr: superfitial factor
    :return sigma: AGMA contact stress [N/m^2]
    """
#    return ze*np.sqrt((wt*ko*kv*ks*kh*zr)/(dp*F*m*zi))
    return ze*np.sqrt((wt*ko*kv*ks*kh*zr)/(dp*F*zi))

class agma:
    """
    AGMA
    """
    
    def __init__(self, power, rpm, z, m, mg, yj, phi=20, ko=1.0, kb=1.0, Qv=6,
                 el=50e+3, yo=1.0, yz=1.0, zr=1.0, zw=1.0, Fwidth=5):
        self.power = power
        self.rpm = n
        self.z = z
        self.m = m
        self.dp = self.m*self.z
        self.F = selectF(self.m, n=Fwidth)
        self.Y = selectY(self.z)
        self.P = 1./m
        self.ks = _ks(self.F, z, self.P)
        self.wt = selectwt(power, rpm, m, z)
        self.kv = _kv(rpm, m, z, Qv)
        self.ko = ko
        self.kb = kb
        self.yj = yj
        self.zi = _zi(phi_t=phi, mg=mg)
        self.ze = _ze()
        self.kh = _kh(cpf=selectcpf(self.F, self.dp), cma=selectcma(self.F))
        self.yo = yo
        self.yz = yz
        self.zr = zr
        self.zw = 1.0
        self.zn = _zn(el*60*rpm)
        self.yn = _yn(el*60*rpm)
        pass


    def sigmab(self):
        """
        """

        sigmabs = sigmaagma(wt=self.wt, ko=self.ko, kv=self.kv, F=self.F,
                           m=self.m, kh=self.kh, ks=self.ks, kb=self.kb, yj=self.yj)
        print 'AGMA bending stress: ' + str(sigmabs) + ' Pa'
        return sigmabs


    def sigmabsf(self, st, sigmabs):
        """
        """
#        (st, yn, sigma, yo=1.0, yz=1.0
        sigmaabssf = _sf(st=st, sigma=sigmabs, yn=self.yn, yo=self.yo,
                      yz=self.yz) 
        print 'AGMA bending stress safety factor: ' + str(sigmaabssf)
        return sigmaabssf

    def sigmac(self):
        """
        ze,wt,ko,kv,ks,kh,dp,F,zi,zr=1.0
        """
        acs = _sigmac(self.ze, self.wt, self.ko, self.kv, self.ks, self.kh,
                      self.dp, self.F, self.zi, self.zr)
        self.acs = acs
        print 'AGMA constact stress: ' + str(acs) + ' Pa'
        return acs

    def sigmacsf(self, sc, acs):
        """
        sc, yo, yz, zn, sigmac, zw=1.0
        """
        acssf = _sfc(sc=sc, sigmac=acs, yo=self.yo, yz=self.yz, zn=self.zn,
                     zw=self.zw)
        self.acssf = acssf
        print 'AGMA contact stress safety factor: ' + str(acssf)
        return acssf
    
    def info(self):
        print 'z: ' + str(self.z) +' teeth, module: ' + str(self.m) + ' m F: ' + str(self.F) + ' m' 
    
        
if __name__ == '__main__':
    # ZG1
    zg1 =agma(power=3728.5, rpm=200, z=56, m=0.004, mg=3.5, yj=0.4)
    zg1.info()
    sigmabs = zg1.sigmab()
    zg1.sigmabsf(st=274.69e+6, sigmabs=sigmabs)
    acs = zg1.sigmac()
    zg1.sigmacsf(sc=742.3e+6, acs=acs)
    ###########################################################################
    #ZP1
    zp1 =agma(power=3728.5, rpm=700, z=16.0, m=0.004, mg=3.5, yj=0.27)
    zp1.info()
    sigmabs = zp1.sigmab()
    zp1.sigmabsf(st=274.69e+6, sigmabs=sigmabs)
    acs = zp1.sigmac()
    zp1.sigmacsf(sc=742.3e+6, acs=acs)
    ###########################################################################
#    ZG2
    zg2 =agma(power=3728.5, rpm=1400, z=32.0, m=0.006, mg=2.0, yj=0.36, Fwidth=3)
    zg2.info()
    sigmabs= zg2.sigmab()
    zg2.sigmabsf(st=274.69e+6, sigmabs=sigmabs)
    acs = zg2.sigmac()
    zg2.sigmacsf(sc=742.3e+6, acs=acs)
    ###########################################################################
    #ZP2
    zp2 =agma(power=3728.5, rpm=700, z=16.0, m=0.006, mg=2.0, yj=0.27, Fwidth=3)
    zp2.info()
    sigmabs= zp2.sigmab()
    zp2.sigmabsf(st=274.69e+6, sigmabs=sigmabs)
    acs = zp2.sigmac()
    zp2.sigmacsf(sc=742.3e+6, acs=acs)