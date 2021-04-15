#####################################################
#Author: Mawson William Sammons                     # 
#Email: mawson.sammons@postgrad.curtin.edu.au       #
#Please cite as : Sammons et al., 2021 (In Prep)    #
#####################################################



import numpy as np
import astropy.units as u
from astropy import constants as const
import scipy.integrate
import scipy.special


def DyerRoeder(z, Dset, eta, cosmo):
    """THIS FUNCTION IS A COMPONENT OF angular_diameter_distance_ZKDR FUNCTION
        
       the 2nd order differential equation to be solved for the ZKDR distance as described in 
       Kayser 1997, A General and Practical Method for Calculating Cosmological Distances
       z= cosmological redshift
       Dset = 1D array containing the current values of D and Dprime (differentiated W.R.T. 
       redshift)
       eta = the fraction of matter that is smooth
       Omega0 = the density parameter for matter at redshift z=0 (now)
       Lambda0 = the cosmological constant at z=0
       cosmo = AstroPy cosmology object containing the users preferred cosmological parameters
       FOR USE WITH SOLVE_IVP
    """
    Omega0=cosmo.Om0
    Lambda0=cosmo.Ode0
    D, Dprime = Dset
    Qval=Q(z, Omega0, Lambda0)
    Qprimeval=Qprime(z, Omega0, Lambda0)
    temp=(-1*(2*Qval/(1+z)+1/2*Qprimeval)*Dprime-3/2*eta*Omega0*(1+z)*D)/Qval
    dDsetdz=[Dprime, temp]
    return dDsetdz

def Qprime(z, Omega0, Lambda0): 
    """THIS FUNCTION IS A COMPONENT OF THE DyerRoeder FUNCTION
 
       Q parameter from Kayser 1997 differentiated w.r.t. z
    """
    return 3*Omega0*(1+z)**2-2*(Omega0+Lambda0-1)*(1+z) 

def Q(z, Omega0, Lambda0):
    """THIS FUNCTION IS A COMPONENT OF THE DyerRoeder FUNCTION
       Q parameter from Kayser 1997
    """
    return Omega0*(1+z)**3-(Omega0+Lambda0-1)*(1+z)**2+Lambda0

def GeneralisedDyerRoeder(z, Dset, eta, cosmo):
    """THIS FUNCTION IS A COMPONENT OF angular_diameter_distance_ZKDR
       E.V. Linder's generalisation of Dyer-Roeder equation from Light Propagation in 
       Generalized Friedmann Universes 1988. Generalisation comes in  the form of incorporating
       a wider range of density parameters, i.e. including relativistic matter Onu and radiation
       Ogamma
       z = cosmological redshift
       Dset = 1D array containing the current values of D and Dprime (differentiated W.R.T. 
       redshift)
       eta = the fraction of matter that is smooth
       cosmo = AstroPy cosmology object containing the users preferred cosmological parameters
       FOR USE WITH SOLVE_IVP
    """
    D, Dprime= Dset
    temp=-1*((3+decelParam(z, cosmo))/(1+z)*Dprime+3/(2*(1+z)**2)*D*(eta*cosmo.Om(z)+4/3*cosmo.Onu(z)+4/3*cosmo.Ogamma(z)))
    dDsetdz=[Dprime, temp]
    return dDsetdz

def decelParam(z, cosmo):
    """THIS FUNCTION IS A COMPONENT OF THE GeneralisedDyerRoeder FUNCTION
       Deceleration parameter (q) in E. Linder's generalised DyerRoeder equation from 
       Light Propagation in Generalized Friedmann Universes, 1988
    """
    temp=1/2*(cosmo.Om0*(1+z)+cosmo.Onu0*2*(1+z)**2+cosmo.Ogamma0*2*(1+z)**2-cosmo.Ode0*2*(1+z)**(-2))/(cosmo.Om0*(1+z)+cosmo.Onu0*(1+z)**2+cosmo.Ogamma0*(1+z)**2+cosmo.Ode0*(1+z)**(-2)-cosmo.Ok0)
    return temp

def angular_diamter_distance_ZKDR(z1, z2, eta, method, cosmo, maxStep=0.01):
    """Solves the 2nd order D.E. for the ZKDR distance from z1 to z2
       z1= start redshift
       z2= end redshift
       eta = fraction of matter which is smoothly distributed
       method = which implementation to use, 0=Kayser 1997, 1=Linder 1988 
                Kayser - Assumes a universe made of cold matter (OmegaM), dark energy (Lambda) 
                         and some curvature (1-OmegaM-Lambda)
                Linder - Assumes a generalised universe made of cold matter (OmegaM), dark 
                         energy (Lambda), relativistic matter (OmegaNu), radiation (OmegaGamma)
                         and some curvature (K)
                - NOTE: In this version these parameters will be extracted from the cosmology
                        object (cosmo.Ox)
       cosmo = AstroPy cosmology object containing the users preferred cosmological parameters
    """ 
    Dset0=[0, cosmo.H0/(cosmo.H(z1)*(1+z1))]
    if method==0: #select Kayser method 
        sol=scipy.integrate.solve_ivp(DyerRoeder, [z1,z2], Dset0, args=[eta,cosmo], max_step=maxStep)
    elif method==1: #select Linder method
        sol=scipy.integrate.solve_ivp(GeneralisedDyerRoeder, [z1,z2], Dset0, args=[eta,cosmo], max_step=maxStep)
    else:
        print('For the method parameter select either 0 : Kayser, or 1 : Linder')
    return (sol.y[0,len(sol.t)-1]*const.c/cosmo.H0).decompose()/(const.pc*1e6)*u.Mpc 
        
        
    
