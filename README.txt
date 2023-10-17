ZKDR Distance equation implemented in Python 
    - Author: Mawson William Sammons (mawson.sammons@postgrad.curtin.edu.au)
    - Please cite corresponding article : https://doi.org/10.1093/mnras/stac3013 
  
This function set provides a Python implementation of the ZKDR distance equation for both 
the model derived by Kayser 1997 in "General and Practical Method for Calculating Cosmological
distances" and Linder 1988 in "Light Propagation in Generalized Friedmann Universes". 
Functionally this code calculates the angular diameter distance at redshift z of a source 
in a universe where a fraction (eta) of the matter is distributed smoothly and the remainder is
in compact clumps far from the line of sight. For details on the motivation and application
of this distance measure see the accompanying paper (PAPER LINK). 

/////////////////
USE and Examples:
/////////////////

This function set has been designed to integrate with existing AstroPy cosmology functions and
requires the user to create their own cosmology object with their preferred parameters.

Once a cosmology object has been generated the functions contained in ZKDRDistance.py have been
imported the ZKDR angular diameter distance can be calculated by calling the relevant function
as shown below:


----Python 3.0 code------
from astropy.cosmology import Planck18 as cosmo
import ZKDRDistance as ZKDR

z1=0      #Start redshift
z2=1.0    #End redshift
eta=1.0   #fraction of matter which is smoothly distributed
method=0  #method to use (1)=Kayser,   (2)=Linder

D_eta=ZKDR.angular_diameter_distance_ZKDR(z1, z2, eta, method, cosmo)
#Answer should be 1698.8580811 Mpc
-------------------------


This gives the angular diameter distance for an observer at redshift zero in a universe where 
matter is distributed completely homogeneously. The angular diameter distance to a fictional 
observer at some non-zero redshift z1 (often required in gravitational lensing calculations)
can be calculated just as above by having a non-zero z1 value.

The answer calculated above is close (relative error ~ 1/1000) to the analytic angular diameter
distance at redshift of 1 for a homogeneous universe; 1697.81723558 Mpc (as given by 
cosmo.angular_diameter_distance(z2)). As mentioned in the accompanying paper this discrepancy 
arises due to cosmo.Om0+cosmo.Ode0 < 1 in the Planck model; the remaining density, which is in 
relativisitic matter and radiation is instead interpretted as curvature in the Kayser model, 
yielding the slightly different value. By forcing the Planck cosmology to be flat this error 
can be reduced, this is done by aggregating relativistic matter and radiation densities into 
the cold matter density parameter, as shown below:


------------------------
from astropy.cosmology import LambdaCDM
cosmo2=lambdaCDM(H0=cosmo.H0, Om0=1-cosmo.Ode0, Ode0=cosmo.Ode0)
D_etaForcedFlat=ZKDR.angular_diameter_distance_ZKDR(z1, z2, eta, method, cosmo2)
#Answer should be 1697.88049949 Mpc
----------------------- 


This answer is much closer (relative error ~1/10000).

Alternatively we can make use of Linder's generalised method which treats both relativisitic
matter and radiation densities directly.


-----------------------
method=1
D_etaLinder=ZKDR.angular_diameter_distance_ZKDR(z1, z2, eta, method, cosmo)
#Answer should be 1696.85930469 Mpc
-----------------------


However the answer has a similar relative error as the first attempt (relative error ~ 1/1000).






