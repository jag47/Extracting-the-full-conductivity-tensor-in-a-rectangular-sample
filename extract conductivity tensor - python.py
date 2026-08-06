import numpy as np
import scipy.optimize as opt
import scipy.special as sp 
from mpmath import mp

#side lengths of the sample. d1 is along the horizontal, d2 along the vertical
d1 = 2.7
d2 = 1.2

#three resistances measurments from configurations 1-3, 
corner_resistance_1 = 0.1636153846153846
corner_resistance_2 = 0.01676923076923077
corner_resistance_3 = -0.055923076923076916

#resistance measurement from configuration 5
midpoint_resistance_5 = 0.04815384615384615

#find r, hypergeometric parameter
def find_r(x):
    return np.log(1/x**2)/np.log(1-x**2)+corner_resistance_2/corner_resistance_1

sol = opt.root(find_r, 0.5)
r = sol.x[0]
print("r= "+str(r))

#calculate hall conductivity
hall_conductivity = (np.log(1-r**2)*(corner_resistance_3*np.log(1-r**2)+corner_resistance_1*np.log(1/r**2-1))/
                    ((corner_resistance_1*np.pi)**2+(corner_resistance_1*np.log(1/r**2-1)+corner_resistance_3*np.log(1-r**2))**2))
print("Hall conductivity= "+str(hall_conductivity))

#calculate geometric mean 
geomean = (-corner_resistance_1*np.log(1-r**2)*np.pi/
                    ((corner_resistance_1*np.pi)**2+(corner_resistance_1*np.log(1/r**2-1)+corner_resistance_3*np.log(1-r**2))**2))
print("Geometric mean= "+str(geomean))



#find a, hypergeometric parameter
midpoint = (1/r**2-(1-r**2)**(-midpoint_resistance_5/corner_resistance_1))/(1-(1-r**2)**(-midpoint_resistance_5/corner_resistance_1))
print("midpoint= "+str(midpoint))

def find_a(x):
    return np.pi/2*mp.hyp2f1(x, 1-x, 1, r**2)*(1-x)-mp.appellf1(1-x, 1-x, x, 2-x, midpoint, midpoint*r**2)*midpoint**(1-x)*mp.sin(np.pi*x)

a = np.float64(mp.findroot(find_a, mp.mpf(0.25)))
print("a= "+str(a))


#find alpha, the principal axes angle
def k(a,r):
    return sp.hyp2f1(a, 1-a, 1, r**2)

alpha = 1/2*np.arctan(2*d1*d2*k(a, r)*k(a, np.sqrt(1-r**2))/(d1**2*k(a, np.sqrt(1-r**2))**2-d2**2*k(a, r)**2)*np.cos(np.pi*a))
if a>1/2:
    while alpha < np.pi/2:
        alpha += np.pi/2
else:
    while alpha < 0:
        alpha += np.pi/2
print("alpha= "+str(alpha))

#get anistropy ratio
anisotropy_ratio = (-1+np.sqrt(1+np.tan(a*np.pi)**2*np.sin(2*alpha)**2))/(np.tan(a*np.pi)*np.sin(2*alpha))
print("Anisotropy ratio= "+str(anisotropy_ratio))

#get anistropic conductivities
sigma_minus = geomean*anisotropy_ratio
sigma_plus = geomean/anisotropy_ratio
print("sigma_-= "+str(sigma_minus))
print("sigma_+= "+str(sigma_plus))