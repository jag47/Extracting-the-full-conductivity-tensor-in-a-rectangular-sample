import numpy as np
import scipy.optimize as opt
import scipy.special as sp
import mpmath as mp

#resistances, enter in a list for configurations 1-4
# corner_resistance_1 = 0.0394754
# corner_resistance_2 = 0.0998718
# corner_resistance_3 = 0.0821355

# midpoint_resistance_1 = 0.17868

corner_resistance_1 = 0.0142/1.3
corner_resistance_2 = 0.3071/1.3
corner_resistance_3 = 0.2586/1.3

midpoint_resistance_1 = 0.438/1.3

#find r, hypergeometric parameter
def find_r(x):
    return np.log(1/x**2)/np.log(1-x**2)+corner_resistance_2/corner_resistance_1

sol = opt.root(find_r, 0.5)
r = sol.x[0]
print("r="+str(r))

#hall conductivity
hall_conductivity = (np.log(1-r**2)*(corner_resistance_3*np.log(1-r**2)+corner_resistance_1*np.log(1/r**2-1))/
                    ((corner_resistance_1*np.pi)**2+(corner_resistance_1*np.log(1/r**2-1)+corner_resistance_3*np.log(1-r**2))**2))
print("Hall conductivity="+str(hall_conductivity))

#geometric mean of anisotropic conductivity
anistropic_conductivity_geomean = (-corner_resistance_1*np.log(1-r**2)*np.pi/
                    ((corner_resistance_1*np.pi)**2+(corner_resistance_1*np.log(1/r**2-1)+corner_resistance_3*np.log(1-r**2))**2))
print("Geometric mean of anisotropic conductivity="+str(anistropic_conductivity_geomean))

#find a, hypergeometric parameter
midpoint_1 = (1/r**2-(1-r**2)**(-midpoint_resistance_1/corner_resistance_1))/(1-(1-r**2)**(-midpoint_resistance_1/corner_resistance_1))
print("midpoint_1="+str(midpoint_1))

def find_a(x):
    return np.pi/2*mp.hyp2f1(x, 1-x, 1, r**2)*(1-x)-mp.appellf1(1-x, 1-x, x, 2-x, midpoint_1, midpoint_1*r**2)*midpoint_1**(1-x)*mp.sin(np.pi*x)

a = np.float64(mp.findroot(find_a, mp.mpf(0.25)))
print("a="+str(a))


#find alpha, principal axes angle
def k(a,r):
    return sp.hyp2f1(a, 1-a, 1, r**2)

alpha = 1/2*np.arctan(2*k(a, r)*k(a, np.sqrt(1-r**2))/(k(a, np.sqrt(1-r**2))**2-k(a, r)**2)*np.cos(np.pi*a))
print("alpha="+str(alpha))

#get anistropy ratio
anisotropy_ratio = (-1+np.sqrt(1+np.tan(a*np.pi)**2*np.sin(2*alpha)**2))/(np.tan(a*np.pi)*np.sin(2*alpha))
print("Anisotropy ratio="+str(anisotropy_ratio))

#get anistropy conductivities
sigma_minus = anistropic_conductivity_geomean*anisotropy_ratio
sigma_plus = anistropic_conductivity_geomean/anisotropy_ratio
print("Anisotropic conductivity 1="+str(sigma_minus))
print("Anisotropic conductivity 2="+str(sigma_plus))