import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
from matplotlib.pyplot import figure, show, rc

def angles(E, v):
    ''' Conversion from E, v to theta according to the formula above '''
    k0 = E
    m = 0.51
    k = (1.-v) * E
    aux = ( (k0 + m) * k - k0 * m ) / ( k0 * k )
    return np.arccos(aux)

gamma = pp.particle.GammaDef.get()
medium = pp.medium.StandardRock(1.0)
cuts_cont = pp.EnergyCutSettings(-1, -1)

param = pp.parametrization.compton.KleinNishina(gamma, medium, cuts_cont, 1.0)

cross = pp.crosssection.ComptonIntegral(param)

fig = figure(figsize=(8, 8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],
                  projection='polar', facecolor='#d5de9c')

energies = [1e-3, 511e-3, 10][:]
#energies = [10000]
sigma_list = []
angle_list = []
v_list = []
for E in energies:
    sigma_list_tmp = []
    angle_list_tmp = []
    vlist_tmp = np.linspace(0, 1. - 1. / (1. + 2 * E / 0.51), 100000 )
    for v in vlist_tmp:
        aux = angles(E, v)
        correct = E * 0.51 / ( 0.51 - E * np.cos(aux) + E )**2 # correction factor
        if not np.isnan(aux):
            sigma_list_tmp.append(param.differential_crosssection(E, v) * correct)
            angle_list_tmp.append(angles(E, v))
    sigma_list_tmp = np.array(sigma_list_tmp)
    angle_list_tmp = np.array(angle_list_tmp)
    sigma_list_tmp = np.append(sigma_list_tmp, sigma_list_tmp)
    angle_list_tmp = np.append(angle_list_tmp, 2*np.pi - angle_list_tmp)
    sigma_list.append(sigma_list_tmp)
    angle_list.append(angle_list_tmp)
    v_list.append(vlist_tmp)


for i, (theta, r) in enumerate(zip(angle_list, sigma_list)):
    ax.plot(theta, r, label=str(energies[i]) + " MeV")

ax.legend(loc='upper left', fontsize=14)

fig.savefig("compton.pdf")