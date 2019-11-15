import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt

gamma = pp.particle.GammaDef.get()
medium_list = [
    #pp.medium.Hydrogen(1.0),
    #pp.medium.Iron(1.0),
    #pp.medium.Lead(1.0),
    #pp.medium.Uranium(1.0),
    pp.medium.StandardRock(1.0),
    pp.medium.Air(1.0)
]
energy_cut_settings = pp.EnergyCutSettings(1e-10, -1)

param_list_photopair = []
param_list_compton = []
param_list_photangle = []

for medium in medium_list:
    photoangle = pp.parametrization.photopair.PhotoAngleNoDeflection(gamma, medium)
    param_list_photangle.append(photoangle)
    param_photopair = pp.parametrization.photopair.Tsai(gamma, medium, 1.0)
    param_list_photopair.append(param_photopair)
    param_compton = pp.parametrization.compton.KleinNishina(gamma, medium, energy_cut_settings, 1.0 )
    param_list_compton.append(param_compton)
    

crosssection_list_photopair = []
crosssection_list_compton = []



for param, photoangle in zip(param_list_photopair, param_list_photangle):
    cross_photopair = pp.crosssection.PhotoPairIntegral(param, photoangle)
    crosssection_list_photopair.append(cross_photopair)

for param in param_list_compton:
    cross_compton = pp.crosssection.ComptonIntegral(param)
    crosssection_list_compton.append(cross_compton)


energy_list = np.logspace(0, 5, 1000)
for cross_photopair, cross_compton, medium in zip(crosssection_list_photopair, crosssection_list_compton, medium_list):
    sigma_photopair = np.vectorize(cross_photopair.calculate_dNdx)(energy_list)
    plt.loglog(energy_list, sigma_photopair, label=r'$e$-Paarbildung')
    sigma_compton = np.vectorize(cross_compton.calculate_dNdx)(energy_list)
    plt.loglog(energy_list, sigma_compton, label="Compton")
    plt.ylabel(r'$\sigma$', fontsize=14)
    plt.xlabel(r'$E \,/\, MeV$', fontsize=14)
    plt.legend()
    plt.grid()
    plt.savefig("compton_comparison_" + str(medium.name) + ".pdf")
    plt.clf()
    
