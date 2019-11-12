import pyPROPOSAL as pp
import numpy as np
import matplotlib.pyplot as plt


particle_def_builder = pp.particle.ParticleDefBuilder()
particle_def_builder.SetParticleDef(pp.particle.GammaDef.get())
particle_def_builder.SetLow(0.1) # switch this line to change the e_low parameter
photon_def = particle_def_builder.build()
#photon_def = pp.particle.GammaDef.get()
#photon_def = pp.particle.EMinusDef.get()

photon_def = pp.particle.GammaDef.get()
prop = pp.Propagator(particle_def=photon_def, config_file="config_rad_length.json")
particle_to_prop = prop.particle
particle_backup = pp.particle.Particle(particle_to_prop)

statistics = 2e5
particle_backup.direction = pp.Vector3D(0, 0, -1)
particle_backup.position = pp.Vector3D(0, 0, 1000000)
particle_backup.propagated_distance = 0
particle_backup.energy = 1e5

X_0 = 36.62 # g/cm^2, PDG 2014 .
rho = 1.205e-3 # g / cm^3

L_theorie = X_0 / rho * (9./7.)

energies = np.logspace(1, 8, 12)
L_mc_list = []
for E in energies:    
    particle_backup.energy = E
    prop_length = []

    for i in range(int(statistics)):
        particle_to_prop.inject_state(particle_backup)
        list_of_secondaries = prop.propagate()
        prop_length.append(prop.particle.propagated_distance) # in cm
    L_mc_list.append(np.mean(prop_length))
    

plt.plot(energies, L_mc_list, 'x', label='PROPOSAL')
plt.axhline(y = L_theorie, label=r'$\overline{\ell}$', color='k', linestyle='dashed', linewidth=1)
plt.xscale('log')
plt.legend()
plt.grid()
plt.xlabel('E \ MeV')
plt.ylabel(r'$\overline{\ell}$ \ cm')
plt.savefig("rad_length.pdf")