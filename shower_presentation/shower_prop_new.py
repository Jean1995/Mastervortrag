from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pyPROPOSAL as pp
from enum import Enum
import sys

photon_def = pp.particle.GammaDef.get()
prop_photon = pp.Propagator(particle_def=photon_def, config_file="shower/config_photon.json")

electron_def = pp.particle.EMinusDef.get()
prop_electron = pp.Propagator(particle_def=electron_def, config_file="shower/config_electron.json")

positron_def = pp.particle.EPlusDef.get()
prop_positron = pp.Propagator(particle_def=positron_def, config_file="shower/config_positron.json")

def propagate(propagator, position, direction, energy, time, propagated_distance):
    particle_backup = pp.particle.Particle(propagator.particle)
    
    ## assign values
    particle_backup.position = position
    particle_backup.direction = direction
    particle_backup.energy = energy
    particle_backup.time = time
    particle_backup.propagated_distance = propagated_distance
    
    propagator.particle.inject_state(particle_backup)
    list_of_secondaries = propagator.propagate()
    
    return (list_of_secondaries, propagator.particle.position)


def propagate_photon(position, direction, energy, time, propagated_distance):    
    return propagate(prop_photon, position, direction, energy, time, propagated_distance)

def propagate_electron(position, direction, energy, time, propagated_distance):    
    return propagate(prop_electron, position, direction, energy, time, propagated_distance)

def propagate_positron(position, direction, energy, time, propagated_distance):    
    return propagate(prop_positron, position, direction, energy, time, propagated_distance)

class Particle(Enum):
    PHOTON = 0
    ELECTRON = 1
    POSITRON = 2
    BREMS = 3
    
class ParticleDef:
    def __init__(self, position, direction, energy, time, propagated_distance, depth, ID):
        self.position = position
        self.direction = direction
        self.energy = energy
        self.time = time
        self.propagated_distance = propagated_distance
        self.depth = depth
        self.ID = ID


def start_shower(energy, max_depth, direction, position):
    time = 0
    propagated_distance = 0
    
    particle_list = []
    particle_list.append(ParticleDef(position, direction, energy, time, propagated_distance, 0, Particle.PHOTON))    
    
    list_position_i = []
    list_position_f = []
    list_ID = []
    list_energy_i = []

    particle_num = 0

    while(len(particle_list)>0):
        CP = particle_list.pop() # get last item and remove it from list
        
        ### Stop propagation if we are in too deep
        if(CP.depth>=max_depth):
            continue
        
	#print(CP.ID)
	#print(CP.energy)
	#print(CP.position)
        #print(CP.depth)
	#print(CP.direction)
	if(CP.ID == Particle.PHOTON or CP.ID == Particle.BREMS):
            secondaries, end_position = propagate_photon(CP.position, CP.direction, CP.energy, CP.time, CP.propagated_distance)
        elif(CP.ID == Particle.ELECTRON):
            secondaries, end_position = propagate_electron(CP.position, CP.direction, CP.energy, CP.time, CP.propagated_distance)
        elif(CP.ID == Particle.POSITRON):
            secondaries, end_position = propagate_positron(CP.position, CP.direction, CP.energy, CP.time, CP.propagated_distance)
        else:
            print("This should never happen")
            

        position_previous = CP.position

        for sec in secondaries:
            if(sec.id == pp.particle.Data.Particle):
		particle_num = particle_num + 1
                if(sec.particle_def.name == "EPlus"):
                    particle_list.append(ParticleDef(sec.position, sec.direction, sec.energy, sec.time, sec.propagated_distance, CP.depth+1, Particle.POSITRON))
                elif(sec.particle_def.name == "EMinus"):
                    particle_list.append(ParticleDef(sec.position, sec.direction, sec.energy, sec.time, sec.propagated_distance, CP.depth+1, Particle.ELECTRON))
                elif(sec.particle_def.name == "Gamma"):
                    particle_list.append(ParticleDef(sec.position, sec.direction, sec.energy, sec.time, sec.propagated_distance, CP.depth+1, Particle.PHOTON))
            elif(sec.id == pp.particle.Data.Brems):
		    particle_num = particle_num + 1
                    particle_list.append(ParticleDef(sec.position, sec.direction, sec.energy, sec.time, sec.propagated_distance, CP.depth+1, Particle.BREMS))

            ##save information here
            list_ID.append(CP.ID)
            list_energy_i.append(sec.parent_particle_energy)
            list_position_i.append(position_previous)
            list_position_f.append(sec.position)
            
            position_previous = sec.position


    print("Shower completed, " + str(particle_num) + " particles have been propagated")
    return (list_position_i, list_position_f, list_ID, list_energy_i)


if(len(sys.argv) !=3 and len(sys.argv)!=4):
    print("Invalid arguments! First argument must be energy, second argument shower depth. Third argument optional filename")
    sys.exit()
print("Start shower propagation with initial energy " + str(sys.argv[1] + " MeV and maximal depth of " + str(sys.argv[2]) + " generations."))

list_position_i, list_position_f, list_ID, list_energy_i = start_shower(float(sys.argv[1]), int(sys.argv[2]), direction=pp.Vector3D(0, 0, -1), position=pp.Vector3D(0, 0, 1000000))

x_list_i = []
y_list_i = []
z_list_i = []

x_list_f = []
y_list_f = []
z_list_f = []

for positions in list_position_i:
    x_list_i.append(positions.x)
    y_list_i.append(positions.y)
    z_list_i.append(positions.z)

for positions in list_position_f:
    x_list_f.append(positions.x)
    y_list_f.append(positions.y)
    z_list_f.append(positions.z)

### workaround because there is an old version of the Enum package on vollmond...

list_ID_tmp = []

if(type(list_ID[0]) != int):
    for element in list_ID:
        list_ID_tmp.append(element.value)
    list_ID = list_ID_tmp

# save data to txt

if(len(sys.argv)==4):
    np.savetxt(sys.argv[3], np.column_stack([x_list_i, y_list_i, z_list_i, x_list_f, y_list_f, z_list_f, list_ID, list_energy_i]), header="x_i y_i z_i x_f y_f z_f ID E_i")
else:
    np.savetxt('data_' + sys.argv[1] + '_' + sys.argv[2] + '.txt', np.column_stack([x_list_i, y_list_i, z_list_i, x_list_f, y_list_f, z_list_f, list_ID, list_energy_i]), header="x_i y_i z_i x_f y_f z_f ID E_i")





