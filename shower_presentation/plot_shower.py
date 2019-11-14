import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib
import matplotlib.cm as cm
from matplotlib.lines import Line2D

fig = plt.figure(figsize=plt.figaspect(1.5)*3.0)
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')
ax.view_init(0, 45)

x_i_list, y_i_list, z_i_list, x_f_list, y_f_list, z_f_list, ID_list, energy_i_list = np.genfromtxt(sys.argv[1], unpack=True)

count_electron = 0
count_positron = 0
count_photon = 0

for x_i, y_i, z_i, x_f, y_f, z_f, ID, energy_i in zip(x_i_list, y_i_list, z_i_list, x_f_list, y_f_list, z_f_list, ID_list, energy_i_list):

	if(energy_i < 100):
		continue

	if(ID == 1):
		ax.plot([x_i, x_f], [y_i, y_f], [z_i ,z_f], c='g', linewidth=0.15)
		count_electron+=1
	elif(ID == 0 or ID==3):
		ax.plot([x_i, x_f], [y_i, y_f], [z_i ,z_f], c='b', linewidth=0.15)
		count_photon+=1
	elif(ID == 2):
		ax.plot([x_i, x_f], [y_i, y_f], [z_i ,z_f], c='r', linewidth=0.15)
		count_positron+=1
	else:
		print("Unknown particle_id")


## custom legend

print("Showing " + str(count_photon) + " Photons, " + str(count_electron) + " electrons and " + str(count_positron) + " positrons")

custom_lines = [Line2D([0], [0], color='b', lw=2),
                Line2D([0], [0], color='g', lw=2),
                Line2D([0], [0], color='r', lw=2)]

ax.legend(custom_lines, ['Photon', 'Elektron', 'Positron'], fontsize=15, loc='best')

ax.set_xlabel('X \ cm')
ax.set_ylabel('Y \ cm')
ax.set_zlabel('Z \ cm')
ax.set_xlim(-100000, 100000)
ax.set_ylim(-100000, 100000)
ax.set_zlim(0, 1000000)

### Enable for propaganda-mode
#plt.axis('off')
#ax.grid(False)
#plt.savefig('demo1.png', transparent=True, dpi=2500)

if(len(sys.argv)==2):
	plt.show()
else:
	plt.savefig(sys.argv[2] + '.pdf')

