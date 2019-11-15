import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import matplotlib
import matplotlib.cm as cm
from matplotlib.lines import Line2D


x_i_list, y_i_list, z_i_list, x_f_list, y_f_list, z_f_list, ID_list, energy_i_list = np.genfromtxt(sys.argv[1], unpack=True)

vec_i = np.array([x_i_list, y_i_list, z_i_list]).T
vec_f = np.array([x_f_list, y_f_list, z_f_list]).T

## choose projection plane here
projection_vec_1 = np.array([1,0,0])
projection_vec_2 = np.array([0,1,0])

proj_1_i_list = np.matmul(vec_i, projection_vec_1)
proj_2_i_list = np.matmul(vec_i, projection_vec_2)

proj_1_f_list = np.matmul(vec_f, projection_vec_1)
proj_2_f_list = np.matmul(vec_f, projection_vec_2)





count_electron = 0
count_positron = 0
count_photon = 0

plt.figure(figsize=(5, 5))

for proj_1_i, proj_2_i, proj_1_f, proj_2_f, ID, energy_i in zip(proj_1_i_list, proj_2_i_list, proj_1_f_list, proj_2_f_list, ID_list, energy_i_list):

	if(energy_i < 100):
		continue

	if(ID == 1):
		plt.plot([proj_1_i,proj_1_f], [proj_2_i,proj_2_f], 'g-', linewidth=0.15)
		count_electron+=1
	elif(ID == 0 or ID==3):
		plt.plot([proj_1_i,proj_1_f], [proj_2_i,proj_2_f], 'b-', linewidth=0.15)
		count_photon+=1
	elif(ID == 2):
		plt.plot([proj_1_i,proj_1_f], [proj_2_i,proj_2_f], 'r-', linewidth=0.15)
		count_positron+=1
	else:
		print("Unknown particle_id")


## custom legend

print("Showing " + str(count_photon) + " Photons, " + str(count_electron) + " electrons and " + str(count_positron) + " positrons")

custom_lines = [Line2D([0], [0], color='b', lw=2),
                Line2D([0], [0], color='g', lw=2),
                Line2D([0], [0], color='r', lw=2)]

plt.legend(custom_lines, ['Photon', 'Elektron', 'Positron'], fontsize=15, loc='best')

plt.grid(True)
plt.xlabel('X \ cm', fontsize=15)
plt.ylabel('Y \ cm', fontsize=15)
plt.xlim(-50000, 50000)
plt.ylim(-50000, 50000)
plt.tight_layout()

### Enable for propaganda-mode
#plt.axis('off')
#ax.grid(False)
#plt.savefig('demo1.png', transparent=True, dpi=2500)

if(len(sys.argv)==2):
	plt.show()
else:
	plt.savefig(sys.argv[2] + '.pdf')

