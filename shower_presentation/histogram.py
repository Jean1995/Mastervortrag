import numpy as np
import matplotlib.pyplot as plt
import sys


def binning(upper, lower, x_up, x_low):
    if(upper < lower):
        upper, lower = lower, upper
    condition =  (x_up >= lower) & (x_low <= lower) | (x_up >= upper) & (x_low <= upper) | (x_up <= upper) & (x_low >= lower) | (x_up >= upper) & (x_low <= lower) 
    return np.sum(condition)

x, y, z, x_pre, y_pre, z_pre, ID, E_i = np.genfromtxt(sys.argv[1], unpack=True)

if(len(sys.argv)<=2):
	bin_num = 100
else:
	bin_num = sys.argv[2]

## for air from http://www-ekp.physik.uni-karlsruhe.de/~jwagner/WS0809/Vorlesung/TP-WS08-6.pdf
E_c = 84 # MeV
X_0 = 36.62 # g / cm^2
rho = 1.205 * 1e-3 # g/cm^3
X = X_0 / rho


if(len(sys.argv)>=3):
	if(sys.argv[3] != 0):
		X_max = X * (np.log(float(sys.argv[3])) - np.log(E_c)) / np.log(2)
		print(X_max)
		plt.axvline(x = 1000000 - X_max, label='Schauermaximum (Heitler)', color='k', linestyle='dashed', linewidth=1)

		

binning=np.vectorize(binning, excluded = ['x_up', 'x_low'])
bin_borders = np.linspace(min( min(z), min(z_pre) ), max( max(z), max(z_pre) ), bin_num)

bin_vals = binning(bin_borders[:-1], bin_borders[1:], x_up = z_pre, x_low = z)


plt.step(bin_borders[:-1], bin_vals, where='post')
plt.ylabel("Teilchenanzahl")
plt.xlabel("Abstand vom Erdboden / cm")
plt.xlim(0, 1000000)
plt.legend(loc = 'best')
plt.grid()

if(len(sys.argv)<=4):
	plt.show()
else:
	plt.savefig(sys.argv[4] + '.pdf')
