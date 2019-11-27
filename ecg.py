#%%
F = 200 #Frecuencia de sampling
ts = 1/200 #Periodo de sampling 
import scipy.io as scp
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.fftpack import fftshift
import math

datamat = scp.loadmat('Datos/103m.mat')
datamat_y = datamat['y']
y = np.array(datamat_y)[0]
y-=y.mean()
unCiclo = y[100:280]

t = np.arange(len(y))
t = t*ts*1000.0

marcas = np.loadtxt('Datos/marcas.txt')
if __name__ == '__main__':
	plt.plot(y)
	plt.title("Muesta completa")
	plt.xlabel("Numero de muestra")
	plt.ylabel("Amplitud")
	plt.savefig('graficos/muestra_completa.png',dpi=300)
	plt.show()
	plt.plot(y)
	plt.xlim([100 ,280])
	plt.ylim([-125,300])
    
	plt.title("Ciclo Completo")
	plt.xlabel("Numero de muestra")
	plt.ylabel("Amplitud")
	plt.savefig('graficos/un_ciclo.png',dpi=300)
	plt.show()

x = np.arange(len(y))
