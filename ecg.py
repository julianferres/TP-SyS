F = 200
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

marcas = np.loadtxt('Datos/marcas.txt')
if __name__ == '__main__':
	plt.plot(y)
	plt.title("Muesta completa")
	plt.show()
	plt.plot(unCiclo)
	plt.title("Ciclo Completo")
	plt.show()
x = np.arange(len(y))



