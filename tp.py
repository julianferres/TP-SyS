F = 200
import scipy.io as scp
import numpy as np
import matplotlib.pyplot as plt

datamat = scp.loadmat('Datos/103m.mat')
datamat_y = datamat['y']
y = np.array(datamat_y)[0]
y-=y.mean()

print(y)
marcas = np.loadtxt('Datos/marcas.txt')
print(len(marcas), marcas)
plt.plot(y)
plt.show()
x = np.arange(len(y))

#Me pasa a ms desde el num de muestra
pasarAms = lambda x: x*5

"""

"""

#1
"""
Hacerlo en draw.io para un periodo
"""

#2
yQRS = y[137:160]
xQRS = x[137:160]

plt.plot(xQRS,yQRS)
plt.show()

unQRS = y[312:329]
fft = fft(unQRS,2048)
freq = np.linspace(-math.pi,math.pi,len(fft))
response = 20* np.log10(np.abs(fftshift(fft/abs(fft).max())))

sinDB = abs(fftshift(fft))

plt.show(freq, sinDB)

plt.show()

plt.figure()

plt.plot(unQRS)

plt.show()

#Saber el uso de filter



