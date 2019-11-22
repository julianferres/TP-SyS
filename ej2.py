# Ejercicio 2
from ecg import *

yQRS = y[137:160]
xQRS = x[137:160]

plt.plot(xQRS,yQRS)
plt.show()

unQRS = y[312:329]
fft = fft(unQRS,2048)
freq = np.linspace(-math.pi,math.pi,len(fft))
response = 20* np.log10(np.abs(fftshift(fft/abs(fft).max())))

sinDB = abs(fftshift(fft))

plt.plot(freq, sinDB)
plt.title("Espectro de frecuencias del complejo QRS")

plt.show()

plt.figure()
