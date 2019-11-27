#%%
# Ejercicio 2
from ecg import *

yQRS = y[137:160]
xQRS = x[137:160]

plt.plot(xQRS,yQRS)
plt.title("QRS de un ciclo")
plt.xlabel("Numero de muestra")
plt.savefig('graficos/un_cicloQRS.png',dpi=300)
plt.show()

unQRS = y[312:329]
fft = fft(unQRS,2048)
freq = np.linspace(-math.pi,math.pi,len(fft))
response = 20* np.log10(np.abs(fftshift(fft/abs(fft).max())))

sinDB = abs(fftshift(fft))

#Sin dB
plt.plot(freq, sinDB)
plt.title("Espectro de frecuencias de un complejo QRS")
plt.ylabel("Modulo")
plt.savefig('graficos/QRSfrecuencia.png',dpi=300)
plt.show()

#Con dB
plt.plot(freq, response, xunits = radians)
plt.title("Espectro de frecuencias de un complejo QRS en dB")
plt.ylabel("Modulo [dB]")
plt.savefig('graficos/QRSfrecuenciadB.png',dpi=300)
plt.show()

plt.figure()
