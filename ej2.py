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
freq = 100/np.pi * np.linspace(-math.pi,math.pi,len(fft))
#Lo dejo en Hz (fs/2 = pi)
response = 20* np.log10(np.abs(fftshift(fft/abs(fft).max())))

sinDB = abs(fftshift(fft))

#Sin dB
plt.plot(freq, sinDB)
plt.title("Espectro de frecuencias de un complejo QRS")
plt.xlabel("Frecuencia [Hz]")
plt.ylabel("Modulo")
plt.xlim(0,100)
plt.savefig('graficos/QRSfrecuencia.png',dpi=300)
plt.show()

#Con dB
plt.plot(freq, response)
plt.title("Espectro de frecuencias de un complejo QRS en dB")
plt.ylabel("Modulo [dB]")
plt.xlabel("Frecuencia [Hz]")
plt.xlim(0,100)

plt.savefig('graficos/QRSfrecuenciadB.png',dpi=300)
plt.show()
