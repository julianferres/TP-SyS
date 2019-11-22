from ecg import y
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft
from scipy.fftpack import fftshift

fs = 200

'''
Espectrograma completo
'''
window = signal.tukey(3200) # ventana de Tukey de 3200 muestras
f, t, Sxx = signal.spectrogram(y, fs, window)
plt.pcolormesh(t,f,np.log10(Sxx))
plt.title("Espectrograma completo con ventana de 3200 muestras.")
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()


# 5 Latidos, con diferentes ventanas
window = signal.tukey(55) # ventana de Tukey de 3200 muestras
f, t, Sxx = signal.spectrogram(y[100:1000], fs, window)
plt.pcolormesh(t,f,Sxx)
plt.title("Espectrograma de 5 latidos con ventana de 55 muestas.")
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

# 5 Latidos, con diferentes ventanas
window = signal.tukey(256) # ventana de Tukey de 3200 muestras
f, t, Sxx = signal.spectrogram(y[100:1000], fs, window)
plt.pcolormesh(t,f,Sxx)
plt.title("Espectrograma de 5 latidos con ventana de 256 muestas.")
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

