# -*- coding: utf-8 -*-
"""
SyS - Curso 1 - FIUBA
Comparamos 3 formas diferentes de realizar filtrado.
Se ejemplifica el uso de:
  - lfilter (FIR e IIR)
  - conv
  - fft
Finalmente, se analiza el filtro IIR utilizado.
"""

from numpy import * # evita poner 'np.' cada vez
from scipy.signal import lfilter, freqz
from scipy.fftpack import fft, ifft
from matplotlib.pyplot import *

#%% Filtrado FIR

# Definimos senial x (una delta desplazada) y filtro h (FIR)
x = zeros(15); x[3] = 1
h = r_[1:6]  # ojo: r_[1:6] genera vector [1,2,3,4,5]
n = r_[0:15] # eje temporal (solo para graficos)

# Probamos diferentes formas de filtrado
y1 = lfilter(h,1,x)
y2 = convolve(h,x); y2 = y2[0:-(len(h)-1)] # cortamos transitorio final
Nfft = len(x) + len(h) - 1 # esto completa con ceros las seniales al hacer fft
y3 = ifft(fft(x,Nfft)*fft(h,Nfft)); y3 = y3[0:-(len(h)-1)] # cortamos transitorio final
y3 = real(y3) # eliminamos posible residuo imaginario numerico

# Graficamos y verificamos que son todas iguales
subplot(411),stem(n,x,'r'),title('Senial x(n)'),ylim(0,5)
subplot(412),stem(n,y1),title('y1 usando lfilter'),ylim(0,5)
subplot(413),stem(n,y2),title('y2 usando conv'),ylim(0,5)
subplot(414),stem(n,y3),title('y3 usando fft'),ylim(0,5)
xlabel('n'),tight_layout(rect=(0,0,1,2)),show()

#%% Filtrado IIR

# En lugar de definir un filtro FIR mediante un vector, lo definimos
# mediante 2 vectores los cuales guardan los coeficientes del numerador y
# denominador de la ecuacion en diferencias. Implementamos:
#   y(n) - 0.3*y(n-1) = 0.5*x(n) + x(n-1)

# Extraemos coeficientes de la ec. en diferencias
b = array([0.5,1])  # coef del numerador
a = array([1,-0.3]) # coef del denominador

# Definimos senial a filtrar x2 (una delta) y filtramos. Para filtrar IIR
# se utiliza lfilter (no sirven conv ni fft a menos que aproximemos h(n))
x2 = zeros(20); x2[0] = 1
y4 = lfilter(b,a,x2)
n = r_[0:20]

# Graficamos. Notar que y4 es la rta al impulso del sistema dado por la ec
# en diferencias (es IIR pero se extingue exponencialmente)
stem(n,y4,'r'),title('y4 usando lfilter IIR'),show()

#%% Analisis de filtro IIR

# Polos y ceros de H(z)
z = roots(b)
p = roots(a)
plot(real(z),imag(z),'ob',markerfacecolor='None')
plot(real(p),imag(p),'xb'),grid()
plot(cos(r_[0:6.5:0.1]),sin(r_[0:6.5:0.1]),color='0.7') # circulo unitario
title('Polos y ceros de H(w)'),show()

# Magnitud y fase de H(Omega)
[w,H] = freqz(b,a)
subplot(211),plot(w/pi,20*log10(abs(H))),ylabel('Modulo [dB]'),grid()
title('Respuesta en frecuencia de H(w)')
subplot(212),plot(w/pi,unwrap(angle(H))),ylabel('Fase [rad]'),grid()
xlabel('Frecuencia normalizada x\pi')
show()

# Nota: para filtros FIR, b es el filtro y poner a=1 (coef del denominador)



