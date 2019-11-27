#%% Imports y funciones auxiliares
from ecg import *
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft
from scipy.fftpack import fftshift


"""Funciones auxiliares"""
def polosYceros (num,den, nombre, archivo,carpeta):
    z = np.roots(num)
    p = np.roots(den)
    plt.plot(np.real(z), np.imag(z), 'ob', markerfacecolor='None',markersize=10)
    plt.plot(np.real(p), np.imag(p), 'xr',markersize=10)
    plt.grid(which='both')
    plt.xlabel("Re(Z)")
    plt.ylabel("Im(Z)")
    plt.plot(np.cos(np.r_[0:6.5:0.1]), np.sin(np.r_[0:6.5:0.1]), color='0.2')  # circulo unitario
    plt.title('Polos y ceros de %s'%nombre)
    plt.axis('equal')
    plt.show();

def rtaFreq (num, den, nombre, archivo, carpeta):
    
    [w, H] = signal.freqz(num, den)
    plt.subplot(2,1,1)
    plt.plot(w / np.pi, 20 * np.log10(np.abs(H)))
    plt.ylabel('Modulo [dB]')
    plt.grid(b=1, which='both', axis='both', linewidth=0.5)
    plt.title('Respuesta en frecuencia de %s'%nombre)
    plt.subplot(2,1,2)
    plt.plot(w / np.pi, np.angle(H) / np.pi)
    plt.ylabel('Fase [rad/pi]') # fase normalizada por sobre pi 1---pi lo mismo en frecuencia
    plt.grid(b=1, which='both', axis='both', linewidth=0.5)
    plt.xlabel('Frecuencia normalizada x/pi')
    plt.show();
    plt.close()

def ver_ancho_de_banda(num,den):
    [w, H] = signal.freqz(num, den)
    H = 20 * np.log10(np.abs(H))
    print(max(H))
    l,r = -1,-1 
    for i in range(len(w)):
        if(l==-1 and H[i]>-3): #Cuando subo de -3dB, tengo el lado izq
            l = i
        if(l!=-1 and H[i]<-3): #Cuando vuelvo a bajar de -3dB, y ya tengo l, tengo r
            r = i
    l = l/len(w) * fs/2; r = r/len(w)*fs/2
    return l,r


def rtaImpulso (num, den, nombre, archivo, cant, carpeta, minx):
    x_in = np.zeros(cant+1);
    x_in[0] = 1
    y_out = signal.lfilter(num, den, x_in)
    print(y_out)
    n = np.r_[0:cant+1]

    fig = plt.figure(1)

    ax = fig.add_subplot(1, 1, 1)

    x_major_ticks = np.arange(0,cant+1,2)
    y_major_ticks = np.arange(0,y_out.max() + 1,minx)
    
    ax.set_xticks(x_major_ticks)
    
    ax.set_yticks(y_major_ticks)

    # And a corresponding grid
    ax.grid(which='both')

    print(y_out)
    plt.stem(n, y_out, 'r')
    plt.title('Respuesta al impulso de %s'%nombre)
    plt.xlim([-1 ,cant+2])
    plt.xlabel('n')
    plt.ylabel('h (n)')
    plt.show();
    plt.close()

def rtaEscalon(num, den, nombre, archivo,cantidad,carpeta):
    
    x_in=np.ones(cantidad)
    y_out = signal.lfilter(num, den ,x_in)
    n = np.r_[0:cantidad]

    fig = plt.figure(1)

    ax = fig.add_subplot(1, 1, 1)

    x_major_ticks = np.arange(-1, cantidad+1, 2)
    y_major_ticks = np.arange(0,41,2)

    ax.set_xticks(x_major_ticks)
    ax.set_yticks(y_major_ticks)

    # And a corresponding grid
    ax.grid(which='both')


    plt.stem(n, y_out, 'r')
    plt.xlabel('n')
    plt.ylabel('s (n)')
    plt.title('Respuesta al escalon de %s'%nombre)
    plt.show();

filter_signal = lambda num,den,entrada: signal.lfilter(num,den,entrada)

compensar_retardo = lambda N,entrada : np.concatenate((entrada[N:len(entrada)], np.zeros(N)))
# CORRO A LA IZQUIERDA N


def qrs_detection(fm):
    
    candidato_int = np.zeros(fm.size)
    noise_peak_int = np.zeros(fm.size)
    n = 0
    esperar = 0
    rr_avg_2_int = 0

    while n <= fm.size:
        if n == 720:
            aux = fm[0:n]
            thr_signal_int = aux.max()
            sig_level_int = thr_signal_int
            aux = aux[aux.nonzero()]
            thr_noise_int = min(aux)
            noise_level_int = thr_noise_int

        if n > 720 and esperar == 0:
            if fm[n-2] > thr_signal_int:
                candidato_int[n-2] = fm[n-2]
                esperar = 72

                aux, = np.where(candidato_int)
                aux_diff = np.diff(aux)
                principio = max(0, aux_diff.size-8-1)
                if aux_diff.size:
                    rr_avg_1_int = aux_diff[principio:aux_diff.size].mean()
                else:
                    rr_avg_1_int = np.nan

                if rr_avg_2_int == 0:
                    rr_avg_2_int = rr_avg_1_int
                elif aux_diff.size > 8:
                    aux_diff_2, = np.where(np.logical_and(aux_diff > 0.92*rr_avg_2_int, aux_diff < 1.16*rr_avg_2_int))
                    if aux_diff_2.size:
                        principio2 = max(0, aux_diff_2.size-1-8)
                        rr_avg_2_int = aux_diff_2[principio2:aux_diff_2.size].mean()

                rr_missed_limit = rr_avg_2_int*1.66

                if aux_diff.size > 8 and aux_diff[aux_diff.size-1] > rr_missed_limit:
                    fm_aux = fm[(aux[aux.size-2]):(aux[aux.size-1])]
                    aux_fm_aux, = aux[aux.size-2] + 1 + np.where(fm_aux > thr_noise_int)

                    if np.size(aux_fm_aux):
                        candidato_int[n-2] = 0
                        candidato_int[aux_fm_aux[1]-1] = fm[aux_fm_aux[1]-1]
                        n = aux_fm_aux[0] + 1
                        sig_level_int = 0.25*fm[n-2] + 0.75*sig_level_int
                else:
                    sig_level_int = 0.125*fm[n-2] + 0.875*sig_level_int

            else:
                noise_level_int = 0.125*fm[n-2] + 0.875*noise_level_int
                noise_peak_int[n - 2] = 1
            thr_signal_int = noise_level_int + 0.25 *(sig_level_int-noise_level_int)
            thr_noise_int = 0.5*thr_signal_int
        elif esperar > 0:
            esperar = esperar - 1
        n = n + 1

    return candidato_int


def desplazar_marcas (marcas_propias,y,rango):
    marcas_desplazadas  = np.zeros(len(marcas_propias))
    for i in range(len(marcas_propias)-rango):
        if(marcas_propias[i]!=0):
            idxMaxLocal, actMax = -1, -float('inf')
            for delta in range(rango): #Tengo mi maximo local para cada marca
                if(y[i+delta]>actMax):
                    idxMaxLocal = i+delta
                    actMax = y[i+delta]
            marcas_desplazadas[idxMaxLocal] = marcas_propias[i]
    return marcas_desplazadas
    

#%%EJERCICIO 5:

"""DEFINO HL (W)"""
num_hl = np.zeros(13)
den_hl = np.zeros(13)

num_hl[0], num_hl[6], num_hl[12] = 1,-2,1
#Menos el signo

den_hl[0], den_hl[1], den_hl[2] = 1, -2, 1


"""Diagrama de polos y ceros de HL"""
polosYceros(num_hl, den_hl,'HL','ej5_diag_pyz_hl','ej5/HL')

"""Respuesta en frecuencia de HL"""
rtaFreq(num_hl,den_hl,'HL(w)','ej5_rta_freq_hl','ej5/HL')

"""RTA AL IMPULSO HL"""
rtaImpulso(num_hl, den_hl,'HL','ej5_rta_impulso_hl',18,'ej5/HL',1)

"""RTA AL ESCALON HL"""
rtaEscalon(num_hl, den_hl, 'HL', 'ej5_rta_escalon_hl',18,'ej5/HL')

"""RETARDO DE HL"""

w, gd = signal.group_delay((num_hl,den_hl))
plt.plot(w/np.pi,gd)
plt.title("Retardo producido por HL")
plt.xlabel ('Pulsacion normalizada/pi ')
plt.ylabel('Retardo en cantidad de muestras')
plt.show();


#%% EJERCICIO 6
"""DEFINO HH (W)"""
num_hh = np.zeros(33, dtype=float)
den_hh = np.zeros(33, dtype=float)

num_hh[0], num_hh[16], num_hh[17], num_hh[32] = -1/32, 1, -1, 1/32
 # MENOS SIGNIFICATIVO

den_hh[0], den_hh[1] = 1, -1


"""DIAG POLOS Y ZEROS HH"""
polosYceros(num_hh, den_hh,'HH','ej5_diag_pyz_hh','ej5/HH')

"""RTA EN FREQ HH"""
rtaFreq(num_hh,den_hh,'HH','ej5_rta_freq_hh','ej5/HH')

"""RTA AL IMPULSO HH"""
rtaImpulso(num_hh,den_hh,'HH','ej5_rta_impulso_hh',40,'ej5/HH',0.2)

"""RTA AL ESCALON HH"""
rtaEscalon(num_hh, den_hh, 'HH', 'ej5_rta_escalon_hh',40,'ej5/HH')

"""RETARDO DE HH"""

w, gd = signal.group_delay((num_hh,den_hh)) # 16 muestras
plt.plot(w/np.pi,gd)
plt.grid(which='both')
plt.title("Retardo producido por HH")
plt.xlabel('Pulsacion normalizada/pi ')
plt.ylabel('Retardo en cantidad de muestras')
plt.show();


#################################################

"""EFECTO SOBRE LA SEÑAL"""
fs = 200;
y = (y - y.mean()) / 200 # ganancia 200

y_out_hl = filter_signal(num_hl,den_hl,y)
y_out_hl = compensar_retardo(5,y_out_hl) 
y_out_hh = filter_signal(num_hh,den_hh,y_out_hl) 
y_out_hh = compensar_retardo(16,y_out_hh)

fig=plt.figure()
ax1=fig.add_subplot(3,1,1)
plt.plot(t, y)
plt.xticks(fontsize=0)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.yticks(np.arange(-1,1.5,0.5))
plt.grid(which='both')

plt.xlim(2000,4000)
plt.ylim(-1,1.5)


ax2=fig.add_subplot(3,1,2)
plt.plot(t,y_out_hl)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.xticks(fontsize=0)
plt.yticks(np.arange(-20, 25, 5))
plt.grid(which='both')

plt.ylabel('Amplitud [mV]')
plt.xlim(2000, 4000)
plt.ylim(-15, 20)

""""""""""""

plt.subplot(3,1,3)
plt.plot(t, y_out_hh)
plt.xlabel('Tiempo [ms]')
plt.yticks(np.arange(-20, 25, 5))
plt.grid(which='both')
plt.xlim(2000, 4000)
plt.ylim(-15, 22)
plt.show();


"""TRANSFORMADAS DE CADA UNO"""

fig = plt.figure()
n = 1024
fs = 200

A = fft(y[137:160],n)
B = fft(y_out_hl[137:160],n)
C = fft(y_out_hh[137:160],n)

freq = np.linspace(-fs/2, fs/2, len(A)) # -fs/2 fs/2 y dsp pueden hacer xlim ([0, fs/2])

response = 20 * np.log10(np.abs(fftshift(A)))
response_hl = 20 * np.log10(np.abs(fftshift(B)))
response_hh = 20 * np.log10(np.abs(fftshift(C)))

ax1=fig.add_subplot(3,1,1)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.plot(freq, response)
plt.xlim([0,fs/2])
plt.title("Carácterísticas en frecuencia QRS:Filtrado pasa-banda")

plt.xticks(np.arange(0, 101, 5))
plt.yticks(np.arange(-100,80,10))
plt.axis([0, 100, -65, 20])
plt.grid(which='both')
plt.show();


ax2 = fig.add_subplot(3, 1, 2)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.plot(freq, response_hl)
plt.xlim([0, fs / 2])
plt.ylabel("Modulo [dB]")
plt.xticks(np.arange(0, 101, 5))
plt.yticks(np.arange(-100, 80, 10))
plt.axis([0, 100, -20, 50])
plt.grid(which='both')

ax3 = fig.add_subplot(3, 1, 3)
plt.plot(freq, response_hh)
plt.xlim([0, fs / 2])
plt.xlabel("Frecuencia [Hz]")
plt.xticks(np.arange(0, 101, 5))
plt.yticks(np.arange(-100, 80, 10))
plt.axis([0, 100, -20, 50])
plt.grid(which='both')
plt.close()


#%%EJERCICIO 8

time_sample = 1/fs
h = 1

num_hd = np.zeros(5)
den_hd = 1

num_hd[0] = 1 / (8*h)
num_hd[1] = 1 / (4*h)
num_hd[3] = -1 / (4 * h)
num_hd[4] = -1 / (8*h) #MENOS SIGNIFICATIVO

rtaFreq(num_hd,den_hd,'HD','ej8rta_freq_hd','ej8/HD')


y_out_hd = filter_signal(num_hd,den_hd,y_out_hh)

print(ver_ancho_de_banda(num_hd,den_hd))

"""LINEAL DE HD"""
[w, H] = signal.freqz(num_hd, den_hd)
plt.subplot(2, 1, 1)
plt.plot(w / np.pi, np.abs(H))
plt.ylabel('Modulo')
plt.xticks(np.arange(0, 1.1, 0.1))
plt.grid(b=1, which='both', axis='both', linewidth=0.5)
plt.title('Respuesta en frecuencia de HD LINEAL')

plt.subplot(2, 1, 2)
plt.plot(w / np.pi, np.angle(H) / np.pi)
plt.ylabel('Fase [rad/pi]')  # fase normalizada por sobre pi 1<->pi lo mismo en frecuencia
plt.xticks(np.arange(0, 1.1, 0.1))
plt.grid(b=1, which='both', axis='both', linewidth=0.5)
plt.xlabel('Frecuencia normalizada x/pi')

plt.show();
plt.close()


"""RETARDO DE HL"""

w, gd = signal.group_delay((num_hd, den_hd))
plt.plot(w/np.pi,gd)
plt.title("Retardo producido por HD")
plt.xlabel ('Pulsacion normalizada/pi ')
plt.ylabel('Retardo en cantidad de muestras')
plt.show();

y_out_hd = compensar_retardo(2,y_out_hd)

plt.figure()

plt.plot(t, y_out_hd)
plt.title('Señal filtrada con HD')
plt.xlabel('Tiempo [ms]')
plt.ylabel('Amplitud [mV]')
plt.xlim(2000, 4000)
plt.ylim(-6, 6)
plt.grid(which='both')
plt.show();

plt.close()




# %% EJERCICIO 9 (Elevar al cuadrado)

y_out_cuad = np.power(y_out_hd,2)

plt.figure()
plt.plot(t,y_out_cuad)
plt.title('Señal al cuadrado')
plt.xlabel('Tiempo [ms]')
plt.ylabel('Amplitud [mV]')
plt.xlim(2000, 4000)
plt.ylim(-1, 35)
plt.grid(which='both')
plt.show()


# Caracterización en frecuencia de la señal elevada al cuadrado
fig = plt.figure()
n = 1024 #Cantidad de puntos de la FFT

freqCuad = fft(y_out_cuad,n)

freq = np.linspace(-fs/2, fs/2, len(freqCuad)) # -fs/2 fs/2 y dsp pueden hacer xlim ([0, fs/2])
response = 20 * np.log10(np.abs(fftshift(freqCuad)))
plt.title("Espectro de frecuencias de la señal al cuadrado")
plt.ylabel("Modulo [dB]")
plt.xlabel("Frecuencia [Hz]")

plt.plot(freq, response)
plt.show();



# %% EJERCICIO 10 (Integración) 1/N *(1-z^-N)/(1-z^-1)

N = 20 #PROBAR 20

num_hi = np.zeros(N+1)
den_hi = np.zeros (2)

num_hi[0], num_hi[N] = 1, -1

den_hi[0], den_hi[1] = N, -N

rtaFreq(num_hi,den_hi,'HI','ej10rta_freq_hi','ej10/HI') # GROUP DELAY DE 15 MUESTRAS

y_out_hi = filter_signal(num_hi,den_hi,y_out_cuad)


"""RETARDO DE HI"""

w, gd = signal.group_delay((num_hi,den_hi))
plt.plot(w/np.pi,gd)
plt.title("Retardo producido por HL")
plt.xlabel ('Pulsacion normalizada/pi ')
plt.ylabel('Retardo en cantidad de muestras')
plt.show();

rtaImpulso(num_hi, den_hi,'HI','ej5_rta_impulso_hl',22,'ej5/HL',1)

y_out_hi = compensar_retardo( N//2, y_out_hi)

plt.plot(t, y_out_hi)
plt.title('Señal integrada')
plt.xlabel('Tiempo [ms]')
plt.ylabel('Amplitud [mV]')
plt.xlim(2000, 4000)
plt.ylim(-1, 10)
plt.grid(which='both')
plt.show();


plt.plot(t, y_out_hi)
plt.title('Señal integrada')
plt.xlabel('Tiempo [ms]')
plt.ylabel('Amplitud [mV]')
plt.xlim(2000, 4000)
plt.ylim(-1, 20)
plt.grid(which='both')
plt.subplot()
plt.plot(t, y_out_cuad,'r') # PARA VER SUPERPUESTA EL CUAD
plt.legend(["Señal integrada","Señal al cuadrado"])
plt.title('Señal integrada y señal al cuadrado')
plt.show();




#%% EJERCICIO 12 (Gold Standard algorithm)

marcas_propias = qrs_detection(y_out_hi)

marcas_desplazadas = desplazar_marcas(marcas_propias,y,30)
plt.figure()
plt.plot (t,marcas_propias,'go', markersize = 5)
plt.plot (t,marcas_desplazadas,'rx', markersize = 5)
plt.ylim(-1,1.5)

plt.subplot()

plt.plot (t,y)
plt.xlim(20000, 22000)
plt.show()

