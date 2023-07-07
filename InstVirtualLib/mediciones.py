# -*- coding: utf-8 -*-
"""
@author: Pablo, Ramiro

Este Módulo contiene la biblioteca de mediciones. Todos los procedimientos de 
medición que se deseen automatizar se implementaran en esta clase.

La idea es que esta clase tome como entrada los vectores (tension, fase, etc)
y calcule los valores solicitados.

Todos los calculos de los 

"""

from math import floor
import numpy as np
import scipy.signal as dsp

class Mediciones():

    def __init__(self):
        pass

    def Vp(self, tiempo, tension):
        """  devuelve el valor pico max """
        return np.max(tension)

    def Vrms(self, tiempo, tension):
        """ retorna el valor RMS de la señal """
        return np.sqrt(np.average(tension**2))

    def Vmed(self, tiempo, tension):
        """ retorna el valor medio de modulo de la señal"""
        return np.average(tension)

    def Indice_MOD(self, tiempo, tension):
        """ retorna el indice de modulacion de una señal modulada en AM"""
        pass

    def Delta_f(self, tiempo, tension, fc):
        """ devuelve el valor de desviacion en frecuencia dada una frecuencia de portadora fc"""
        pass

    def THD(self,time,voltage):
        """Calculo de la distorsion armonica."""
        # Calculo la fft, quedandome solo con el espectro positivo y saco la continua.
        
        yf = np.fft.fft(voltage)
        yf = yf[1:floor(len(yf)/2)]
        yf = np.abs(yf) # Obtengo el modulo

        # Calculo el indice donde esta la frecuencia fundamental
        f0_index = np.argmax(yf)
        
        # Armo vector con los indices de los armonicos
        harmonics__index = np.arange(f0_index,len(yf)-1,f0_index)
        
        # Creo vector con valores de los harmonicos
        harmonics_values = yf[harmonics__index]

        # Calculo thd
        thd = np.sqrt(np.sum(harmonics_values[1:]**2))/harmonics_values[0]
 
        return thd
    
    def calculo_Capacitor(self, valor_r, tiempo, tension_r, tension_gen, modo = "TIEMPO"):
        '''
        Calculo de capacitor por distintos metodos
        
        Parameters
        ----------
        valor_r : INT
            valor de resistencia en ohms.
        tiempo : vector
            Vector de tiempos de muestras.
        tension_r : Vector
            Vector de valores de tension [V] en la resistencia.
        tension_gen : Vector
            Vector de valores de tension [V] en el generador
        modo : {TIEMPO; LISSAJ; POT; FFT} string en mayusculas
           Se elije de que modo se calculara el capacitor (Default = TIEMPO)

        Returns
        -------
        valor del capacitor calculado por el metodo indicado

        '''
        # Determino que metodo de medicion utilizo
        if modo == "FFT": 
            
            # obtengo la frecuencia de muestreo
            fs=1/(tiempo[1]-tiempo[0]) 
            
            # Rango del espectro de fourier
            fcia=np.linspace(0,fs/2,len(tiempo)//2)
            
            #Aplico una ventana tipo flat top para mayor presicion
            window = dsp.flattop(len(fcia)) # Ventana flat top
            tension_r *= window
            tension_gen *= window
            
            # obtengo modulo de transformada de fourier de las señales en cuestion
            fft_r   = np.abs(np.fft.fft(tension_r)) 
            fft_r   = fft_r[0:len(fft_r)//2] # elimino espectro repetido
            fft_r   /= len(fcia) # desnormalizo la fft
            
            fft_gen = np.abs(np.fft.fft(tension_gen))
            fft_gen = fft_gen[0:len(fft_gen)//2] # elimino espectro repetido
            fft_gen /= len(fcia) # desnormalizo la fft
            
            #obtengo valores maximos de amplitud y sus fases
            valor_max_fft_r     = np.max(fft_r) # valor pico de la señal
            posicion_max_fft_r  = np.where(fft_r == np.max(fft_r))
            angulo_fft_r =np.angle(np.fft.fft(tension_r)[posicion_max_fft_r])
            
            valor_max_fft_gen     = np.max(fft_gen) # valor pico de la señal
            posicion_max_fft_gen  = np.where(fft_gen == np.max(fft_gen))
            angulo_fft_gen =np.angle(np.fft.fft(tension_gen)[posicion_max_fft_gen])
            
            #calculo del capacitor
            
            frecuencia_Señal = fcia[posicion_max_fft_r] # Frecuencia del pico (de la senoidal)
            itot= valor_max_fft_r/valor_r
            zt = (valor_max_fft_gen/itot) #No se multiplica las ventanas por el coef de correccion ya que se cancelan

            angle_total = angulo_fft_r - angulo_fft_gen 
            xc = zt * np.sin(angle_total) 
            valor_cap = 1/(2 * np.pi * frecuencia_Señal *xc) #valor calculado

        elif modo == "POT":     # Falta implementar
            valor_cap = 0
            pass
        elif modo == "LISSAJ":
            valor_cap = calculo_rc_lissaj(tension_r,tension_gen,tiempo,valor_r)
            pass
        elif modo == "TIEMPO":  # Falta implementar
            valor_cap = 0
            pass

        return valor_cap

def get_zeros_signal(signal_in, delta_min = 5):
      zero_cross = list()
      cross_slope = list()
      i_last = 0
      for i in range(1, signal_in.shape[0]):

        if (signal_in[i-1] <= 0) and (signal_in[i] > 0):
          if (i_last == 0) or ((i-i_last)>delta_min):
            zero_cross.append(i)
            cross_slope.append(1)
            i_last = i
        if (signal_in[i-1] >= 0) and (signal_in[i] < 0):
          if (i_last == 0) or ((i-i_last)>delta_min):
            zero_cross.append(i)
            cross_slope.append(-1)
            i_last = i

      return zero_cross, cross_slope


def calculo_rc_lissaj(vr,vg,t,r):
    ###########################################################################
    # Paso 1: Le saco el valor medio y filtro la señal
    ###########################################################################
    vg = vg - np.average(vg)
    vr = vr - np.average(vr)
    vg = dsp.savgol_filter(vg, 59, 4) # 59 = ventana ; 4 = grado de polinomio
    vr = dsp.savgol_filter(vr, 59, 4)

    ###############################################################################
    # Paso 2: Separar en ciclos a la señal del generador
    ###############################################################################
    # Busco los cruces por cero, de la señal de tensón para dividir los ciclos
    cruces_list, slopes_list = get_zeros_signal(vg)
    limits_ciclos = list()
    zeros_mid_list = []
    trigger_slope = slopes_list[0] # Primera pendiente
    current_limits = np.zeros(2, dtype=np.int32)
    current_limits[0] = cruces_list[0] # Primer cruce
    for cruce, slope in zip(cruces_list, slopes_list):
      # Ignoro el primero
      if current_limits[0] == cruce:
        continue
      # Si la pendiente es igual a la inicial quiere decir que arrancÃ³ de nuevo
      if (trigger_slope == slope):
        # Agrego el final del ciclo
        current_limits[1] = cruce
        limits_ciclos.append(current_limits)
        # Donde termina una empiza el siguiente...
        current_limits = np.zeros(2, dtype=np.int32)
        current_limits[0] = cruce
      # Quiere decir q es un cruce por cero a mitad de ciclo
      else:
        #print(cruce)
        zeros_mid_list.append(cruce)
    #Elimino ultimo cruce por cero
    zeros_mid_list.pop()
    zeros_mid_list = np.array(zeros_mid_list)
    
    ###############################################################################
    # Paso 3: Calculo la W del generador
    ###############################################################################
    ciclos = len(limits_ciclos)
    posMin = limits_ciclos[0][0]
    posMax = limits_ciclos[-1][1]
    tiempoTotal = t[posMax]-t[posMin]
    T0 = tiempoTotal / ciclos
    f = 1/T0
    w = 2*np.pi*f
    #print(f"Ciclos completos: {ciclos}")
    #print(f"Tiempo entre ciclos: {tiempoTotal}")
    #print(f"T: {T0}")
    #print(f"F: {f}")
    #print(f"W: {w}")
    
    ###############################################################################
    # Paso 4: Calcular A y B para cada figura de lisajouse de cada ciclo
    ###############################################################################
    ciclos_vg_list = list()
    ciclos_vr_list = list()
    zeros_vg = list()
    ciclo = 0
    for lim_ciclo in limits_ciclos:
      ciclos_vg_list.append(np.copy(vg[lim_ciclo[0]:lim_ciclo[1]]))
      ciclos_vr_list.append(np.copy(vr[lim_ciclo[0]:lim_ciclo[1]]))
      zeros_vg.append(zeros_mid_list[ciclo] - lim_ciclo[0])
      ciclo+=1
    
    A = 0
    B = 0
    for vg_this, vr_this, zero in zip(ciclos_vg_list, ciclos_vr_list, zeros_vg):
      A += np.max(vr_this) - np.min(vr_this)
      B += np.abs( vr_this[0] - vr_this[zero] )
    A /= ciclos
    B /= ciclos
    #print("A=",A)
    #print("B=",B)
    
    ###############################################################################
    # Paso 5: Calculo la fase
    ###############################################################################
    fase = np.abs( np.arcsin(B/A) )
    #print("Fase = ", fase*180/np.pi)
    
    ###############################################################################
    # Paso 6: Calculo el capacitor
    ###############################################################################
    #f = atg(XC/R))
    #tg(f) = XC/R
    #tg(f) R = 1/WC
    #C = 1/( tg(f)RW )
    C = 1/(np.tan(fase) * r * w)
    return C

