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
from matplotlib import pyplot as plt
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
        elif modo == "LISSAJ":  # Falta implementar
            valor_cap = 0
            pass
        
        elif modo == "TIEMPO":  # Falta implementar
            
            Med = Mediciones()

            Vpico_RC = Med.Vp(tiempo, tension_gen)
            Vpico_R = Med.Vp(tiempo, tension_r)

            ### Datos circuito
            f_gen = 602.85      # Frecuencia asociada al RC

            # Modulo de z

            Ipico_RC = Vpico_R/valor_r

            ### |Z| = |V| / |I|

            z = Vpico_RC/Ipico_RC

            ### Calculo de Angulos

            for i in range(0, 4000):        # Recorro todo el csv
                aux1 = tension_gen[i]
                if aux1 == 0:
                    time1 = tiempo[i]
                    
                aux2 = tension_r[i]
                if aux2 == 0:
                    time2 = tiempo[i]

            delta_t = time1 - time2        # Delta t entre ceros de la senoidal

            periodo_T = 1/f_gen            # Periodo

            alfa = delta_t * 360 / periodo_T ## Medido en grados

            ### Reactancia capacitiva: Xc= |Z| * sen(alfa) = 1/ (2* pi *f * C)

            Xc = z * np.sin(alfa*np.pi/180)

            valor_cap = 1/(2* np.pi * f_gen * Xc)
        
            pass
        
        return valor_cap