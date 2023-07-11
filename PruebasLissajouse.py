"""
Pruebas para MÃ©todo Lisajouse
"""

import numpy as np
from scipy import signal as dsp

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

# Cargar los datos desde el archivo CSV
vg = np.genfromtxt('tension1.csv', delimiter=',') 
vr = np.genfromtxt('tension2.csv', delimiter=',')
t = np.genfromtxt('tiempo1.csv', delimiter=',')
c = calculo_rc_lissaj(vr,vg,t,1200)
print(f"C = {c*1000*1000*1000}nF")