# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 09:57:19 2018

@author: Ramiro
"""

# Traemos la libreria VISA
import pyvisa as visa
# Traemos matplotlib para poder graficar
import matplotlib.pyplot as plt
# Agreamos el path de las librerias
import sys
sys.path.insert(0,'InstVirtualLib')
import platform
# Traemos todos los osciloscopios
from InstVirtualLib.osciloscopios import GW_Instek
from InstVirtualLib.osciloscopios import rigol
from InstVirtualLib.osciloscopios import Tektronix_DSO_DPO_MSO_TDS
# Traemos el operador
import operador
import numpy as np
import csv ##usado para guardar csv de otra forma
# Definimos una funcion para poder ejecutar un mensaje de error
def excepthook(type, value, traceback):
    print(value)

sys.excepthook = excepthook

plt.close('all')

# Seteamos el tipo de osciloscio a utilizar
OSCILOSCOPIOS = 0	# 0: GW_Instek
        			# 1: Rigol
			        # 2: Tektronix_DSO_DPO_MSO_TDS

USE_DEVICE = 1



# Abrimos el instrumento
platforma = platform.platform()
print(platforma)
rm=visa.ResourceManager()
# El handle puede controlar osciloscopios, analizadores, etc.
# Es generico
instrument_handler=rm.open_resource(rm.list_resources()[USE_DEVICE])

if OSCILOSCOPIOS == 0:
	MiOsciloscopio = GW_Instek(instrument_handler)
elif OSCILOSCOPIOS == 1:
	MiOsciloscopio = rigol(instrument_handler)
elif OSCILOSCOPIOS == 2:
	MiOsciloscopio = Tektronix_DSO_DPO_MSO_TDS(instrument_handler)
else:
	raise ValueError('Tipo de osciloscopio fuera de lista.')


# Informamos el modelo del osciloscopio conectado
print("Esta conectado un %s"%MiOsciloscopio.INSTR_ID)


# Pedimos el trazo de cada canal, la salida es en ([seg.],[volt])
# BUG!!! si no se agrega el VERBOSE=False no anda el GW Instek
print('-------------')
tiempo1,tension1=MiOsciloscopio.get_trace("1",VERBOSE=False)
print('-------------')
tiempo2,tension2=MiOsciloscopio.get_trace("2",VERBOSE=False)
print('-------------')

# Ploteamos los canales
fig_RC= plt.figure(1)
ax_RC, ax_R= fig_RC.subplots(2)

fig_RC.sca(ax_RC)
ax_RC.plot(tiempo1,tension1, label='Tension sobre RC')
plt.legend()
fig_RC.sca(ax_R)
ax_R.plot(tiempo2,tension2, color='red', label='Tension sobre R')
plt.legend()
plt.show()
print('LINEA81')

# TODO: hacer que no imprima de vuelta los datos del canal
# Generamos un operador y pedimos el valor RMS actual
operador_1 = operador.Operador_osciloscopio(MiOsciloscopio,"Workbench_I")

val_RMS = operador_1.medir_Vrms(canal = 1, VERBOSE = False)

print('Vrms = %0.5f'%val_RMS)

np.savetxt('tension1.csv', tension1, delimiter=',')
np.savetxt('tension2.csv', tension2, delimiter=',')
np.savetxt('tiempo1.csv', tiempo1, delimiter=',')
np.savetxt('tiempo2.csv', tiempo2, delimiter=',')

# f=open("datos.csv",'w')
# with f:
#     writer=csv.writer(f)
#     for row in tension1:
#         writer.writerow(row)



MiOsciloscopio.close()



### Datos circuito
R= 1200
Rg= 50
C= 220 * 10**-9 #220 nF
f_gen= 999.917 # Medida por el osciloscopio, el generador sacaba 1kHz

### Lecturas osciloscopio
fdv= 500*10**-3
k_punta= 1

n_div_tension= 2
n_div_corriente= 1.4

vg_osc= n_div_tension * fdv * k_punta
ig_osc= n_div_corriente * fdv * k_punta / R

### |Z|= |V| / |I|

z= vg_osc/ig_osc

delta_N_vg= (0.1/n_div_tension) * 100
delta_N_ig= (0.1/n_div_corriente) * 100
delta_fbt= 0.03 * 100
delta_punta= 0.02* 100
delta_R= 0.1 * 100

uc_vg_relativa= np.sqrt(delta_N_vg**2 + delta_fbt**2 + delta_punta**2)
uc_ig_relativa= np.sqrt(delta_N_ig**2 + delta_fbt**2 + delta_punta**2 + delta_R**2)

uc_z_relativa= np.sqrt(uc_vg_relativa**2 + uc_ig_relativa**2)

print(f"Z= {round(z,2)} Ω +- {round(uc_z_relativa,2)} % @68% de confianza")

### Calculo de Angulos

# Formula:  delta_t= n_div_horiz * fdh
fdh= 250 * 10**-6 #250 us/div
delta_t= 1 * fdh
periodo_T= 8.3 * fdh

alfa= delta_t * 360 / periodo_T ## Medido en grados

print(f"Fase de la impedancia: α= {round(alfa,2)}")

### Reactancia capacitiva: Xc= |Z| * sen(alfa) = 1/ (2* pi *f * C)

Xc= z * np.sin(alfa*np.pi/180)

c_calculado= 1/(2* np.pi * f_gen * Xc)

print(f"\nValor del capacitor: C= {round(c_calculado*10**9, 2)} nF")




