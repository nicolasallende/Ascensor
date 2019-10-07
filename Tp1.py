"""
TP1 ANALISIS NUMERICO
Problemas de busqueda de raices
Creando un ascensor
"""

#Imports
import math
import matplotlib
from matplotlib import pyplot

matplotlib.use('Agg')

import numpy as np #Manejo de arrays
import matplotlib.pylab as plt #Rutinas gráficas

##Constantes
fuerza = 280 #Fuerza medida en newtons

masa_cabina = 350 #peso cabina

masa_persona = 75 #peso prom personas

numero_personas = 6 #a determinar

masa = masa_cabina + masa_persona*numero_personas

altura = 3 #altura de un piso medida en metros

acel_maxima = fuerza/masa #aceleracion maxima posible

tiempo_final = math.sqrt(18*masa/fuerza)

tiempoEn30Porciento = ( (fuerza / masa) - acel_maxima * 0.3 ) * (masa/fuerza) * (tiempo_final / 2)


print('-------------------------------------')
print('La fuerza que voy a usar es ', fuerza, '[N]')
print('------------------- ------------------')
print('La masa que voy a usar es ', masa, '[Kg]')
print('-------------------------------------')
print('La altura que voy a usar es ', altura, '[m]')
print('-------------------------------------')
print('La aceleracion maxima es ', acel_maxima, '[m/s^2]')
print('-------------------------------------')

###Funciones posicion (f) y velocidad (f derivada ) y aceleracion (f derivada segunda)
def posicion(t):
    x = (fuerza/masa) * (((t**2) /2) - ((t**3) / (tiempo_final * 3)))
    return x

def velocidad(t):
    v = (fuerza/masa) * (t  - ((t**2) / tiempo_final))
    return v

def aceleracion(t):
    a = ( (fuerza/masa) * (1  - ((t*2) / tiempo_final)))
    return a


#grafico de posicion

x = np.linspace(0, tiempo_final )
y = posicion(x)
plt.figure(figsize = (10,7))
plt.plot(x, y)
pyplot.axhline(0, color="black")
pyplot.axvline(0, color="black")
plt.xlabel('tiempo')
plt.ylabel('altura')
plt.title('Funcion de la posicion')
plt.grid(True)
plt.savefig('Grafico Posicion')

#grafico de velocidad

x = np.linspace(0, tiempo_final )
y = velocidad(x)
plt.figure(figsize = (10,7))
plt.plot(x, y)
pyplot.axhline(0, color="black")
pyplot.axvline(0, color="black")
plt.xlabel('tiempo')
plt.ylabel('vel')
plt.title('Funcion de la velocidad')
plt.grid(True)
plt.savefig('Grafico Velocidad')

#grafico de aceleracion

x = np.linspace(0, tiempo_final )
y = aceleracion(x)
plt.figure(figsize = (10,7))
plt.plot(x, y)
pyplot.axhline(0, color="black")
pyplot.axvline(0, color="black")
plt.xlabel('tiempo')
plt.ylabel('aceleracion')
plt.title('Funcion de la aceleracion')
plt.grid(True)
plt.savefig('Grafico Aceleracion')



def f(x):
    fa = posicion(x) - posicion(tiempoEn30Porciento) #este es el X0
    df = velocidad(x)
    return fa, df
    
TOL = 5E-16    # tolerancia
MAXITER_B = 2  # maxima cantidad de iteraciones para biseccion
MAXITER_N = 60  # maxima cantidad de iteraciones para newton

def Biseccion(a,b):
   formato_salida_B = '%03d %+.17E %+.3E'
   fa = f(a)[0]
   fb = f(b)[0]
   print('-------------------------------------')
   print("Método de bisección")
   if abs(fa) < TOL:
       xo = a
       print(formato_salida_B % (0, xo, fa))
   elif abs(fb) < TOL:
        xo = b
        print(formato_salida_B % (0, xo, fb))
   elif fb * fa < 0:
       k = 0
       fin = False
       while k < MAXITER_B and not fin:
           xo = (a + b) / 2
           f_val = f(xo)[0]
           print(formato_salida_B % (k, xo, f_val))
           k += 1
           if abs(f_val) < TOL:
               fin = True
           elif f(a)[0] * f_val < 0:
               b = xo
           else:
               a = xo
   else:
       print("El intervalo no contiene raíces")
       xo = a
   return xo




def Newton():
    xo = Biseccion(0,tiempo_final)
    formato_salida_N = '%03d %+.17E %+.3E %+.3E'
    print('-------------------------------------')
    print("Método de Newton")
    k = 0
    fin = False
    while k < MAXITER_N and not fin:
        fo = f(xo)[0]
        dfo = f(xo)[1]
        xk = -fo / dfo + xo
        e_u = (xo - xk) / xk
        fin = abs(fo) < TOL and abs(e_u) < TOL
        print(formato_salida_N % (k, xo, fo, e_u))
        xo = xk
        k += 1
    return None

Newton()
