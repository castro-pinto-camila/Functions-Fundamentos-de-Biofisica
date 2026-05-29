import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint, quad
from scipy.stats import norm

def energia_cinetica(masa, velocidad):
    energia = 0.5 * masa * (velocidad ** 2)
    return energia 


def ajustar_curva(x_datos, y_datos, modelo, p0):
   
    params, covarianza = curve_fit(modelo, x_datos, y_datos, p0=p0)
    errores = np.sqrt(np.diag(covarianza))
    y_pred = modelo(x_datos, *params)
    return params, errores, y_pred

def resolver_edo(f, y0, t):
    
    y = odeint(f, y0, t).flatten()
    return t, y
 

def error_propagacion(f, variables, valores, errores, h=1e-6):
    
    valores = np.array(valores, dtype=float)
    errores = np.array(errores, dtype=float)
 
    f_valor = f(valores)
 
    # Derivadas parciales numéricas (diferencia central)
    derivadas = np.zeros(len(valores))
    for i in range(len(valores)):
        vals_mas = valores.copy()
        vals_menos = valores.copy()
        vals_mas[i] += h
        vals_menos[i] -= h
        derivadas[i] = (f(vals_mas) - f(vals_menos)) / (2 * h)
 
    # Propagación cuadrática
    error_absoluto = np.sqrt(np.sum((derivadas * errores) ** 2))
    error_relativo = abs(error_absoluto / f_valor) if f_valor != 0 else float('inf')
 
    return f_valor, error_absoluto, error_relativo
 
 
def probabilidad_normal(media, desviacion, a, b):

    z_a = (a - media) / desviacion
    z_b = (b - media) / desviacion
    probabilidad = norm.cdf(z_b) - norm.cdf(z_a)
    return probabilidad, z_a, z_b
 
 
def michaelis_menten(S, Vmax, Km):
    
    S = np.asarray(S, dtype=float)
    v = (Vmax * S) / (Km + S)
    return v