"""
prueba.py
=========
Script de pruebas para mimodulobiofisica.py
Cada prueba está basada en los ejercicios del curso de Fundamentos de Biofísica.
"""

import numpy as np
import matplotlib.pyplot as plt
from mimodulobiofisica import (
    energia_cinetica,
    ajustar_curva,
    resolver_edo,
    error_propagacion,
    probabilidad_normal,
    michaelis_menten,
)

print("=" * 55)
print("       PRUEBAS - MÓDULO DE BIOFÍSICA")
print("=" * 55)


# =============================================================
# PRUEBA 1: energia_cinetica()
# =============================================================
print("\n--- Prueba 1: energia_cinetica() ---")
# Un objeto de 10 kg se mueve a 3 m/s
masa = 10       # kg
vel  = 3        # m/s
Ek = energia_cinetica(masa, vel)
print(f"  Masa     : {masa} kg")
print(f"  Velocidad: {vel} m/s")
print(f"  Resultado: Ek = {Ek} J")
assert Ek == 45.0, "ERROR en energia_cinetica"
print("  [OK]")


# =============================================================
# PRUEBA 2: resolver_edo()
# Ejercicio 3 del examen:
#   dy/dt = -k * y^2   con y(0) = y0, k constante
#   Solución analítica: y(t) = y0 / (1 + k * y0 * t)
# =============================================================
print("\n--- Prueba 2: resolver_edo() ---")
k  = 0.001
y0 = 1000.0
t  = np.linspace(0, 10, 200)

f_bacterias = lambda y, t: -k * y**2

t_sol, y_sol = resolver_edo(f_bacterias, y0, t)

# Solución analítica para comparar
y_analitica = y0 / (1 + k * y0 * t_sol)
error_max = np.max(np.abs(y_sol - y_analitica))

print(f"  k              : {k}")
print(f"  y(0)           : {y0} bacterias")
print(f"  y(t=10) num.   : {y_sol[-1]:.4f}")
print(f"  y(t=10) anali. : {y_analitica[-1]:.4f}")
print(f"  Error máximo   : {error_max:.6f}")
assert error_max < 0.01, "ERROR en resolver_edo"
print("  [OK]")

plt.figure()
plt.plot(t_sol, y_sol, label="Solución numérica")
plt.plot(t_sol, y_analitica, "--", label="Solución analítica")
plt.xlabel("Tiempo (s)")
plt.ylabel("Número de bacterias")
plt.title("Prueba 2: resolver_edo() — Crecimiento bacteriano")
plt.legend()
plt.tight_layout()
plt.savefig("prueba2_edo.png")
print("  Gráfica guardada: prueba2_edo.png")


# =============================================================
# PRUEBA 3: error_propagacion()
# Ejercicio 4 del examen:
#   Volumen de esfera V = (4/3) * pi * r^3
#   r = 1 m,  Δr = 0.02 m
#   Error relativo esperado = 3 * (Δr/r) = 0.06 = 6%
# =============================================================
print("\n--- Prueba 3: error_propagacion() ---")
f_esfera = lambda vals: (4/3) * np.pi * vals[0]**3

V, delta_V, delta_V_rel = error_propagacion(
    f=f_esfera,
    variables=["r"],
    valores=[1.0],
    errores=[0.02]
)

print(f"  Radio          : 1.0 m  ±  0.02 m")
print(f"  Volumen V      : {V:.6f} m³")
print(f"  Error absoluto : ΔV = {delta_V:.6f} m³")
print(f"  Error relativo : ΔV/V = {delta_V_rel:.4f}  ({delta_V_rel*100:.2f}%)")
assert abs(delta_V_rel - 0.06) < 1e-4, "ERROR en error_propagacion"
print("  [OK]")


# =============================================================
# PRUEBA 4: probabilidad_normal()
# Ejercicio 5 del examen:
#   Truchas: media = 28 cm, sigma = 2 cm
#   P(24 <= longitud <= 26) = ?
# =============================================================
print("\n--- Prueba 4: probabilidad_normal() ---")
prob, z_a, z_b = probabilidad_normal(
    media=28,
    desviacion=2,
    a=24,
    b=26
)

print(f"  Media          : 28 cm")
print(f"  Desv. estándar : 2 cm")
print(f"  Intervalo      : [24, 26] cm")
print(f"  z_a = {z_a:.2f},  z_b = {z_b:.2f}")
print(f"  P(24 ≤ x ≤ 26) = {prob:.4f}  ({prob*100:.2f}%)")
assert 0.13 < prob < 0.16, "ERROR en probabilidad_normal"
print("  [OK]")

x_plot = np.linspace(20, 36, 300)
from scipy.stats import norm as sp_norm
y_plot = sp_norm.pdf(x_plot, 28, 2)

plt.figure()
plt.plot(x_plot, y_plot, color="steelblue", label="Distribución normal")
x_fill = np.linspace(24, 26, 100)
plt.fill_between(x_fill, sp_norm.pdf(x_fill, 28, 2),
                 alpha=0.4, color="steelblue", label=f"P = {prob:.4f}")
plt.xlabel("Longitud (cm)")
plt.ylabel("Densidad de probabilidad")
plt.title("Prueba 4: probabilidad_normal() — Truchas")
plt.legend()
plt.tight_layout()
plt.savefig("prueba4_normal.png")
print("  Gráfica guardada: prueba4_normal.png")


# =============================================================
# PRUEBA 5: michaelis_menten()
# Datos típicos de cinética enzimática
# =============================================================
print("\n--- Prueba 5: michaelis_menten() ---")
Vmax = 1.0    # mmol/min
Km   = 2.0    # mM
S_vals = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0])

v_vals = michaelis_menten(S_vals, Vmax, Km)

print(f"  Vmax = {Vmax} mmol/min,  Km = {Km} mM")
print(f"  {'[S] (mM)':<12} {'v (mmol/min)':<15}")
for s, v in zip(S_vals, v_vals):
    print(f"  {s:<12.1f} {v:<15.4f}")

# A S = Km la velocidad debe ser exactamente Vmax/2
v_en_Km = michaelis_menten(Km, Vmax, Km)
assert abs(v_en_Km - Vmax / 2) < 1e-10, "ERROR en michaelis_menten"
print(f"  Verificación: v(S=Km) = {v_en_Km:.4f} = Vmax/2 ✓")
print("  [OK]")

plt.figure()
S_cont = np.linspace(0.01, 20, 300)
plt.plot(S_cont, michaelis_menten(S_cont, Vmax, Km), label="Modelo M-M")
plt.scatter(S_vals, v_vals, color="red", zorder=5, label="Puntos de prueba")
plt.axhline(Vmax, linestyle="--", color="gray", label=f"Vmax = {Vmax}")
plt.axvline(Km,   linestyle=":",  color="orange", label=f"Km = {Km} mM")
plt.xlabel("[S] (mM)")
plt.ylabel("v (mmol/min)")
plt.title("Prueba 5: michaelis_menten()")
plt.legend()
plt.tight_layout()
plt.savefig("prueba5_michaelis.png")
print("  Gráfica guardada: prueba5_michaelis.png")


# =============================================================
# PRUEBA 6: ajustar_curva()
# Ajuste de Michaelis-Menten con datos experimentales con ruido
# =============================================================
print("\n--- Prueba 6: ajustar_curva() ---")
np.random.seed(42)
S_exp = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
v_exp = michaelis_menten(S_exp, Vmax=1.0, Km=2.0) + np.random.normal(0, 0.02, len(S_exp))

def modelo_mm(S, Vmax, Km):
    return (Vmax * S) / (Km + S)
params, errores, y_pred = ajustar_curva(S_exp, v_exp, modelo_mm, p0=[0.8, 1.0])

print(f"  Vmax ajustado : {params[0]:.4f} ± {errores[0]:.4f}")
print(f"  Km ajustado   : {params[1]:.4f} ± {errores[1]:.4f}")
print(f"  Valores reales: Vmax=1.0,  Km=2.0")
assert abs(params[0] - 1.0) < 0.1, "ERROR en ajustar_curva (Vmax)"
assert abs(params[1] - 2.0) < 0.3, "ERROR en ajustar_curva (Km)"
print("  [OK]")

plt.figure()
S_cont = np.linspace(0.01, 20, 300)
plt.scatter(S_exp, v_exp, color="red", zorder=5, label="Datos experimentales")
plt.plot(S_cont, modelo_mm(S_cont, *params), label="Curva ajustada")
plt.xlabel("[S] (mM)")
plt.ylabel("v (mmol/min)")
plt.title("Prueba 6: ajustar_curva() — Ajuste Michaelis-Menten")
plt.legend()
plt.tight_layout()
plt.savefig("prueba6_ajuste.png")
print("  Gráfica guardada: prueba6_ajuste.png")


# =============================================================
print("\n" + "=" * 55)
print("   TODAS LAS PRUEBAS PASARON CORRECTAMENTE")
print("=" * 55)
