# MiModuloBiofisica

## Descripción del proyecto

MiModuloBiofisica es una librería desarrollada en Python que contiene funciones útiles para cálculos matemáticos, biofísicos y estadísticos.
El módulo permite resolver problemas relacionados con energía cinética, ajuste de curvas, ecuaciones diferenciales, propagación de errores, distribuciones normales y cinética enzimática.

---

# Funciones implementadas

## Función 1: energia_cinetica()

Recibe masa (kg) y velocidad (m/s).
Calcula la energía cinética del cuerpo.
Devuelve la energía cinética en Joules.

Resultado esperado:

| Función            | Resultado                             |
| ------------------ | ------------------------------------- |
| energia_cinetica() | Energía cinética del objeto en Joules |

---

## Función 2: ajustar_curva()

Recibe datos experimentales (x, y) y un modelo matemático.
Ajusta la curva por regresión no lineal.

Devuelve:

* params: valores óptimos de los parámetros
* errores: incertidumbre de cada parámetro
* y_pred: valores predichos por el modelo

Resultado esperado:

| Función         | Resultado                              |
| --------------- | -------------------------------------- |
| ajustar_curva() | Parámetros óptimos y valores ajustados |

---

## Función 3: resolver_edo()

Recibe una ecuación diferencial dy/dt y una condición inicial y(0).
Resuelve la EDO numéricamente en un rango de tiempo.

Devuelve:

* y(t): valores de la variable en cada instante
* t: arreglo de puntos de tiempo

Resultado esperado:

| Función        | Resultado                                    |
| -------------- | -------------------------------------------- |
| resolver_edo() | Solución numérica de la ecuación diferencial |

---

## Función 4: error_propagacion()

Recibe una función f, sus variables y sus errores absolutos.
Calcula el error usando derivadas parciales.

Devuelve:

* Δf: error absoluto de la función
* Δf/f: error relativo

Resultado esperado:

| Función             | Resultado                 |
| ------------------- | ------------------------- |
| error_propagacion() | Error absoluto y relativo |

---

## Función 5: probabilidad_normal()

Recibe media μ, desviación estándar σ e intervalo [a, b].
Calcula el área bajo la curva normal en ese intervalo.

Devuelve:

* P(a ≤ x ≤ b): probabilidad en el intervalo dado

Resultado esperado:

| Función               | Resultado                    |
| --------------------- | ---------------------------- |
| probabilidad_normal() | Probabilidad en el intervalo |

---

## Función 6: michaelis_menten()

Recibe concentración de sustrato [S], Vmax y Km.
Calcula la velocidad de reacción enzimática.

Devuelve:

* v: velocidad de reacción en unidades de Vmax

Resultado esperado:

| Función            | Resultado                        |
| ------------------ | -------------------------------- |
| michaelis_menten() | Velocidad de reacción enzimática |

---

# Ejemplos de uso

```python
import mimodulodebiofisica as mb

energia = mb.energia_cinetica(10, 5)

print(energia)
```

---

# Integrantes del grupo

* Bruno Alonso Subauste Tokumura
    * Andrés Céspedes Prada
     . Diego Fabricio Juro Matamoros
     .Castro Pinto, Camila Cielo
