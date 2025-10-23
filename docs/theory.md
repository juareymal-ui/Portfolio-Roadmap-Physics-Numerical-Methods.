# 📐 Fundamentos Teóricos

## Tabla de Contenidos

- [Ecuaciones Diferenciales Ordinarias](#ecuaciones-diferenciales-ordinarias)
- [Métodos de Integración Numérica](#métodos-de-integración-numérica)
- [Análisis de Errores](#análisis-de-errores)
- [Estabilidad Numérica](#estabilidad-numérica)

---

## Ecuaciones Diferenciales Ordinarias

### Método de Euler

El método más simple para resolver EDOs de la forma:

```
dy/dt = f(t, y)
```

**Fórmula de iteración:**

```
y(n+1) = y(n) + h * f(t(n), y(n))
```

donde `h` es el tamaño del paso.

**Propiedades:**
- **Orden:** O(h²) error local, O(h) error global
- **Ventajas:** Simple, fácil de implementar
- **Desventajas:** Baja precisión, puede ser inestable

**Ejemplo físico:** Caída libre

```
dv/dt = -g
v(n+1) = v(n) - g*h
```

---

### Método de Runge-Kutta de 4to Orden (RK4)

Método de mayor precisión que usa 4 evaluaciones por paso:

**Fórmula:**

```
k₁ = h * f(t, y)
k₂ = h * f(t + h/2, y + k₁/2)
k₃ = h * f(t + h/2, y + k₂/2)
k₄ = h * f(t + h, y + k₃)

y(n+1) = y(n) + (k₁ + 2k₂ + 2k₃ + k₄) / 6
```

**Propiedades:**
- **Orden:** O(h⁵) error local, O(h⁴) error global
- **Ventajas:** Muy preciso, ampliamente usado
- **Desventajas:** 4 evaluaciones por paso

---

## Métodos de Integración Numérica

### Regla del Trapecio

Aproxima el área bajo la curva usando trapecios:

```
∫ₐᵇ f(x)dx ≈ h/2 * [f(a) + 2∑f(xᵢ) + f(b)]
```

**Error:** O(h²)

### Regla de Simpson

Usa interpolación parabólica:

```
∫ₐᵇ f(x)dx ≈ h/3 * [f(a) + 4∑f(x₂ᵢ₋₁) + 2∑f(x₂ᵢ) + f(b)]
```

**Error:** O(h⁴)

### Cuadratura de Gauss

Elige puntos óptimos para máxima precisión con mínimas evaluaciones.

**Para n=5 puntos:** Exacto para polinomios de grado ≤ 9

---

## Análisis de Errores

### Error de Truncamiento

Error introducido al truncar la serie de Taylor:

```
Error local = O(hᵖ⁺¹)
Error global = O(hᵖ)
```

donde `p` es el orden del método.

### Error de Redondeo

Error debido a precisión finita de punto flotante.

**Regla práctica:**
- `double` (64 bits): ~16 dígitos decimales
- Evitar sumas de números muy diferentes en magnitud

---

## Estabilidad Numérica

### Criterio de Estabilidad

Para ecuación del calor:

```
∂u/∂t = α ∂²u/∂x²
```

**Criterio CFL:**

```
α * Δt / Δx² ≤ 0.5
```

Si se viola, la solución puede volverse inestable y diverger.

### Métodos Implícitos vs Explícitos

**Explícitos:** Más simples, pero restricciones de estabilidad
**Implícitos:** Incondicionalmente estables, pero requieren resolver sistemas

---

## Problemas Físicos Clásicos

### Péndulo Simple

**Ecuación:**

```
d²θ/dt² + (g/L)sin(θ) = 0
```

**Aproximación lineal (ángulos pequeños):**

```
d²θ/dt² + ω²θ = 0, donde ω² = g/L
```

**Período:** T = 2π√(L/g)

### Oscilador Armónico Amortiguado

**Ecuación:**

```
d²x/dt² + 2γ dx/dt + ω²x = 0
```

**Tipos de amortiguamiento:**
- **Subamortiguado:** γ < ω (oscila)
- **Crítico:** γ = ω
- **Sobreamortiguado:** γ > ω (no oscila)

### Movimiento Planetario

**Segunda Ley de Newton + Gravitación:**

```
F = -GMm/r² r̂
ma = F
```

**Conservación:**
- Energía total: E = K + U = constante
- Momento angular: L = r × p = constante

---

## Referencias

1. **Numerical Recipes** - Press et al.
2. **Computational Physics** - Mark Newman
3. **An Introduction to Computer Simulation Methods** - Gould & Tobochnik

---

## Notación

- `h, Δt`: Paso temporal
- `Δx`: Paso espacial
- `O(hⁿ)`: Orden del error
- `f(t,y)`: Función que define la EDO
- `α`: Difusividad térmica
- `γ`: Coeficiente de amortiguamiento
- `ω`: Frecuencia angular

