#!/bin/bash

# Continuación: Creando archivos de ejemplos, tests y docs

cd Physics_Numerical_Methods

# ==================== examples/pendulum_simulation.cpp ====================
cat > examples/pendulum_simulation.cpp << 'EOF'
#include "physics_simulations.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== Simulación de Péndulo Simple ===\n\n";
    
    // Crear péndulo: g=9.81 m/s², L=1.0 m
    physics::SimplePendulum pendulum(9.81, 1.0);
    
    // Condiciones iniciales
    double theta0 = M_PI / 4.0;  // 45 grados
    double omega0 = 0.0;          // Velocidad inicial = 0
    
    std::cout << "Parámetros:\n";
    std::cout << "  g = 9.81 m/s²\n";
    std::cout << "  L = 1.0 m\n";
    std::cout << "  θ₀ = " << theta0 * 180 / M_PI << "°\n";
    std::cout << "  ω₀ = " << omega0 << " rad/s\n\n";
    
    // Simular 10 segundos
    pendulum.simulate(theta0, omega0, 10.0, 1000, "pendulum_data.csv");
    
    std::cout << "\n✓ Datos guardados en pendulum_data.csv\n";
    std::cout << "\nPara visualizar, ejecuta:\n";
    std::cout << "  python3 visualize_pendulum.py\n";
    
    return 0;
}
EOF

# ==================== examples/harmonic_oscillator.cpp ====================
cat > examples/harmonic_oscillator.cpp << 'EOF'
#include "physics_simulations.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== Oscilador Armónico Amortiguado ===\n\n";
    
    // Parámetros
    double omega = 2.0 * M_PI;  // ω = 2π rad/s (periodo de 1 segundo)
    double gamma = 0.1;          // γ = 0.1 s⁻¹ (amortiguamiento débil)
    
    physics::DampedOscillator oscillator(omega, gamma);
    
    // Condiciones iniciales
    double x0 = 1.0;   // Posición inicial: 1 metro
    double v0 = 0.0;   // Velocidad inicial: 0
    
    std::cout << "Parámetros:\n";
    std::cout << "  ω = " << omega << " rad/s\n";
    std::cout << "  γ = " << gamma << " s⁻¹\n";
    std::cout << "  x₀ = " << x0 << " m\n";
    std::cout << "  v₀ = " << v0 << " m/s\n\n";
    
    // Tipo de amortiguamiento
    double discriminant = gamma * gamma - omega * omega;
    if (discriminant < 0) {
        std::cout << "Tipo: Subamortiguado (oscila)\n\n";
    } else if (discriminant > 0) {
        std::cout << "Tipo: Sobreamortiguado (no oscila)\n\n";
    } else {
        std::cout << "Tipo: Críticamente amortiguado\n\n";
    }
    
    // Simular 10 segundos
    oscillator.simulate(x0, v0, 10.0, 1000, "oscillator_data.csv");
    
    std::cout << "✓ Datos guardados en oscillator_data.csv\n";
    std::cout << "\nPara visualizar, ejecuta:\n";
    std::cout << "  python3 visualize_oscillator.py\n";
    
    return 0;
}
EOF

# ==================== examples/planetary_motion.cpp ====================
cat > examples/planetary_motion.cpp << 'EOF'
#include "ode_solvers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>

class PlanetaryMotion {
private:
    double G;  // Constante gravitacional
    double M;  // Masa del sol
    
public:
    PlanetaryMotion(double mass_sun) : G(6.67430e-11), M(mass_sun) {}
    
    // Estado: [x, y, vx, vy]
    numerical::StateVector derivatives(double t, const numerical::StateVector& state) {
        double x = state[0], y = state[1];
        double vx = state[2], vy = state[3];
        
        double r = std::sqrt(x*x + y*y);
        double r3 = r * r * r;
        
        double ax = -G * M * x / r3;
        double ay = -G * M * y / r3;
        
        return {vx, vy, ax, ay};
    }
    
    void simulate(const numerical::StateVector& initial, 
                  double t_max, int steps, const std::string& filename) {
        
        auto f = [this](double t, const numerical::StateVector& s) {
            return this->derivatives(t, s);
        };
        
        auto result = numerical::runge_kutta_4_system(f, initial, 0.0, t_max, steps);
        
        std::ofstream file(filename);
        file << "time,x,y,vx,vy,r,energy\n";
        
        double dt = t_max / steps;
        for (int i = 0; i <= steps; i++) {
            double t = i * dt;
            double x = result[i][0];
            double y = result[i][1];
            double vx = result[i][2];
            double vy = result[i][3];
            
            double r = std::sqrt(x*x + y*y);
            double v_sq = vx*vx + vy*vy;
            double energy = 0.5 * v_sq - G * M / r;
            
            file << t << "," << x << "," << y << "," 
                 << vx << "," << vy << "," << r << "," << energy << "\n";
        }
        
        file.close();
    }
};

int main() {
    std::cout << "=== Movimiento Planetario (Tierra alrededor del Sol) ===\n\n";
    
    // Masa del Sol
    double M_sun = 1.989e30;  // kg
    
    PlanetaryMotion simulation(M_sun);
    
    // Condiciones iniciales de la Tierra
    double AU = 1.496e11;           // 1 Unidad Astronómica en metros
    double v_earth = 29783;         // Velocidad orbital en m/s
    
    numerical::StateVector initial = {AU, 0.0, 0.0, v_earth};
    
    std::cout << "Parámetros:\n";
    std::cout << "  Masa del Sol: " << M_sun << " kg\n";
    std::cout << "  Distancia inicial: " << AU/1e9 << " Gm\n";
    std::cout << "  Velocidad orbital: " << v_earth << " m/s\n\n";
    
    // Simular 1 año (31,536,000 segundos)
    double year = 365.25 * 24 * 3600;
    
    std::cout << "Simulando 1 año terrestre...\n";
    simulation.simulate(initial, year, 10000, "orbit_data.csv");
    
    std::cout << "\n✓ Datos guardados en orbit_data.csv\n";
    std::cout << "\nPara visualizar la órbita:\n";
    std::cout << "  python3 visualize_orbit.py\n";
    
    return 0;
}
EOF

# ==================== examples/heat_equation.cpp ====================
cat > examples/heat_equation.cpp << 'EOF'
#include "physics_simulations.h"
#include <iostream>

int main() {
    std::cout << "=== Ecuación del Calor 1D ===\n\n";
    
    // Parámetros
    double alpha = 0.01;  // Difusividad térmica (m²/s)
    double L = 1.0;       // Longitud del dominio (m)
    int nx = 100;         // Número de puntos espaciales
    
    std::cout << "Parámetros:\n";
    std::cout << "  α (difusividad) = " << alpha << " m²/s\n";
    std::cout << "  L (longitud) = " << L << " m\n";
    std::cout << "  Puntos espaciales = " << nx << "\n";
    std::cout << "  Condición inicial: Pulso gaussiano\n";
    std::cout << "  Condiciones de frontera: T = 0 en los bordes\n\n";
    
    physics::HeatEquation1D heat(alpha, L, nx);
    
    // Resolver para 1 segundo
    double t_max = 1.0;
    int nt = 1000;
    
    std::cout << "Resolviendo ecuación del calor...\n";
    std::cout << "Tiempo de simulación: " << t_max << " s\n";
    std::cout << "Pasos temporales: " << nt << "\n\n";
    
    heat.solve(t_max, nt, "heat_data.csv");
    
    std::cout << "\n✓ Datos guardados en heat_data.csv\n";
    std::cout << "\nPara visualizar la evolución térmica:\n";
    std::cout << "  python3 visualize_heat.py\n";
    
    return 0;
}
EOF

# ==================== tests/unit_tests.cpp ====================
cat > tests/unit_tests.cpp << 'EOF'
#include "ode_solvers.h"
#include "integration.h"
#include "numerical_analysis.h"
#include <iostream>
#include <cmath>
#include <cassert>

void test_euler_method() {
    std::cout << "Testing Euler Method... ";
    
    // dy/dt = -2y, y(0) = 1
    // Solución exacta: y(t) = e^(-2t)
    auto f = [](double t, double y) { return -2.0 * y; };
    
    auto result = numerical::euler_method(f, 1.0, 0.0, 1.0, 1000);
    double exact = std::exp(-2.0);
    double error = std::abs(result.back() - exact);
    
    assert(error < 0.01);  // Error tolerable para Euler
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_runge_kutta_4() {
    std::cout << "Testing Runge-Kutta 4... ";
    
    auto f = [](double t, double y) { return -2.0 * y; };
    
    auto result = numerical::runge_kutta_4(f, 1.0, 0.0, 1.0, 100);
    double exact = std::exp(-2.0);
    double error = std::abs(result.back() - exact);
    
    assert(error < 1e-6);  // RK4 es mucho más preciso
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_trapezoidal_rule() {
    std::cout << "Testing Trapezoidal Rule... ";
    
    // Integrar x² de 0 a 1 (resultado exacto = 1/3)
    auto f = [](double x) { return x * x; };
    
    double result = numerical::trapezoidal_rule(f, 0.0, 1.0, 1000);
    double exact = 1.0 / 3.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_simpson_rule() {
    std::cout << "Testing Simpson Rule... ";
    
    auto f = [](double x) { return x * x; };
    
    double result = numerical::simpson_rule(f, 0.0, 1.0, 1000);
    double exact = 1.0 / 3.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-10);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_newton_raphson() {
    std::cout << "Testing Newton-Raphson... ";
    
    // Encontrar √2 (raíz de x² - 2 = 0)
    auto f = [](double x) { return x*x - 2.0; };
    auto df = [](double x) { return 2.0 * x; };
    
    double result = numerical::newton_raphson(f, df, 1.0, 1e-10);
    double exact = std::sqrt(2.0);
    double error = std::abs(result - exact);
    
    assert(error < 1e-10);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_bisection_method() {
    std::cout << "Testing Bisection Method... ";
    
    // Encontrar raíz de x³ - x - 2 = 0 en [1, 2]
    auto f = [](double x) { return x*x*x - x - 2.0; };
    
    double result = numerical::bisection_method(f, 1.0, 2.0, 1e-6);
    double exact = 1.5213797068;  // Raíz conocida
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_gauss_quadrature() {
    std::cout << "Testing Gauss Quadrature... ";
    
    // Integrar sin(x) de 0 a π (resultado = 2)
    auto f = [](double x) { return std::sin(x); };
    
    double result = numerical::gauss_quadrature(f, 0.0, M_PI);
    double exact = 2.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-4);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_numerical_derivative() {
    std::cout << "Testing Numerical Derivative... ";
    
    // Derivada de x² en x=3 (resultado = 6)
    auto f = [](double x) { return x * x; };
    
    double result = numerical::numerical_derivative(f, 3.0);
    double exact = 6.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════╗\n";
    std::cout << "║  Physics Numerical Methods - Unit Tests ║\n";
    std::cout << "╚══════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    try {
        test_euler_method();
        test_runge_kutta_4();
        test_trapezoidal_rule();
        test_simpson_rule();
        test_newton_raphson();
        test_bisection_method();
        test_gauss_quadrature();
        test_numerical_derivative();
        
        std::cout << "\n";
        std::cout << "════════════════════════════════════════════\n";
        std::cout << "  ✓ ALL TESTS PASSED (" << 8 << "/8)\n";
        std::cout << "════════════════════════════════════════════\n";
        std::cout << "\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ TEST FAILED: " << e.what() << "\n";
        return 1;
    }
}
EOF

# ==================== docs/theory.md ====================
cat > docs/theory.md << 'EOF'
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

EOF

echo "✅ Examples, tests y docs creados"

# ==================== Scripts de visualización Python ====================

cat > visualize_pendulum.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('pendulum_data.csv')

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Posición angular vs tiempo
axes[0, 0].plot(data['time'], data['theta'] * 180/np.pi)
axes[0, 0].set_xlabel('Tiempo (s)')
axes[0, 0].set_ylabel('Ángulo (°)')
axes[0, 0].set_title('Posición Angular vs Tiempo')
axes[0, 0].grid(True, alpha=0.3)

# Velocidad angular vs tiempo
axes[0, 1].plot(data['time'], data['omega'], color='orange')
axes[0, 1].set_xlabel('Tiempo (s)')
axes[0, 1].set_ylabel('Velocidad Angular (rad/s)')
axes[0, 1].set_title('Velocidad Angular vs Tiempo')
axes[0, 1].grid(True, alpha=0.3)

# Espacio de fases
axes[1, 0].plot(data['theta'] * 180/np.pi, data['omega'])
axes[1, 0].set_xlabel('Ángulo (°)')
axes[1, 0].set_ylabel('Velocidad Angular (rad/s)')
axes[1, 0].set_title('Espacio de Fases')
axes[1, 0].grid(True, alpha=0.3)

# Energía vs tiempo
axes[1, 1].plot(data['time'], data['energy'], color='green')
axes[1, 1].set_xlabel('Tiempo (s)')
axes[1, 1].set_ylabel('Energía (J)')
axes[1, 1].set_title('Energía Total vs Tiempo')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pendulum_analysis.png', dpi=300, bbox_inches='tight')
print("✓ Gráfico guardado: pendulum_analysis.png")
plt.show()
EOF

chmod +x visualize_pendulum.py

echo "✅ Scripts de visualización creados"

# ==================== Makefile alternativo ====================
cat > Makefile << 'EOF'
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Iinclude
LDFLAGS =

SRC_DIR = src
INC_DIR = include
EXAMPLES_DIR = examples
TEST_DIR = tests
BIN_DIR = bin

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

EXAMPLES = pendulum_simulation harmonic_oscillator planetary_motion heat_equation
EXAMPLE_BINS = $(addprefix $(BIN_DIR)/, $(EXAMPLES))

.PHONY: all clean examples tests

all: $(BIN_DIR) $(OBJECTS) examples tests

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

examples: $(EXAMPLE_BINS)

$(BIN_DIR)/%: $(EXAMPLES_DIR)/%.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LDFLAGS)

tests: $(BIN_DIR)/unit_tests
	./$(BIN_DIR)/unit_tests

$(BIN_DIR)/unit_tests: $(TEST_DIR)/unit_tests.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)/*.o $(BIN_DIR)/$(EXAMPLES) $(BIN_DIR)/unit_tests
	rm -f *.csv *.png

run-all: examples
	@echo "Ejecutando todas las simulaciones..."
	@./$(BIN_DIR)/pendulum_simulation
	@./$(BIN_DIR)/harmonic_oscillator
	@./$(BIN_DIR)/planetary_motion
	@./$(BIN_DIR)/heat_equation

EOF

echo ""
echo "════════════════════════════════════════════════════════"
echo "  ✅ ESTRUCTURA COMPLETA CREADA"
echo "════════════════════════════════════════════════════════"
echo ""
echo "📁 Directorios creados:"
echo "  ✓ include/     - Headers (.h)"
echo "  ✓ src/         - Implementaciones (.cpp)"
echo "  ✓ examples/    - Programas de ejemplo"
echo "  ✓ tests/       - Tests unitarios"
echo "  ✓ docs/        - Documentación teórica"
echo ""
echo "🚀 Para compilar y ejecutar:"
echo ""
echo "  cd Physics_Numerical_Methods"
echo "  mkdir build && cd build"
echo "  cmake .."
echo "  make"
echo "  ./bin/pendulum_simulation"
echo ""
echo "📝 O usando Makefile:"
echo ""
echo "  make all"
echo "  make run-all"
echo ""
echo "🧪 Para ejecutar tests:"
echo ""
echo "  make tests"
echo ""
echo "════════════════════════════════════════════════════════"
