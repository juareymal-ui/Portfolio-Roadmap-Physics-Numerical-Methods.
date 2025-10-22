#  Métodos Numéricos y Computación Científica en C++

<div align="center">

![C++](https://img.shields.io/badge/C++-00599C?style=for-the-badge&logo=cplusplus&logoColor=white)
![Physics](https://img.shields.io/badge/Physics-Computing-purple?style=for-the-badge)
![Numerical](https://img.shields.io/badge/Numerical-Methods-orange?style=for-the-badge)

**Implementaciones computacionales de métodos numéricos para resolver problemas físicos**

[Instalación](#-instalación) • [Métodos](#-métodos-implementados) • [Ejemplos](#-ejemplos) • [Contribuir](#-contribuir)

</div>

---

##  Tabla de Contenidos

- [Introducción](#-introducción)
- [Requisitos](#-requisitos)
- [Instalación](#-instalación)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Métodos Implementados](#-métodos-implementados)
  - [Ecuaciones Diferenciales](#1-ecuaciones-diferenciales-ordinarias-odes)
  - [Métodos de Integración](#2-métodos-de-integración-numérica)
  - [Simulaciones Físicas](#3-simulaciones-físicas)
  - [Análisis Numérico](#4-análisis-numérico)
- [Ejemplos de Uso](#-ejemplos-de-uso)
- [Documentación](#-documentación)
- [Contribuir](#-contribuir)

---

##  Introducción

Este repositorio contiene implementaciones en **C++** de métodos numéricos fundamentales para resolver problemas en física computacional. Los algoritmos están optimizados para rendimiento y precisión, siendo útiles tanto para investigación como para aprendizaje.

### ¿Por qué C++?

-  **Alto rendimiento** para cálculos intensivos
-  **Control preciso** de memoria y recursos
-  **Amplia compatibilidad** con librerías científicas
-  **Estándar** en computación científica de alto rendimiento

---

##  Requisitos

### Software Necesario

```bash
# Compilador C++ (GCC 9+ o Clang 10+)
g++ --version  # Debe ser >= 9.0

# CMake (opcional, para gestión de proyectos)
cmake --version  # Recomendado >= 3.15
```

### Librerías Requeridas

- **Eigen3** - Álgebra lineal (opcional pero recomendado)
- **Boost** - Utilidades matemáticas (opcional)
- **Matplotlib-cpp** - Visualización (opcional)

### Instalación de dependencias

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install g++ cmake libeigen3-dev libboost-all-dev
```

**macOS:**
```bash
brew install gcc cmake eigen boost
```

**Windows:**
- Instalar [MinGW](https://www.mingw-w64.org/) o [Visual Studio](https://visualstudio.microsoft.com/)
- Descargar [Eigen](https://eigen.tuxfamily.org/)

---

##  Instalación

### Clonar el repositorio

```bash
git clone https://github.com/juareymal-ui/Physics_Numerical_Methods.git
cd Physics_Numerical_Methods
```

### Compilación básica

```bash
# Compilar un programa individual
g++ -std=c++17 -O3 -o euler_method euler_method.cpp

# Ejecutar
./euler_method
```

### Compilación con CMake (Recomendado)

```bash
mkdir build && cd build
cmake ..
make
./bin/nombre_programa
```

---

## 📁 Estructura del Proyecto

```
Physics_Numerical_Methods/
├── README.md
├── CMakeLists.txt
├── include/
│   ├── ode_solvers.h          # Solucionadores de EDOs
│   ├── integration.h          # Métodos de integración
│   ├── physics_simulations.h  # Simulaciones físicas
│   └── numerical_analysis.h   # Herramientas de análisis
├── src/
│   ├── ode_solvers.cpp
│   ├── integration.cpp
│   ├── physics_simulations.cpp
│   └── numerical_analysis.cpp
├── examples/
│   ├── pendulum_simulation.cpp
│   ├── harmonic_oscillator.cpp
│   ├── planetary_motion.cpp
│   └── heat_equation.cpp
├── tests/
│   └── unit_tests.cpp
└── docs/
    └── theory.md
```

---

##  Métodos Implementados

### 1. Ecuaciones Diferenciales Ordinarias (ODEs)

####  Método de Euler

**Teoría:** Aproximación de primer orden para resolver EDOs.

**Ecuación:** 
```
y(n+1) = y(n) + h * f(t(n), y(n))
```

**Implementación:**

```cpp
#include <iostream>
#include <vector>
#include <cmath>

// Método de Euler
std::vector<double> euler_method(
    double (*f)(double, double),  // Función dy/dt = f(t, y)
    double y0,                     // Condición inicial
    double t0,                     // Tiempo inicial
    double tf,                     // Tiempo final
    int n_steps                    // Número de pasos
) {
    double h = (tf - t0) / n_steps;  // Tamaño del paso
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        y[i + 1] = y[i] + h * f(t, y[i]);
        t += h;
    }
    
    return y;
}

// Ejemplo: dy/dt = -2y (decaimiento exponencial)
double exponential_decay(double t, double y) {
    return -2.0 * y;
}

int main() {
    auto solution = euler_method(exponential_decay, 1.0, 0.0, 5.0, 100);
    
    std::cout << "Solución final: " << solution.back() << std::endl;
    return 0;
}
```

**Uso:**
```bash
g++ -std=c++17 -o euler euler_method.cpp
./euler
```

---

####  Método de Runge-Kutta (RK4)

**Teoría:** Método de cuarto orden, más preciso que Euler.

**Implementación:**

```cpp
#include <iostream>
#include <vector>

std::vector<double> runge_kutta_4(
    double (*f)(double, double),
    double y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        double k1 = h * f(t, y[i]);
        double k2 = h * f(t + h/2, y[i] + k1/2);
        double k3 = h * f(t + h/2, y[i] + k2/2);
        double k4 = h * f(t + h, y[i] + k3);
        
        y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        t += h;
    }
    
    return y;
}

// Ejemplo: Oscilador armónico - d²x/dt² = -ω²x
// Convertimos a sistema de primer orden:
// dx/dt = v
// dv/dt = -ω²x

struct State {
    double x;  // Posición
    double v;  // Velocidad
};

State harmonic_oscillator(double t, State s, double omega) {
    return {s.v, -omega * omega * s.x};
}
```

---

### 2. Métodos de Integración Numérica

####  Regla del Trapecio

**Teoría:** Aproxima el área bajo la curva usando trapecios.

```cpp
#include <iostream>
#include <cmath>

double trapezoidal_rule(
    double (*f)(double),  // Función a integrar
    double a,             // Límite inferior
    double b,             // Límite superior
    int n                 // Número de subdivisiones
) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; i++) {
        sum += f(a + i * h);
    }
    
    return h * sum;
}

// Ejemplo: Integrar x²
double square(double x) {
    return x * x;
}

int main() {
    // Integrar x² de 0 a 1 (resultado exacto = 1/3)
    double result = trapezoidal_rule(square, 0.0, 1.0, 1000);
    std::cout << "Integral: " << result << std::endl;
    std::cout << "Error: " << std::abs(result - 1.0/3.0) << std::endl;
    return 0;
}
```

---

####  Regla de Simpson

**Teoría:** Método de tercer orden, más preciso que el trapecio.

```cpp
double simpson_rule(
    double (*f)(double),
    double a,
    double b,
    int n  // Debe ser par
) {
    if (n % 2 != 0) n++;  // Asegurar que n es par
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2 * f(x) : 4 * f(x);
    }
    
    return (h / 3.0) * sum;
}
```

---

####  Cuadratura de Gauss

**Teoría:** Método de alta precisión usando puntos y pesos óptimos.

```cpp
#include <vector>

double gauss_quadrature(
    double (*f)(double),
    double a,
    double b,
    int n_points = 5  // Número de puntos de Gauss
) {
    // Puntos y pesos de Gauss-Legendre (n=5)
    std::vector<double> points = {
        -0.9061798459, -0.5384693101, 0.0, 
         0.5384693101,  0.9061798459
    };
    std::vector<double> weights = {
        0.2369268851, 0.4786286705, 0.5688888889,
        0.4786286705, 0.2369268851
    };
    
    double sum = 0.0;
    double c = (b - a) / 2.0;
    double d = (b + a) / 2.0;
    
    for (int i = 0; i < n_points; i++) {
        sum += weights[i] * f(c * points[i] + d);
    }
    
    return c * sum;
}
```

---

### 3. Simulaciones Físicas

####  Péndulo Simple

**Ecuación:** θ'' + (g/L) sin(θ) = 0

```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class SimplePendulum {
private:
    double g;  // Gravedad
    double L;  // Longitud
    
public:
    SimplePendulum(double gravity = 9.81, double length = 1.0)
        : g(gravity), L(length) {}
    
    // Ecuación del movimiento
    std::pair<double, double> equation(double theta, double omega) {
        double d_theta = omega;
        double d_omega = -(g / L) * std::sin(theta);
        return {d_theta, d_omega};
    }
    
    // Simular con RK4
    void simulate(double theta0, double omega0, double t_max, int steps) {
        double dt = t_max / steps;
        double theta = theta0;
        double omega = omega0;
        
        std::ofstream file("pendulum_data.csv");
        file << "time,theta,omega,energy\n";
        
        for (int i = 0; i <= steps; i++) {
            double t = i * dt;
            double energy = 0.5 * L * L * omega * omega + 
                           g * L * (1 - std::cos(theta));
            
            file << t << "," << theta << "," << omega << "," 
                 << energy << "\n";
            
            // RK4
            auto [k1_theta, k1_omega] = equation(theta, omega);
            auto [k2_theta, k2_omega] = equation(
                theta + 0.5*dt*k1_theta, omega + 0.5*dt*k1_omega);
            auto [k3_theta, k3_omega] = equation(
                theta + 0.5*dt*k2_theta, omega + 0.5*dt*k2_omega);
            auto [k4_theta, k4_omega] = equation(
                theta + dt*k3_theta, omega + dt*k3_omega);
            
            theta += (dt/6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
            omega += (dt/6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);
        }
        
        file.close();
        std::cout << "Simulación completada. Datos guardados en pendulum_data.csv\n";
    }
};

int main() {
    SimplePendulum pendulum(9.81, 1.0);
    
    // Condiciones iniciales: 45° (π/4 rad), velocidad = 0
    double theta0 = M_PI / 4.0;
    double omega0 = 0.0;
    
    pendulum.simulate(theta0, omega0, 10.0, 1000);
    
    return 0;
}
```

**Compilar y ejecutar:**
```bash
g++ -std=c++17 -O3 -o pendulum pendulum.cpp
./pendulum
```

---

####  Movimiento Planetario (Problema de 2 Cuerpos)

**Ecuaciones:** Ley de gravitación de Newton

```cpp
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
    std::array<double, 4> derivatives(const std::array<double, 4>& state) {
        double x = state[0], y = state[1];
        double vx = state[2], vy = state[3];
        
        double r = std::sqrt(x*x + y*y);
        double r3 = r * r * r;
        
        double ax = -G * M * x / r3;
        double ay = -G * M * y / r3;
        
        return {vx, vy, ax, ay};
    }
    
    void simulate(std::array<double, 4> initial, double dt, int steps) {
        std::ofstream file("orbit_data.csv");
        file << "x,y,vx,vy\n";
        
        auto state = initial;
        
        for (int i = 0; i <= steps; i++) {
            file << state[0] << "," << state[1] << "," 
                 << state[2] << "," << state[3] << "\n";
            
            // RK4
            auto k1 = derivatives(state);
            
            std::array<double, 4> temp;
            for (int j = 0; j < 4; j++)
                temp[j] = state[j] + 0.5*dt*k1[j];
            auto k2 = derivatives(temp);
            
            for (int j = 0; j < 4; j++)
                temp[j] = state[j] + 0.5*dt*k2[j];
            auto k3 = derivatives(temp);
            
            for (int j = 0; j < 4; j++)
                temp[j] = state[j] + dt*k3[j];
            auto k4 = derivatives(temp);
            
            for (int j = 0; j < 4; j++)
                state[j] += (dt/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        
        file.close();
    }
};
```

---

####  Ecuación del Calor (1D)

**Ecuación:** ∂u/∂t = α ∂²u/∂x²

```cpp
#include <iostream>
#include <vector>
#include <fstream>

class HeatEquation {
private:
    double alpha;  // Difusividad térmica
    double L;      // Longitud del dominio
    int nx;        // Puntos espaciales
    
public:
    HeatEquation(double diffusivity, double length, int n_points)
        : alpha(diffusivity), L(length), nx(n_points) {}
    
    void solve(double t_max, int nt) {
        double dx = L / (nx - 1);
        double dt = t_max / nt;
        
        // Verificar estabilidad (criterio CFL)
        double stability = alpha * dt / (dx * dx);
        if (stability > 0.5) {
            std::cerr << "Advertencia: Inestabilidad numérica posible\n";
        }
        
        std::vector<double> u(nx), u_new(nx);
        
        // Condición inicial: pulso gaussiano en el centro
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            u[i] = std::exp(-50 * std::pow(x - L/2, 2));
        }
        
        std::ofstream file("heat_data.csv");
        
        // Iterar en el tiempo
        for (int n = 0; n < nt; n++) {
            // Escribir cada 10 pasos
            if (n % 10 == 0) {
                for (int i = 0; i < nx; i++)
                    file << u[i] << (i < nx-1 ? "," : "\n");
            }
            
            // Método de diferencias finitas explícito
            for (int i = 1; i < nx - 1; i++) {
                u_new[i] = u[i] + alpha * dt / (dx * dx) * 
                          (u[i+1] - 2*u[i] + u[i-1]);
            }
            
            // Condiciones de frontera (Dirichlet: temperatura = 0)
            u_new[0] = 0;
            u_new[nx-1] = 0;
            
            u = u_new;
        }
        
        file.close();
        std::cout << "Simulación de ecuación del calor completada\n";
    }
};

int main() {
    HeatEquation heat(0.01, 1.0, 100);
    heat.solve(1.0, 1000);
    return 0;
}
```

---

### 4. Análisis Numérico

####  Método de Newton-Raphson

**Encontrar raíces de f(x) = 0**

```cpp
#include <iostream>
#include <cmath>

double newton_raphson(
    double (*f)(double),       // Función
    double (*df)(double),      // Derivada
    double x0,                 // Estimación inicial
    double tolerance = 1e-6,
    int max_iter = 100
) {
    double x = x0;
    
    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);
        double dfx = df(x);
        
        if (std::abs(dfx) < 1e-10) {
            std::cerr << "Derivada muy pequeña\n";
            return x;
        }
        
        double x_new = x - fx / dfx;
        
        if (std::abs(x_new - x) < tolerance) {
            std::cout << "Convergencia en " << i+1 << " iteraciones\n";
            return x_new;
        }
        
        x = x_new;
    }
    
    std::cerr << "No convergió en " << max_iter << " iteraciones\n";
    return x;
}

// Ejemplo: Encontrar √2 (raíz de x² - 2 = 0)
double f(double x) { return x*x - 2; }
double df(double x) { return 2*x; }

int main() {
    double root = newton_raphson(f, df, 1.0);
    std::cout << "√2 ≈ " << root << std::endl;
    std::cout << "Error: " << std::abs(root - std::sqrt(2)) << std::endl;
    return 0;
}
```

---

##  Ejemplos de Uso

### Ejemplo Completo: Oscilador Armónico Amortiguado

```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class DampedHarmonicOscillator {
private:
    double omega;  // Frecuencia natural
    double gamma;  // Coeficiente de amortiguamiento
    
public:
    DampedHarmonicOscillator(double w, double g) : omega(w), gamma(g) {}
    
    // Ecuación: d²x/dt² + 2γ dx/dt + ω²x = 0
    // Sistema: dx/dt = v, dv/dt = -2γv - ω²x
    std::pair<double, double> equation(double x, double v) {
        return {v, -2*gamma*v - omega*omega*x};
    }
    
    void simulate(double x0, double v0, double t_max, int steps) {
        double dt = t_max / steps;
        double x = x0, v = v0;
        
        std::ofstream file("oscillator.csv");
        file << "time,position,velocity,energy\n";
        
        for (int i = 0; i <= steps; i++) {
            double t = i * dt;
            double energy = 0.5 * v*v + 0.5 * omega*omega * x*x;
            
            file << t << "," << x << "," << v << "," << energy << "\n";
            
            // RK4
            auto [k1x, k1v] = equation(x, v);
            auto [k2x, k2v] = equation(x + 0.5*dt*k1x, v + 0.5*dt*k1v);
            auto [k3x, k3v] = equation(x + 0.5*dt*k2x, v + 0.5*dt*k2v);
            auto [k4x, k4v] = equation(x + dt*k3x, v + dt*k3v);
            
            x += (dt/6.0) * (k1x + 2*k2x + 2*k3x + k4x);
            v += (dt/6.0) * (k1v + 2*k2v + 2*k3v + k4v);
        }
        
        file.close();
    }
};

int main() {
    // Crear oscilador: ω = 2π rad/s, γ = 0.1 s⁻¹
    DampedHarmonicOscillator osc(2*M_PI, 0.1);
    
    // Simular 10 segundos con x0=1m, v0=0
    osc.simulate(1.0, 0.0, 10.0, 1000);
    
    std::cout << "Datos guardados en oscillator.csv\n";
    return 0;
}
```

**Visualizar con Python:**

```python
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('oscillator.csv')

plt.figure(figsize=(12, 4))

plt.subplot(131)
plt.plot(data['time'], data['position'])
plt.xlabel('Tiempo (s)')
plt.ylabel('Posición (m)')
plt.title('Posición vs Tiempo')
plt.grid(True)

plt.subplot(132)
plt.plot(data['position'], data['velocity'])
plt.xlabel('Posición (m)')
plt.ylabel('Velocidad (m/s)')
plt.title('Espacio de Fases')
plt.grid(True)

plt.subplot(133)
plt.plot(data['time'], data['energy'])
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (J)')
plt.title('Energía vs Tiempo')
plt.grid(True)

plt.tight_layout()
plt.savefig('oscillator_analysis.png', dpi=300)
plt.show()
```

---

##  Documentación

### Compilación Optimizada

```bash
# Optimización máxima
g++ -std=c++17 -O3 -march=native -o programa programa.cpp

# Con depuración
g++ -std=c++17 -g -Wall -Wextra -o programa programa.cpp

# Con OpenMP (paralelización)
g++ -std=c++17 -O3 -fopenmp -o programa programa.cpp
```

### Buenas Prácticas

1. **Usa `const` cuando sea posible**
2. **Pasa objetos grandes por referencia**
3. **Reserva memoria con `reserve()`** para vectores
4. **Usa tipos apropiados** (`double` para física, `float` si necesitas velocidad)
5. **Valida entrada de datos**
6. **Maneja errores apropiadamente**

### Análisis de Errores

```cpp
// Calcular error relativo
double relative_error(double exact, double approx) {
    return std::abs((exact - approx) / exact);
}

// Calcular error absoluto
double absolute_error(double exact, double approx) {
    return std::abs(exact - approx);
}
```

---

##  Testing

### Ejemplo de prueba unitaria

```cpp
#include <cassert>
#include <cmath>

void test_trapezoidal_rule() {
    auto f = [](double x) { return x*x; };
    double result = trapezoidal_rule(f, 0.0, 1.0, 1000);
    double exact = 1.0/3.0;
    
    assert(std::abs(result - exact) < 1e-4);
    std::cout << "✓ Regla del trapecio: Test pasado\n";
}

void test_newton_raphson() {
    auto f = [](double x) { return x*x - 2; };
    auto df = [](double x) { return 2*x; };
    double root = newton_raphson(f, df, 1.0);
    
    assert(std::abs(root - std::sqrt(2)) < 1e-6);
    std::cout << "✓ Newton-Raphson: Test pasado\n";
}

int main() {
    test_trapezoidal_rule();
    test_newton_raphson();
    std::cout << "\nTodos los tests pasaron ✓\n";
    return 0;
}
```

---

##  Recursos Adicionales

### Libros Recomendados

- **"Numerical Recipes in C++"** - Press et al.
- **"Computational Physics"** - Mark Newman
- **"An Introduction to Computational Physics"** - Tao Pang

### Tutoriales Online

- [CPlusPlus.com](https://www.cplusplus.com/)
- [LearnCpp.com](https://www.learncpp.com/)
- [Computational Physics with C++](https://www.youtube.com/playlist?list=PLQVvvaa0QuDdxtSl0GQpSmIyFrtRSlkF7)

### Librerías Científicas

- [Eigen](https://eigen.tuxfamily.org/) - Álgebra lineal
- [GSL](https://www.gnu.org/software/gsl/) - GNU Scientific Library
- [Armadillo](http://arma.sourceforge.net/) - C++ linear algebra library

---

##  Contribuir

¡Las contribuciones son bienvenidas! Para contribuir:

1. **Fork** el repositorio
2. Crea una rama: `git checkout -b feature/nueva-funcionalidad`
3. Commit: `git commit -am 'Añadir nueva funcionalidad'`
4. Push: `git push origin feature/nueva-funcionalidad`
5. Abre un **Pull Request**

### Guidelines

- Sigue el estándar **C++17** o superior
- Documenta tu código con comentarios claros
- Incluye ejemplos de uso
- Añade tests cuando sea posible
- Mantén el estilo de código consistente

---

##  Licencia

Este proyecto está bajo la licencia MIT. Ver [LICENSE](LICENSE) para más detalles.

---

##  Autor

**Juan Maldonado**

- GitHub: [@juareymal-ui](https://github.com/juareymal-ui)
- LinkedIn: [Juan Maldonado](https://www.linkedin.com/in/juan-de-jes%C3%BAs-reyes-maldonado-b7475a2b9/)
- Email: lessjuareymal@gmail.com

---

##  Roadmap

### Versión Actual (v1.0)
- ✅ Métodos básicos de EDOs (Euler, RK4)
- ✅ Integración numérica (Trapecio, Simpson, Gauss)
- ✅ Simulaciones físicas básicas
- ✅ Análisis numérico fundamental

### Próximas Características (v2.0)

- [ ] **Métodos avanzados de EDOs**
  - Adams-Bashforth
  - Adams-Moulton
  - Métodos implícitos (BDF)
  
- [ ] **EDPs (Ecuaciones en Derivadas Parciales)**
  - Ecuación de onda 2D
  - Ecuación de Laplace
  - Método de elementos finitos
  
- [ ] **Optimización**
  - Gradient descent
  - Método de Newton multivariable
  - Algoritmos genéticos
  
- [ ] **Álgebra Lineal Numérica**
  - Factorización LU
  - Método QR
  - Valores y vectores propios
  
- [ ] **Monte Carlo**
  - Integración Monte Carlo
  - Simulación de procesos estocásticos
  - MCMC (Markov Chain Monte Carlo)

---

## 🔬 Ejemplos Avanzados

### Simulación N-Cuerpos (Sistema Solar)

```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

struct Body {
    double mass;
    double x, y, z;      // Posición
    double vx, vy, vz;   // Velocidad
    double ax, ay, az;   // Aceleración
    
    Body(double m, double px, double py, double pz,
         double velx, double vely, double velz)
        : mass(m), x(px), y(py), z(pz),
          vx(velx), vy(vely), vz(velz),
          ax(0), ay(0), az(0) {}
};

class NBodySimulation {
private:
    std::vector<Body> bodies;
    double G;  // Constante gravitacional
    
public:
    NBodySimulation() : G(6.67430e-11) {}
    
    void add_body(const Body& body) {
        bodies.push_back(body);
    }
    
    void compute_forces() {
        // Reiniciar aceleraciones
        for (auto& body : bodies) {
            body.ax = body.ay = body.az = 0.0;
        }
        
        // Calcular fuerzas entre todos los pares
        for (size_t i = 0; i < bodies.size(); i++) {
            for (size_t j = i + 1; j < bodies.size(); j++) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                double f = G / (r * r * r);
                
                // Fuerza sobre el cuerpo i
                bodies[i].ax += f * bodies[j].mass * dx;
                bodies[i].ay += f * bodies[j].mass * dy;
                bodies[i].az += f * bodies[j].mass * dz;
                
                // Fuerza sobre el cuerpo j (tercera ley de Newton)
                bodies[j].ax -= f * bodies[i].mass * dx;
                bodies[j].ay -= f * bodies[i].mass * dy;
                bodies[j].az -= f * bodies[i].mass * dz;
            }
        }
    }
    
    void update_positions(double dt) {
        // Método de Verlet (más estable para órbitas)
        for (auto& body : bodies) {
            // Actualizar velocidades (medio paso)
            body.vx += 0.5 * body.ax * dt;
            body.vy += 0.5 * body.ay * dt;
            body.vz += 0.5 * body.az * dt;
            
            // Actualizar posiciones
            body.x += body.vx * dt;
            body.y += body.vy * dt;
            body.z += body.vz * dt;
        }
        
        // Recalcular fuerzas
        compute_forces();
        
        // Actualizar velocidades (segundo medio paso)
        for (auto& body : bodies) {
            body.vx += 0.5 * body.ax * dt;
            body.vy += 0.5 * body.ay * dt;
            body.vz += 0.5 * body.az * dt;
        }
    }
    
    void simulate(double t_max, double dt, const std::string& filename) {
        std::ofstream file(filename);
        file << "time,body_id,x,y,z\n";
        
        int steps = static_cast<int>(t_max / dt);
        
        for (int step = 0; step <= steps; step++) {
            double t = step * dt;
            
            // Guardar posiciones cada 100 pasos
            if (step % 100 == 0) {
                for (size_t i = 0; i < bodies.size(); i++) {
                    file << t << "," << i << ","
                         << bodies[i].x << ","
                         << bodies[i].y << ","
                         << bodies[i].z << "\n";
                }
            }
            
            update_positions(dt);
        }
        
        file.close();
        std::cout << "Simulación completada. Datos en " << filename << "\n";
    }
    
    double total_energy() const {
        double kinetic = 0.0, potential = 0.0;
        
        // Energía cinética
        for (const auto& body : bodies) {
            kinetic += 0.5 * body.mass * 
                      (body.vx*body.vx + body.vy*body.vy + body.vz*body.vz);
        }
        
        // Energía potencial
        for (size_t i = 0; i < bodies.size(); i++) {
            for (size_t j = i + 1; j < bodies.size(); j++) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dz = bodies[j].z - bodies[i].z;
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                potential -= G * bodies[i].mass * bodies[j].mass / r;
            }
        }
        
        return kinetic + potential;
    }
};

int main() {
    NBodySimulation sim;
    
    // Sistema Sol-Tierra-Luna (valores aproximados)
    // Masas en kg, distancias en m, velocidades en m/s
    
    // Sol (en el origen)
    sim.add_body(Body(1.989e30, 0, 0, 0, 0, 0, 0));
    
    // Tierra (1 UA = 1.496e11 m)
    sim.add_body(Body(5.972e24, 1.496e11, 0, 0, 0, 29783, 0));
    
    // Luna (384,400 km de la Tierra)
    sim.add_body(Body(7.342e22, 1.496e11 + 3.844e8, 0, 0, 0, 29783 + 1022, 0));
    
    // Simular 1 año (31,536,000 segundos)
    double year = 365.25 * 24 * 3600;
    double dt = 3600;  // 1 hora
    
    std::cout << "Energía inicial: " << sim.total_energy() << " J\n";
    
    sim.simulate(year, dt, "solar_system.csv");
    
    std::cout << "Energía final: " << sim.total_energy() << " J\n";
    
    return 0;
}
```

---

### Ecuación de Schrödinger (1D) - Partícula en una Caja

```cpp
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>

using Complex = std::complex<double>;

class QuantumParticle {
private:
    int nx;           // Puntos espaciales
    double L;         // Longitud de la caja
    double dx;        // Paso espacial
    double hbar;      // Constante de Planck reducida
    double m;         // Masa de la partícula
    
    std::vector<Complex> psi;     // Función de onda
    std::vector<double> V;         // Potencial
    
public:
    QuantumParticle(int n_points, double length, double mass)
        : nx(n_points), L(length), m(mass), hbar(1.054571817e-34) {
        
        dx = L / (nx - 1);
        psi.resize(nx);
        V.resize(nx);
        
        // Inicializar potencial (caja infinita)
        for (int i = 0; i < nx; i++) {
            V[i] = 0.0;  // Potencial = 0 dentro de la caja
        }
    }
    
    void set_initial_state(int n) {
        // Estado propio n-ésimo de una caja infinita
        double A = std::sqrt(2.0 / L);
        
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            psi[i] = A * std::sin(n * M_PI * x / L);
        }
    }
    
    void set_gaussian_wavepacket(double x0, double sigma, double k0) {
        // Paquete de onda gaussiano
        double A = std::pow(1.0 / (sigma * std::sqrt(M_PI)), 0.5);
        
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            double gauss = std::exp(-0.5 * std::pow((x - x0) / sigma, 2));
            psi[i] = A * gauss * std::exp(Complex(0, k0 * x));
        }
        
        normalize();
    }
    
    void normalize() {
        double norm = 0.0;
        for (const auto& val : psi) {
            norm += std::norm(val) * dx;
        }
        norm = std::sqrt(norm);
        
        for (auto& val : psi) {
            val /= norm;
        }
    }
    
    // Evolución temporal usando método de Crank-Nicolson
    void evolve(double dt, int steps, const std::string& filename) {
        std::ofstream file(filename);
        file << "x,real,imag,probability\n";
        
        Complex I(0, 1);
        double alpha = hbar * dt / (2.0 * m * dx * dx);
        
        for (int step = 0; step <= steps; step++) {
            // Guardar cada 10 pasos
            if (step % 10 == 0) {
                for (int i = 0; i < nx; i++) {
                    double x = i * dx;
                    double prob = std::norm(psi[i]);
                    file << x << "," << psi[i].real() << "," 
                         << psi[i].imag() << "," << prob << "\n";
                }
                file << "\n";  // Separador para diferentes tiempos
            }
            
            // Método de Crank-Nicolson (simplificado para potencial = 0)
            std::vector<Complex> psi_new = psi;
            
            for (int i = 1; i < nx - 1; i++) {
                psi_new[i] = psi[i] + I * alpha * 
                            (psi[i+1] - 2.0*psi[i] + psi[i-1]);
            }
            
            // Condiciones de frontera (ψ = 0 en los bordes)
            psi_new[0] = 0;
            psi_new[nx-1] = 0;
            
            psi = psi_new;
            normalize();
        }
        
        file.close();
    }
    
    double expectation_position() const {
        double x_avg = 0.0;
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            x_avg += x * std::norm(psi[i]) * dx;
        }
        return x_avg;
    }
    
    double expectation_momentum() const {
        // <p> usando diferencias finitas
        double p_avg = 0.0;
        for (int i = 1; i < nx - 1; i++) {
            Complex dpsi_dx = (psi[i+1] - psi[i-1]) / (2.0 * dx);
            p_avg += std::real(std::conj(psi[i]) * (-I * hbar * dpsi_dx)) * dx;
        }
        return p_avg;
    }
};

int main() {
    // Crear partícula en caja de 1 nm
    QuantumParticle particle(1000, 1e-9, 9.109e-31);  // electrón
    
    // Estado inicial: paquete gaussiano
    particle.set_gaussian_wavepacket(0.5e-9, 0.1e-9, 1e10);
    
    std::cout << "Posición inicial: " << particle.expectation_position() << " m\n";
    std::cout << "Momento inicial: " << particle.expectation_momentum() << " kg⋅m/s\n";
    
    // Evolucionar 1 picosegundo
    particle.evolve(1e-17, 1000, "quantum_evolution.csv");
    
    std::cout << "Posición final: " << particle.expectation_position() << " m\n";
    std::cout << "Momento final: " << particle.expectation_momentum() << " kg⋅m/s\n";
    
    return 0;
}
```

---

### Método de Monte Carlo - Integración

```cpp
#include <iostream>
#include <random>
#include <cmath>
#include <functional>

class MonteCarlo {
private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> uniform;
    
public:
    MonteCarlo() : rng(std::random_device{}()), uniform(0.0, 1.0) {}
    
    // Integración Monte Carlo en 1D
    double integrate_1d(
        std::function<double(double)> f,
        double a, double b,
        int n_samples
    ) {
        double sum = 0.0;
        
        for (int i = 0; i < n_samples; i++) {
            double x = a + (b - a) * uniform(rng);
            sum += f(x);
        }
        
        return (b - a) * sum / n_samples;
    }
    
    // Integración Monte Carlo en múltiples dimensiones
    double integrate_nd(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& a,
        const std::vector<double>& b,
        int n_samples
    ) {
        int dim = a.size();
        double volume = 1.0;
        
        for (int d = 0; d < dim; d++) {
            volume *= (b[d] - a[d]);
        }
        
        double sum = 0.0;
        std::vector<double> x(dim);
        
        for (int i = 0; i < n_samples; i++) {
            for (int d = 0; d < dim; d++) {
                x[d] = a[d] + (b[d] - a[d]) * uniform(rng);
            }
            sum += f(x);
        }
        
        return volume * sum / n_samples;
    }
    
    // Estimar π usando Monte Carlo
    double estimate_pi(int n_samples) {
        int inside_circle = 0;
        
        for (int i = 0; i < n_samples; i++) {
            double x = uniform(rng);
            double y = uniform(rng);
            
            if (x*x + y*y <= 1.0) {
                inside_circle++;
            }
        }
        
        return 4.0 * inside_circle / n_samples;
    }
    
    // Calcular error estándar
    double standard_error(
        std::function<double(double)> f,
        double a, double b,
        int n_samples
    ) {
        double sum = 0.0;
        double sum_sq = 0.0;
        
        for (int i = 0; i < n_samples; i++) {
            double x = a + (b - a) * uniform(rng);
            double fx = f(x);
            sum += fx;
            sum_sq += fx * fx;
        }
        
        double mean = sum / n_samples;
        double variance = (sum_sq / n_samples) - mean * mean;
        
        return (b - a) * std::sqrt(variance / n_samples);
    }
};

int main() {
    MonteCarlo mc;
    
    // Ejemplo 1: Integrar x² de 0 a 1
    auto f1 = [](double x) { return x * x; };
    double result1 = mc.integrate_1d(f1, 0.0, 1.0, 1000000);
    double error1 = mc.standard_error(f1, 0.0, 1.0, 1000000);
    
    std::cout << "∫₀¹ x² dx = " << result1 << " ± " << error1 << std::endl;
    std::cout << "Valor exacto: " << 1.0/3.0 << std::endl;
    std::cout << "Error: " << std::abs(result1 - 1.0/3.0) << "\n\n";
    
    // Ejemplo 2: Estimar π
    double pi_estimate = mc.estimate_pi(10000000);
    std::cout << "Estimación de π: " << pi_estimate << std::endl;
    std::cout << "Valor real: " << M_PI << std::endl;
    std::cout << "Error: " << std::abs(pi_estimate - M_PI) << "\n\n";
    
    // Ejemplo 3: Integral multidimensional
    // Volumen de una esfera unitaria en 3D
    auto sphere = [](const std::vector<double>& x) {
        double r_sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return (r_sq <= 1.0) ? 1.0 : 0.0;
    };
    
    std::vector<double> lower = {-1, -1, -1};
    std::vector<double> upper = {1, 1, 1};
    
    double volume = mc.integrate_nd(sphere, lower, upper, 1000000);
    std::cout << "Volumen de esfera unitaria: " << volume << std::endl;
    std::cout << "Valor exacto (4π/3): " << 4.0 * M_PI / 3.0 << std::endl;
    
    return 0;
}
```

---

##  Proyectos Propuestos

### Nivel Principiante

1. **Calculadora de Trayectorias**
   - Implementar el movimiento parabólico
   - Incluir resistencia del aire
   - Graficar trayectorias

2. **Simulador de Resortes**
   - Sistema masa-resorte simple
   - Sistema de resortes acoplados
   - Análisis de frecuencias

3. **Integrador Universal**
   - Comparar diferentes métodos
   - Análisis de convergencia
   - Visualización de errores

### Nivel Intermedio

4. **Simulador de Fluidos (Lattice Boltzmann)**
   - Ecuación de Navier-Stokes simplificada
   - Visualización en tiempo real
   - Diferentes condiciones de frontera

5. **Cadenas de Osciladores**
   - Modos normales de vibración
   - Transformada de Fourier
   - Análisis espectral

6. **Problema de Tres Cuerpos**
   - Órbitas caóticas
   - Puntos de Lagrange
   - Conservación de energía

### Nivel Avanzado

7. **Simulación de Spin Glasses (Modelo de Ising)**
   - Algoritmo de Metropolis
   - Transiciones de fase
   - Análisis estadístico

8. **Propagación de Ondas Electromagnéticas**
   - Ecuaciones de Maxwell (FDTD)
   - Reflexión y refracción
   - Visualización 3D

9. **Dinámica Molecular**
   - Potenciales de Lennard-Jones
   - Ensemble termodinámicos
   - Cálculo de propiedades

---

## 💡 Tips y Trucos

### Optimización de Código

```cpp
// ❌ Evitar: Cálculos repetidos
for (int i = 0; i < n; i++) {
    double result = expensive_function(x) * i;
}

// ✅ Mejor: Calcular una vez
double cached = expensive_function(x);
for (int i = 0; i < n; i++) {
    double result = cached * i;
}

// ❌ Evitar: Realocación constante
std::vector<double> data;
for (int i = 0; i < n; i++) {
    data.push_back(i);  // Puede realocar memoria
}

// ✅ Mejor: Reservar memoria
std::vector<double> data;
data.reserve(n);
for (int i = 0; i < n; i++) {
    data.push_back(i);
}

// ❌ Evitar: Copias innecesarias
void process(std::vector<double> data) {  // Copia el vector
    // ...
}

// ✅ Mejor: Pasar por referencia
void process(const std::vector<double>& data) {  // No copia
    // ...
}
```

### Debugging

```cpp
// Macro útil para debugging
#define DEBUG_PRINT(var) \
    std::cout << #var << " = " << (var) << std::endl

// Uso
double x = 3.14;
DEBUG_PRINT(x);  // Imprime: x = 3.14

// Verificar valores NaN o infinitos
#include <cmath>

if (std::isnan(value)) {
    std::cerr << "Error: NaN detectado\n";
}

if (std::isinf(value)) {
    std::cerr << "Error: Infinito detectado\n";
}
```

---

##  Benchmarking

### Medir tiempo de ejecución

```cpp
#include <chrono>

auto start = std::chrono::high_resolution_clock::now();

// Código a medir
for (int i = 0; i < 1000000; i++) {
    // ... operaciones ...
}

auto end = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

std::cout << "Tiempo: " << duration.count() << " ms\n";
```

### Comparar métodos

```cpp
void benchmark_integration_methods() {
    auto f = [](double x) { return std::sin(x); };
    
    auto start = std::chrono::high_resolution_clock::now();
    double result1 = trapezoidal_rule(f, 0, M_PI, 10000);
    auto end = std::chrono::high_resolution_clock::now();
    auto time1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    start = std::chrono::high_resolution_clock::now();
    double result2 = simpson_rule(f, 0, M_PI, 10000);
    end = std::chrono::high_resolution_clock::now();
    auto time2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "Trapecio: " << result1 << " (" << time1.count() << " μs)\n";
    std::cout << "Simpson: " << result2 << " (" << time2.count() << " μs)\n";
}
```

---

##  FAQ

**Q: ¿Por qué usar C++ en lugar de Python para física computacional?**

A: C++ ofrece:
- 10-100x más velocidad que Python
- Mejor control de memoria
- Ideal para simulaciones grandes o en tiempo real
- Aprendizaje de programación de bajo nivel

**Q: ¿Necesito conocimientos avanzados de C++?**

A: No. Con conocimientos básicos de:
- Variables y tipos
- Bucles y condicionales
- Funciones
- Arrays/vectores

Ya puedes empezar. Los ejemplos están diseñados para ser educativos.

**Q: ¿Cómo visualizo los resultados?**

A: Usa Python con matplotlib:
```bash
python3 visualize.py datos.csv
```

O integra matplotlib-cpp directamente en C++.

**Q: ¿Los métodos son precisos para investigación real?**

A: Los métodos básicos (Euler, RK4) son educativos. Para investigación seria, usa librerías especializadas como:
- SUNDIALS (EDOs)
- PETSc (EDPs)
- LAPACK (álgebra lineal)

---

##  Showcase

Comparte tus proyectos basados en este repositorio:

- Crea un issue con tag `showcase`
- Incluye capturas o videos
- Describe tu implementación

Los mejores proyectos serán destacados aquí.

---

##  Agradecimientos

- Comunidad de [r/cpp](https://reddit.com/r/cpp)
- [Computational Physics](http://www-personal.umich.edu/~mejn/cp/) - Mark Newman
- [Numerical Recipes](http://numerical.recipes/)

---

##  Soporte

¿Tienes preguntas? ¿Encontraste un bug?

-  Email: lessjuareymal@gmail.com

---

<div align="center">

**Si este proyecto te fue útil, considera darle una estrella**

[![Star](https://img.shields.io/github/stars/juareymal-ui/Physics_Numerical_Methods?style=social)](https://github.com/juareymal-ui/Portfolio-Roadmap-Physics-Numerical-Methods.)
Made with by [Juan Maldonado](https://github.com/juareymal-ui)

</div>
