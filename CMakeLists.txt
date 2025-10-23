cmake_minimum_required(VERSION 3.15)
project(PhysicsNumericalMethods VERSION 1.0 LANGUAGES CXX)

# C++17 estándar
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Opciones de compilación
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Directorios de include
include_directories(${PROJECT_SOURCE_DIR}/include)

# Librería principal
add_library(numerical_methods
    src/ode_solvers.cpp
    src/integration.cpp
    src/physics_simulations.cpp
    src/numerical_analysis.cpp
)

# Ejemplos
add_executable(pendulum_simulation examples/pendulum_simulation.cpp)
target_link_libraries(pendulum_simulation numerical_methods)

add_executable(harmonic_oscillator examples/harmonic_oscillator.cpp)
target_link_libraries(harmonic_oscillator numerical_methods)

add_executable(planetary_motion examples/planetary_motion.cpp)
target_link_libraries(planetary_motion numerical_methods)

add_executable(heat_equation examples/heat_equation.cpp)
target_link_libraries(heat_equation numerical_methods)

# Tests
add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests numerical_methods)

# Directorios de salida
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
