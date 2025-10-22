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
