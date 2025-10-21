import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("wave_spherical.csv")
r = [i/200 for i in range(201)]

plt.figure()
for i in range(0, len(df), 50):  # plot every 50th time step
    plt.plot(r, df.iloc[i,1:], label=f"t={df.iloc[i,0]:.3f}")
plt.xlabel("r")
plt.ylabel("u(r,t)")
plt.title("Spherically symmetric wave")
plt.legend()
plt.show()
