import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # enables 3D plotting

# Load CSV
data = np.loadtxt("wave3d_snapshot.csv", delimiter=",", skiprows=0)
r = data[:,0]
theta = data[:,1]
phi = data[:,2]
u = data[:,3]

# Get grid dimensions
Nr = len(np.unique(r))
Ntheta = len(np.unique(theta))
Nphi = len(np.unique(phi))

# Reshape to 3D arrays
u3d = u.reshape((Nr, Ntheta, Nphi))
r3d = r.reshape((Nr, Ntheta, Nphi))
theta3d = theta.reshape((Nr, Ntheta, Nphi))
phi3d = phi.reshape((Nr, Ntheta, Nphi))

# Convert spherical to Cartesian coordinates
X = r3d * np.sin(theta3d) * np.cos(phi3d)
Y = r3d * np.sin(theta3d) * np.sin(phi3d)
Z = r3d * np.cos(theta3d)

# For simplicity, plot a single phi slice (phi=0)
phi_index = 0
X_slice = X[:,:,phi_index]
Y_slice = Y[:,:,phi_index]
Z_slice = Z[:,:,phi_index]
U_slice = u3d[:,:,phi_index]

# Create figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Surface plot
surf = [ax.plot_surface(X_slice, Y_slice, Z_slice, facecolors=plt.cm.viridis(U_slice),
                        rstride=1, cstride=1, linewidth=0, antialiased=False)]

ax.set_zlim(-1, 1)

# Animation function
def update(frame):
    ax.clear()
    U_frame = u3d[:,:,frame % Nphi]  # cycle over phi slices as example
    surf = ax.plot_surface(X_slice, Y_slice, Z_slice,
                           facecolors=plt.cm.viridis(U_frame),
                           rstride=1, cstride=1, linewidth=0, antialiased=False)
    ax.set_zlim(-1,1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"Time step {frame}")
    return surf,

ani = animation.FuncAnimation(fig, update, frames=Nphi, interval=100, blit=False)

plt.show()
