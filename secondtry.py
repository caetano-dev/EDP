import numpy as np

# Constants
Ti = 100  # Initial temperature (°F)
Ts = 300  # Surface temperature (°F)
alpha = 0.1  # Diffusivity (ft^2/hour)
tf = 0.1  # Final time (hour)
dx = 0.05  # Spatial step size (ft)
#dt = 0.01  # Time step size (hour)
dt = 0.05  # Time step size (hour)

# Calculate number of spatial and time steps
nx = int(1 / dx) + 1
nt = int(tf / dt) + 1

# Initialize temperature array
u = np.full((nt, nx), Ti)

# Apply boundary conditions
u[:, 0] = Ts
u[:, -1] = Ts

# Perform time-stepping iterations
for n in range(1, nt):
    for i in range(1, nx - 1):
        u[n, i] = u[n-1, i] + alpha * dt / dx**2 * (u[n-1, i+1] - 2*u[n-1, i] + u[n-1, i-1])

# Print the approximate solution at the final time
print(f"Approximate solution at t = {tf} hour:")
for i in range(nx):
    print(f"x = {i * dx:.2f} ft: T = {u[-1, i]:.2f} °F")
