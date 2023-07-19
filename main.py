import numpy as np
  # dx Spatial step in ft
  # dt Time step in hours
def heatmap(dx, dt):
    # Given parameters
    Ti = 100.0  # Initial temperature in °F
    Ts = 300.0  # Surface temperature in °F
    alpha = 0.1  # Diffusivity in ft^2/hour
    tf = 0.1  # Final time in hours

    # Calculate the number of spatial and temporal points
    N = int(1 / dx) + 1
    M = int(tf / dt) + 1

    # Initialize temperature grid
    u = np.zeros((N, M))
    u[:, 0] = Ti  # Initial condition

    # Update boundary conditions
    u[0, :] = Ts
    u[-1, :] = Ts

    # Solve the heat equation iteratively
    for n in range(M - 1):
        for i in range(1, N - 1):
            u[i, n + 1] = u[i, n] + alpha * dt / dx ** 2 * (u[i + 1, n] - 2 * u[i, n] + u[i - 1, n])

    # Extract the temperature distribution at the final time
    temperature_final = u[:, -1]

    # Print the approximate temperature distribution at the final time
    print("Approximate Temperature Distribution at t = 0.1 hour:")
    for i, temp in enumerate(temperature_final):
        print("Position x = {:.2f} ft: Temperature = {:.2f} °F".format(i * dx, temp))

print("Letra A")
heatmap(0.05, 0.01)
print("-------------")
print("Letra B")
heatmap(0.05, 0.05)


def exact_solution(x, t, Ti, Ts, L, alpha, m_max):
    exact_temp = Ts + 2 * (Ti + Ts) * np.sum((1 - (-1) ** np.arange(1, m_max+1)) / (np.pi * np.arange(1, m_max+1)) * np.exp(-(np.pi * np.arange(1, m_max+1) / L) ** 2 * alpha * t) * np.sin(np.pi * np.arange(1, m_max+1) * x / L))
    return exact_temp

# Given parameters
Ti = 100.0  # Initial temperature in °F
Ts = 300.0  # Surface temperature in °F
alpha = 0.1  # Diffusivity in ft^2/hour
L = 1.0  # Length of the wall in ft
t = 0.1  # Time in hours
x_vals = np.linspace(0, L, 21)  # Positions within the wall

# Calculate the exact solution for each position
exact_temps = []
for x in x_vals:
    exact_temp = exact_solution(x, t, Ti, Ts, L, alpha, 100)
    exact_temps.append(exact_temp)

# Print the exact solution
print("Exact Temperature Distribution at t = 0.1 hour:")
for x, temp in zip(x_vals, exact_temps):
    print("Position x = {:.2f} ft: Temperature = {:.2f} °F".format(x, temp))
