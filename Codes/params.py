import numpy as np

# Options
PLOT_POSITIONS = True
PLOT_FIGURE = True
PLOT_VIDEO = False # If set to True, will take a huge amount of time
SAVE_FIGURE = False
SAVE_RESULTS = True
final_time = 20.0 # Try 10 or 15 if want to check/debug and save time

# setting up test params
n_elem = 50
start = np.zeros((3,))
direction = np.array([-1.0, 0.0, 0.0])
normal = np.array([0.0, 0.0, 1.0])
base_length = 0.13
base_radius = 0.01416
base_area = np.pi*(base_radius**2)
density = 1180.0
nu = 10 # Try 7 instead for Task Space Control
youngs_modulus = 3.79e6
poisson_ratio = 0.5
shear_modulus = youngs_modulus/(2*(1+poisson_ratio))

dl = base_length / n_elem
dt = (0.002) * dl
total_steps = int(final_time / dt)
step_skip = (
    60
    if PLOT_VIDEO
    else (int(total_steps / 10) if PLOT_FIGURE else int(total_steps / 200))
)