import sys
sys.path.append("../../")
import numpy as np
import math
# import control
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb
import matplotlib.animation as manimation
from tqdm import tqdm
# import wrappers
from elastica.wrappers import BaseSystemCollection, Constraints, Forcing, CallBacks, Connections
from elastica._linalg import _batch_product_i_k_to_ik, _batch_matvec

# import rod class and forces to be applied
from elastica.rod.cosserat_rod import CosseratRod
from elastica.external_forces import GravityForces, UniformForces, UniformTorques, NoForces

from elastica.boundary_conditions import ConstraintBase, OneEndFixedBC
from elastica.joint import FixedJoint
# import timestepping functions
from elastica.timestepper.symplectic_steppers import PositionVerlet
from elastica.timestepper import integrate

# import call back functions
from elastica.callback_functions import CallBackBaseClass
from collections import defaultdict
# Import parameter values
from params import *

# Creating a class with all the wrappers needed for simulation
class SwingingFlexiblePendulumSimulator(BaseSystemCollection, Connections, Constraints, Forcing, CallBacks):
    pass



# For 10 elements, the prefac is  0.0007
pendulum_sim = SwingingFlexiblePendulumSimulator()


rod_1 = CosseratRod.straight_rod(
    n_elem,
    start,
    direction,
    normal,
    base_length,
    base_radius,
    density,
    nu,
    youngs_modulus,
    shear_modulus = youngs_modulus/(2*(1+poisson_ratio)),
)
pendulum_sim.append(rod_1)

# ------ One End Fixed BC ------
pendulum_sim.constrain(rod_1).using(
    OneEndFixedBC, 
    constrained_position_idx=(0,), 
    constrained_director_idx=(0,)
)

# ------ Creating Rod 2 and connecting with Rod 1 in Series --------
start_2 = np.zeros((3,))
start_2[0] = -base_length
rod_2 = CosseratRod.straight_rod(
    n_elem,
    start_2,
    direction,
    normal,
    base_length,
    base_radius,
    density,
    nu,
    youngs_modulus,
    shear_modulus = youngs_modulus/(2*(1+poisson_ratio)),
)
pendulum_sim.append(rod_2)
pendulum_sim.connect(rod_1, rod_2, first_connect_idx=-1, second_connect_idx=0).using(
    FixedJoint, 
    k = 1e5,
    nu = 0,
    kt = 5e3
)

# ------ Creating Rod 3 and connecting with Rod 2 in Series --------
start_3 = np.zeros((3,))
start_3[0] = -2*base_length
rod_3 = CosseratRod.straight_rod(
    n_elem,
    start_3,
    direction,
    normal,
    base_length,
    base_radius,
    density,
    nu,
    youngs_modulus,
    shear_modulus = youngs_modulus/(2*(1+poisson_ratio)),
)
pendulum_sim.append(rod_3)
pendulum_sim.connect(rod_2, rod_3, first_connect_idx=-1, second_connect_idx=0).using(
    FixedJoint, 
    k = 1e5,
    nu = 0,
    kt = 5e3
)

# ------------- FEEDBACK FUNCTION ------------
# -------- JOINT SPACE PD CONTROL ------------
def JointFeedBack(system, time, rod_number):
    tang = system.tangents.copy()
    s = system.lengths.copy()
    curv = system.kappa.copy()
    ang_vel = system.omega_collection.copy()
    q_0 = s[:-1].dot(curv[1,:].T)
    # print(q_0)
    theta = np.arctan2(-1.0*tang[2, 0], tang[0, 0])
    if(theta < -3.055): # 5 degrees threshold
        theta = -theta
    theta_dot = ang_vel[1, 0]
    # print(theta_dot)
    q_0_dot = ang_vel[1, -1] - ang_vel[1, 0]
    # -------- If nu = 10 ---------
    # Kp_q_0 = 1.7275
    # -------- If nu = 7 ----------
    Kp_q_0 = 3.3275
    # Kd_q_0 = 0.001
    Kd_q_0 = 0.0075
    if rod_number == 1:
        q_0_des = np.pi/8*math.sin(2*np.pi/2.5*time)
    elif rod_number == 2:
        q_0_des = np.pi/4*math.cos(2*np.pi/5*time)
    elif rod_number == 3:
        q_0_des = np.pi/2*math.sin(2*np.pi/4*time + np.pi/9)    
    theta_des = np.pi
    theta_dot_des = 0
    k = 0.942 # Stiffness (N-m/rad)
    beta = 0.07 # Damping (N-m-s/rad)
    torque = Kp_q_0*(q_0_des - q_0) - Kd_q_0*q_0_dot
    return torque + k*(q_0)

# ---------------- Point Torques ------------------------
class StartPointTorque(NoForces):
    def __init__(self, torque, direction=np.array([0.0, 0.0, 0.0]), rod_number=1):
        """
        Parameters
        ----------
        torque: float
            Torque magnitude applied to a rod-like object.
        direction: numpy.ndarray
            1D (dim) array containing data with 'float' type.
            Direction in which torque applied.
        n: Node number 
        """
        super(StartPointTorque, self).__init__()
        self.torque = direction
        self.rod_number = rod_number

    def apply_torques(self, system, time: np.float64 = 0.0):
        system.external_torques[..., 0] += JointFeedBack(system, time, self.rod_number)*self.torque
        
class EndPointTorque(NoForces):
    def __init__(self, torque, direction=np.array([0.0, 0.0, 0.0]), rod_number=1):
        """
        Parameters
        ----------
        torque: float
            Torque magnitude applied to a rod-like object.
        direction: numpy.ndarray
            1D (dim) array containing data with 'float' type.
            Direction in which torque applied.
        n: Node number 
        """
        super(EndPointTorque, self).__init__()
        self.torque = direction
        self.rod_number = rod_number

    def apply_torques(self, system, time: np.float64 = 0.0):
        system.external_torques[..., -1] += JointFeedBack(system, time, self.rod_number)*self.torque
        
# ----------- ENDPOINT TORQUES --------

pendulum_sim.add_forcing_to(rod_1).using(
    EndPointTorque,
    0.125,
    direction=np.array([0.0, 1.0, 0.0]),
    rod_number=1
)

pendulum_sim.add_forcing_to(rod_2).using(
    StartPointTorque,
    -0.025,
    direction=np.array([0.0, -1.0, 0.0]),
    rod_number=2
)

pendulum_sim.add_forcing_to(rod_2).using(
    EndPointTorque,
    -0.025,
    direction=np.array([0.0, 1.0, 0.0]),
    rod_number=2
)

pendulum_sim.add_forcing_to(rod_3).using(
    StartPointTorque,
    0.025,
    direction=np.array([0.0, -1.0, 0.0]),
    rod_number=3
)

pendulum_sim.add_forcing_to(rod_3).using(
    EndPointTorque,
    0.025,
    direction=np.array([0.0, 1.0, 0.0]),
    rod_number=3
)

# -------------- UNCOMMENT TO ADD GRAVITATIONAL FORCES -----------
# gravitational_acc = -9.80665
# pendulum_sim.add_forcing_to(rod_1).using(
#     GravityForces, acc_gravity=np.array([gravitational_acc, 0.0, 0.0])
# )
# pendulum_sim.add_forcing_to(rod_2).using(
#     GravityForces, acc_gravity=np.array([gravitational_acc, 0.0, 0.0])
# )
# pendulum_sim.add_forcing_to(rod_3).using(
#     GravityForces, acc_gravity=np.array([gravitational_acc, 0.0, 0.0])
# )

# ---------------- ADD CALLBACK FUNCTION -------------
class PendulumCallBack(CallBackBaseClass):
    """
    Call back function for continuum snake
    """

    def __init__(self, step_skip: int, callback_params: dict):
        CallBackBaseClass.__init__(self)
        self.every = step_skip
        self.callback_params = callback_params

    def make_callback(self, system, time, current_step: int):
        if current_step % self.every == 0:
            self.callback_params["time"].append(time)
            self.callback_params["position"].append(system.position_collection.copy())
            self.callback_params["directors"].append(system.director_collection.copy())
            self.callback_params["curvature"].append(system.kappa.copy())
            self.callback_params["lengths"].append(system.lengths.copy())
            self.callback_params["tangents"].append(system.tangents.copy())
            if time > 0.0:
                self.callback_params["internal_stress"].append(
                    system.internal_stress.copy()
                )
                self.callback_params["internal_couple"].append(
                    system.internal_couple.copy()
                )
        return




print("Total steps", total_steps)
recorded_history_1 = defaultdict(list)
recorded_history_2 = defaultdict(list)
recorded_history_3 = defaultdict(list)

# ------------ RECORD HISTORY FOR EACH ROD ----------------
pendulum_sim.collect_diagnostics(rod_1).using(
    PendulumCallBack, step_skip=step_skip, callback_params=recorded_history_1
)
pendulum_sim.collect_diagnostics(rod_2).using(
    PendulumCallBack, step_skip=step_skip, callback_params=recorded_history_2
)
pendulum_sim.collect_diagnostics(rod_3).using(
    PendulumCallBack, step_skip=step_skip, callback_params=recorded_history_3
)
pendulum_sim.finalize()

# ---------- TIME STEPPING SCHEME: POSITION-VERLET --------------
timestepper = PositionVerlet()

# -------------- INTEGRATING THE PDE ---------------
# (Numerical Solver like ode45 in MATLAB or odeint in Python BUT FOR PDEs)
integrate(timestepper, pendulum_sim, final_time, total_steps)


# ------------- SAVE HISTORY IN .dat files ----------
if SAVE_RESULTS:
    import pickle as pickle

    file_name = input("Enter filename: ")
    filename = file_name + "_rod_1.dat"
    with open(filename, "wb") as file:
        pickle.dump(recorded_history_1, file)
    filename = file_name + "_rod_2.dat"
    with open(filename, "wb") as file:
        pickle.dump(recorded_history_2, file)
    filename = file_name + "_rod_3.dat"
    with open(filename, "wb") as file:
        pickle.dump(recorded_history_3, file)