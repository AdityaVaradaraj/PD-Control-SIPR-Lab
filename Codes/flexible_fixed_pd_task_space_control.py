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
# --------- TASK SPACE PD ------------
def TaskFeedBack(rod_number):
    pos_1 = rod_1.position_collection.copy()
    vel_1 = rod_1.velocity_collection.copy()
    tang_1 = rod_1.tangents.copy()
    s_1 = rod_1.lengths.copy()
    curv_1 = rod_1.kappa.copy()
    ang_vel_1 = rod_1.omega_collection.copy()
    q_0_1 = s_1[:-1].dot(curv_1[1,:].T)
    # print(np.shape(pos_1))
    # print(q_0)
    theta_1 = np.arctan2(-1.0*tang_1[2, 0], tang_1[0, 0])
    if(theta_1 < -3.055): # 5 degrees threshold
        theta_1 = -theta_1
    theta_1_dot = ang_vel_1[1, 0]
    # print(theta_dot)
    q_0_1_dot = ang_vel_1[1, -1] - ang_vel_1[1, 0]

    pos_2 = rod_2.position_collection.copy()
    vel_2 = rod_2.velocity_collection.copy()
    tang_2 = rod_2.tangents.copy()
    s_2 = rod_2.lengths.copy()
    curv_2 = rod_2.kappa.copy()
    ang_vel_2 = rod_2.omega_collection.copy()
    q_0_2 = s_2[:-1].dot(curv_2[1,:].T)
    # print(q_0)
    theta_2 = np.arctan2(-1.0*tang_2[2, 0], tang_2[0, 0])
    if(theta_2 < -3.055): # 5 degrees threshold
        theta_2 = -theta_2
    theta_2_dot = ang_vel_2[1, 0]
    # print(theta_dot)
    q_0_2_dot = ang_vel_2[1, -1] - ang_vel_2[1, 0]

    pos_3 = rod_3.position_collection.copy()
    vel_3 = rod_3.velocity_collection.copy()
    tang_3 = rod_3.tangents.copy()
    s_3 = rod_3.lengths.copy()
    curv_3 = rod_3.kappa.copy()
    ang_vel_3 = rod_3.omega_collection.copy()
    q_0_3 = s_3[:-1].dot(curv_3[1,:].T)
    # print(q_0)
    theta_3 = np.arctan2(-1.0*tang_3[2, 0], tang_3[0, 0])
    if(theta_3 < -3.055): # 5 degrees threshold
        theta_3 = -theta_3
    theta_3_dot = ang_vel_3[1, 0]
    # print(theta_dot)
    q_0_3_dot = ang_vel_3[1, -1] - ang_vel_3[1, 0]
    # Jacobian
    if abs(q_0_1) < 0.035:
        if q_0_1_dot >= 0:
            q_0_1 = 0.035 # 2 degrees
        else:
            q_0_1 = -0.035
    
    if abs(q_0_2) < 0.035:
        if q_0_2_dot >= 0:
            q_0_2 = 0.035
        else:
            q_0_2 = -0.035
    
    if abs(q_0_3) < 0.035:
        if q_0_3_dot >= 0:
            q_0_3 = 0.035
        else:
            q_0_3 = -0.035
    
    if q_0_1 != 0 and q_0_2 !=0 and q_0_3 != 0 and q_0_1 is not np.nan and q_0_2 is not np.nan and q_0_3 is not np.nan:
        J_s = np.array([[(base_length*math.cos(q_0_1/2)*math.cos(q_0_1/2 - np.pi/2))/q_0_1 - (2*base_length*math.cos(q_0_1/2 - np.pi/2)*math.sin(q_0_1/2))/q_0_1**2 - (base_length*math.sin(q_0_1/2)*math.sin(q_0_1/2 - np.pi/2))/q_0_1 - (2*base_length*math.sin(q_0_2/2)*math.sin(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2 - (2*base_length*math.sin(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3,
            (base_length*math.cos(q_0_2/2)*math.cos(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2 - (2*base_length*math.sin(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3 - (2*base_length*math.sin(q_0_2/2)*math.cos(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2**2 - (base_length*math.sin(q_0_2/2)*math.sin(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2,  
            (base_length*math.cos(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.cos(q_0_3/2))/q_0_3 - (base_length*math.sin(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3 - (2*base_length*math.cos(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3**2],
            [(2*base_length*math.sin(q_0_1/2)*math.sin(q_0_1/2 - np.pi/2))/q_0_1**2 - (base_length*math.cos(q_0_1/2 - np.pi/2)*math.sin(q_0_1/2))/q_0_1 - (base_length*math.cos(q_0_1/2)*math.sin(q_0_1/2 - np.pi/2))/q_0_1 - (2*base_length*math.sin(q_0_2/2)*math.cos(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2 - (2*base_length*math.cos(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3, 
            (2*base_length*math.sin(q_0_2/2)*math.sin(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2**2 - (base_length*math.sin(q_0_2/2)*math.cos(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2 - (base_length*math.cos(q_0_2/2)*math.sin(q_0_1 + q_0_2/2 - np.pi/2))/q_0_2 - (2*base_length*math.cos(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3, 
            (2*base_length*math.sin(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3**2 - (base_length*math.cos(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2)*math.sin(q_0_3/2))/q_0_3 - (base_length*math.cos(q_0_3/2)*math.sin(q_0_1 + q_0_2 + q_0_3/2 - np.pi/2))/q_0_3]])

        # print(q_0_1, q_0_2, q_0_3)
        # print(J_s)
        X = np.array([[pos_3[2, -1]],[pos_3[0, -1]]])
        X_des = np.array([[0.5*base_length],[-2.75*base_length]])

        X_dot = np.array([[vel_3[2, -1]],[vel_3[0, -1]]])
        X_des_dot = np.array([[0.0],[0.0]]) 
         
        K_p = 3.3275
        K_d = 0.075
        k = 0.942 # Stiffness (N-m/rad)
        torque = np.dot(np.transpose(J_s), (-K_p*(X - X_des) - K_d* (X_dot-X_des_dot)))
        torque += k*np.array([[q_0_1],[q_0_2],[q_0_3]])
        
        # print(np.shape(torque))
    else:
        torque = 0.0*np.ones((3,1))
    return torque[rod_number-1]


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
        system.external_torques[..., 0] += TaskFeedBack(self.rod_number)*self.torque

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
        system.external_torques[..., -1] += TaskFeedBack(self.rod_number)*self.torque

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
