from cProfile import label
import pickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import math
from params import PLOT_VIDEO, PLOT_POSITIONS

file_name = input("Enter filename: ")
file = open(file_name + "_rod_1.dat", "rb")
recorded_history_1 = pickle.load(file)
file.close()
file = open(file_name + "_rod_2.dat", "rb")
recorded_history_2 = pickle.load(file)
file.close()
file = open(file_name + "_rod_3.dat", "rb")
recorded_history_3 = pickle.load(file)
file.close()
def plot_bending_angle(plot_params: dict, rod_number: int):
    curvature = np.array(plot_params["curvature"])
    time = np.array(plot_params["time"])
    lengths = np.array(plot_params["lengths"])
    print(lengths.shape[0], time.shape)
    q = np.zeros((lengths.shape[0],1))
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(50, 8), frameon=True, dpi=150)
    ax = fig.add_subplot(111)
    for i in range(lengths.shape[0]):
        q[i] = lengths[i, :-1].dot(curvature[i,1,:].T)
    if rod_number == 1:
        q_des = [np.pi/8*math.sin(2*np.pi/2.5*t) for t in time]
    elif rod_number == 2:
        q_des = [np.pi/4*math.cos(2*np.pi/5*t) for t in time]
    else:
        q_des = [np.pi/2*math.sin(2*np.pi/4*t + np.pi/9) for t in time]
    ax.plot(time, q,  "k", label = "q_" + str(rod_number))
    ax.plot(time, q_des, "r", label = "q_des")
    ax.set_ylabel("q_" + str(rod_number), fontsize=16)
    ax.set_xlabel("t", fontsize=16)
    ax.set_xlim(0, time[-1])
    # ax.set_ylim(-2*np.pi, 2*np.pi)
    # fig.legend(prop={"size": 16})
    ax.legend(loc='lower right')
    plt.show()

    plt.close(plt.gcf())
plot_bending_angle(recorded_history_1, 1)
plot_bending_angle(recorded_history_2, 2)
plot_bending_angle(recorded_history_3, 3)

def plot_position_end(recorded_history_1, recorded_history_2, recorded_history_3):
    # fig = plt.figure(figsize=(10, 8), frameon=True, dpi=150)
    # ax = fig.add_subplot(111)
    # ax.set_aspect("equal", adjustable="box")
    # Should give a (n_time, 3, n_elem) array
    time = np.array(recorded_history_3["time"])
    print(np.shape(time))
    positions = np.array(recorded_history_1["position"])
    positions_2 = np.array(recorded_history_2["position"])
    positions_3 = np.array(recorded_history_3["position"])
    print(positions_3[:, 2, -1])
    plt.plot(time, positions_3[:, 2, -1], 'k', label='Z')
    plt.plot(time, 0.5*0.13*np.ones(np.shape(time)), 'r', label='Z_des')
    plt.ylabel('Z')
    plt.xlabel('t')
    # fig.show()
    # ax1 = fig.add_subplot(111)
    # ax1.plot(time,  positions_3[:, 0, -1])
    # fig.show()
    plt.legend()
    plt.show()
    plt.plot(time, positions_3[:, 0, -1], 'k', label='X')
    plt.plot(time, -2.75*0.13*np.ones(np.shape(time)), 'r', label='X_des')
    plt.ylabel('X')
    plt.xlabel('t')
    plt.legend()
    plt.show()
    # plt.close(plt.gcf())

plot_position_end(recorded_history_1, recorded_history_2, recorded_history_3)



if PLOT_VIDEO:
    def plot_video(
        plot_params_1: dict,
        plot_params_2: dict,
        plot_params_3: dict,
        video_name="video.mp4",
        margin=0.2,
        fps=60,
        step=1,
        *args,
        **kwargs
    ):  # (time step, x/y/z, node)
        import matplotlib.animation as manimation

        plt.rcParams.update({"font.size": 22})

        # Should give a (n_time, 3, n_elem) array
        positions_1 = np.array(plot_params_1["position"])
        positions_2 = np.array(plot_params_2["position"])
        positions_3 = np.array(plot_params_3["position"])

        print("plot video")
        FFMpegWriter = manimation.writers["ffmpeg"]
        metadata = dict(
            title="Movie Test", artist="Matplotlib", comment="Movie support!"
        )
        writer = FFMpegWriter(fps=fps, metadata=metadata)
        dpi = 300
        fig = plt.figure(figsize=(10, 8), frameon=True, dpi=dpi)
        ax = fig.add_subplot(111)
        ax.set_aspect("equal", adjustable="box")
        # plt.axis("square")
        i = 0
        (rod_line_1,) = ax.plot(positions_1[i, 2], positions_1[i, 0], lw=3.0)
        (tip_line_1,) = ax.plot(positions_1[:i, 2, -1], positions_1[:i, 0, -1], "k--")
        (rod_line_2,) = ax.plot(positions_2[i, 2], positions_2[i, 0], lw=3.0)
        (tip_line_2,) = ax.plot(positions_2[:i, 2, -1], positions_2[:i, 0, -1], "k--")
        (rod_line_3,) = ax.plot(positions_3[i, 2], positions_3[i, 0], lw=3.0)
        (tip_line_3,) = ax.plot(positions_3[:i, 2, -1], positions_3[:i, 0, -1], "k--")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim([-0.1 - margin, 0.1 + margin])
        ax.set_ylim([-0.1 - margin, 0.0 + margin])
        with writer.saving(fig, video_name, dpi):
            with plt.style.context("seaborn-white"):
                for i in range(0, positions_1.shape[0], int(step)):
                    rod_line_1.set_xdata(positions_1[i, 2])
                    rod_line_1.set_ydata(positions_1[i, 0])
                    tip_line_1.set_xdata(positions_1[:i, 2, -1])
                    tip_line_1.set_ydata(positions_1[:i, 0, -1])
                    rod_line_2.set_xdata(positions_2[i, 2])
                    rod_line_2.set_ydata(positions_2[i, 0])
                    tip_line_2.set_xdata(positions_2[:i, 2, -1])
                    tip_line_2.set_ydata(positions_2[:i, 0, -1])
                    rod_line_3.set_xdata(positions_3[i, 2])
                    rod_line_3.set_ydata(positions_3[i, 0])
                    tip_line_3.set_xdata(positions_3[:i, 2, -1])
                    tip_line_3.set_ydata(positions_3[:i, 0, -1])
                    writer.grab_frame()

    plot_video(recorded_history_1, recorded_history_2, recorded_history_3, "swinging_flexible_pendulum_fixed.mp4")

if PLOT_POSITIONS:
    fig = plt.figure(figsize=(10, 8), frameon=True, dpi=150)
    ax = fig.add_subplot(111)
    ax.set_aspect("equal", adjustable="box")
    # Should give a (n_time, 3, n_elem) array
    positions_1 = np.array(recorded_history_1["position"])
    positions_2 = np.array(recorded_history_2["position"])
    positions_3 = np.array(recorded_history_3["position"])
    for i in range(positions_1.shape[0]):
        ax.plot(positions_1[i, 2], positions_1[i, 0], lw=2.0)
        ax.plot(positions_2[i, 2], positions_2[i, 0], lw=2.0)
        ax.plot(positions_3[i, 2], positions_3[i, 0], lw=2.0)
    fig.show()
    plt.show()